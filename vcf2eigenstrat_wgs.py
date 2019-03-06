#!/bin/python
from numpy import where
from pysam import VariantFile, FastaFile

base_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
vcf2eigen_geno_dic = {(0, 0): '2', (0, 1): '1', (1, 0): '1', (1, 1): '0', (-1, -1): '9'}

def allele_combination(alleles):
    '''
    argument:
    ---------
    (a0, a1)
    
    return:
    -------
    alleles combitation: list
    [(a0, a1), (a0_c, a1_c), (a1, a0), (a1_c, a0_c)]
    '''
    alleles_c = tuple([base_complement[a] for a in alleles])
    return [alleles, alleles_c, alleles [::-1], alleles_c[::-1]]

def any_match(a0, a_list):
    '''
    compare a0 to a list of alleles
    if there is matching pairs in a0 == a_list[i]
    
    e.g 
    list(zip(('A','0'), ('A','C') )) 
            [('A', 'A'), ('0', 'C')] is TRUE
    '''
    match = []
    for a in a_list:
        match.append(any(x == y for x, y in zip(a0, a)))
    return match

def strip_target(target):
    target = target.strip().split()
    target_chrom, target_pos = target[1], int(target[3])
    target_alleles = (target[4], target[5])
    allele_cb = allele_combination(target_alleles)
    
    return target_chrom, target_pos, target_alleles, allele_cb

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', type = str, default = "", help = "vcf input")
parser.add_argument('-f', '--reference', type = str, default = "", help = "fasta reference file")
parser.add_argument('-t', '--target', type = str, default = "", help = "eigenstrat target snp")
parser.add_argument('-o', '--output', type = str, default = "", help = "output prefix")
args = parser.parse_args()


if os.path.exists(args.out + '.geno'):
    raise Exception('outputs have already exsited.')

vcf = VariantFile(args.vcf)
ref_fa = FastaFile(args.reference)
    
nsample = len((vcf.header.samples))


with open(args.target, 'r') as target, open(args.out + '.snp', 'a') as snp, open(args.out + '.geno', 'a') as geno:
    
    for t in target:
        target_chrom, target_pos, target_alleles, allele_cb = strip_target(t)
       
        for rec in vcf.fetch(target_chrom, target_pos -1, target_pos):
            # record exisits in vcf
            if not rec.chrom:
                fasta_ref = ref_fa.fetch(target_chrom, target_pos -1, target_pos)
                if fasta_ref in list(base_complementent.keys()):
                    # find matched alt allele
                    matched_alleles = allele_cb[where(any_match((fasta_ref, '0'), allele_cb))[0][0]]
                    snp.write('\t'.join(['.', target_chrom, '0', str(target_pos), fasta_ref, matched_alleles[1]]) + '\n')
                    geno.write('2' * nsample + '\n')
            else:
                # skip if no matched alleles
                alleles = (rec.ref, rec.alts[0])
                if alleles in allele_cb:
                    gt = [vcf2eigen_geno_dic[s['GT']] for s in rec.samples.values()]
                    # deal with no variant ID
                    vid = (rec.chrom + str(rec.pos)) if rec.id is None else rec.id
                    snp.write('\t'.join([vid, rec.chrom, '0', str(rec.pos), rec.ref, rec.alts[0]]) + '\n')
                    geno.write(''.join(gt) + '\n')
