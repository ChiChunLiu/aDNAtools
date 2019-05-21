#!/bin/python
import os
import argparse
from numpy import where
from pysam import VariantFile, FastaFile
from utils import allele_combination, any_match, base_complement


vcf2eigen_geno_dic = {(0, 0): '2', (0, 1): '1', (1, 0): '1', (1, 1): '0', (-1, -1): '9', (None, None): '9'}

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


if os.path.exists(args.output + '.geno'):
    raise Exception('outputs have already exsited.')

vcf = VariantFile(args.vcf)
ref_fa = FastaFile(args.reference)
    
nsample = len((vcf.header.samples))


with open(args.target, 'r') as target, open(args.output + '.snp', 'a') as snp, open(args.output + '.geno', 'a') as geno:
    
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
                    vid = (rec.chrom + ':' + str(rec.pos)) if rec.id is None else rec.id
                    snp.write('\t'.join([vid, rec.chrom, '0', str(rec.pos), rec.ref, rec.alts[0]]) + '\n')
                    geno.write(''.join(gt) + '\n')
