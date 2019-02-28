#!/bin/python3

import numpy as np
import argparse
import sys
import gzip

'''
Pull pileup file down into pseudo-haploid data in eigenstrat format.
Bcftools mpileup takes a while (~1hr for array data), and pull-down
takes only a few minutes. The script streams through both vcf and 
target SNP files, so it's memory efficient.

Usage:  

bcftools mpileup --ignore-RG -B -q30 -Q30 \
-R file_name.pos \
-f hs37d5.fa \
-b bam.txt \
-a AD,DP | \
bcftools norm -f hs37d5.fa -Oz -o file_name.mpileup.vcf.gz

python random_draw.py \
-s some_path/file_name.pos \
-p some_path/file_name.mpileup.vcf.gz \
-o dest/file_name.ancient

REMEMBER:
 - B flag disables probabilistic realignment for the computation
   of base alignment quality (BAQ). 
'''

base_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

def progressBar(bar_length, value, total):
    '''
    Prints a simple progress bar denoting progress in a for loop
    From Alan. Alan did not do the weird face
    '''

    frac = value/total
    prog = 'ʕ·ᴥ·ʔ ' * int(round(frac * bar_length))
    space = '      ' * (bar_length - len(prog)//6)

    sys.stdout.write('\rProgress: [{0}] {1:2.2f}%'.format('ʕ·ᴥ·ʔ ' + prog + space, frac*100))
    sys.stdout.flush()

    if frac >= 1:
        print('\n')
    
def sum_line(file):
    with open(file) as f:
        return sum(1 for _ in f)

def vcf_open(filename):
    if filename.endswith('.gz'):
        open_vcf = gzip.open
    else:
        open_vcf = open
    return open_vcf


def get_samples(filename):
    open_vcf = vcf_open(filename)
    
    with open_vcf(filename, 'r') as vcf:
        for line in vcf:
            if line.startswith(b'#CHROM'):
                samples_str = line.decode('utf-8').split('FORMAT')[1]
                samples = samples_str.strip().split('\t')
                
                return samples


def next_vcf_line(vcf, skip_ambiguous_snp = True):
    try:
        line = next(vcf)
        while line.startswith(b'#'):
            line = next(vcf)

        entry = line.decode('utf-8').strip().split('\t')
        REF, ALT= entry[3:5]

        if skip_ambiguous_snp:
            while base_complement[REF] == ALT:
                
                line = next(vcf)
                entry = line.decode('utf-8').strip().split('\t')
                REF, ALT= entry[3:5]
        return entry
    except:
        return 'EOF'


def random_draw(pileup, target, output_file, progress = True, remove_ambiguous_snp = True):
    '''
    Draw a random allele representation from reads at target position,
    and write it as eigenstrat format.
    
    arguments
    ---------
    pileup: str
        bcftools mpileup with vcf or vcf.gz format
    target: str
        target SNP file with 5 columns - ID, CHROM, POS, A0, A1
    output_file: str
        prefix for ouput file        
    '''
    open_vcf = vcf_open(pileup)
    samples = get_samples(pileup)
    n = len(samples)
    nline_target = sum_line(target)

    with open(output_file + '.geno', 'a') as geno, open(output_file + '.snp', 'a') as snp, open(output_file + '.ind', 'a') as ind:
        with open_vcf(pileup, 'r') as vcf, open(target, 'r') as snps:
            
            entry = next_vcf_line(vcf)
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = entry[:8]
            CHROM, POS = int(CHROM), int(POS)
            
            snp_counter = 0
            for s in snps:
                
                # progress bar
                if progress:
                    snp_counter += 1
                    if snp_counter%1000 == 0:
                        progressBar(20, snp_counter, nline_target)
                        
                target_id, target_chr, target_pos, target_a0, target_a1 = s.strip().split('\t')
                # check ambiguos strand: TRUE -> next target SNP
                if remove_ambiguous_snp:
                    if base_complement[target_a0] == target_a1:
                        continue
                target_chr, target_pos = int(target_chr), int(target_pos)
                
                if (CHROM == target_chr) & (POS == target_pos):
                    FORMAT = entry[-n-1]
                    records = entry[-n:]

                    match_allele = {target_a0: target_a1,
                                    target_a1: target_a0,
                                    base_complement[target_a0]: base_complement[target_a1],
                                    base_complement[target_a1]: base_complement[target_a0]}
                    matched_alt = match_allele[REF]

                    # locate alternative allele in the target file
                    # add one for REF
                    if matched_alt in ALT.split(','):
                        matched_alt_index = ALT.split(',').index(matched_alt) + 1
                    else:
                        matched_alt_index = 1

                    # check allele swap
                    if (REF == target_a0) | (REF == base_complement[target_a0]):
                        swap = False
                    else:
                        swap = True
                    
                    GT_pseudo = []
                    for r in records:
                        # parse record using AD location in FORMAT
                        AD = r.split(':')[FORMAT.split(':').index('AD')]
                        AD = np.array(AD.split(',')).astype(int)

                        ADsum = AD[0] + AD[matched_alt_index]
                        if ADsum > 0:
                            rdraw = 2 - 2 * np.random.binomial(1, AD[matched_alt_index]/ADsum)
                            if swap:
                                rdraw = 2 - rdraw
                            GT_pseudo.append(str(rdraw))                        
                        else:
                            GT_pseudo.append('9')

                    geno.write(''.join(GT_pseudo) + '\n')
                    snp.write('\t'.join([target_id, str(target_chr), '0', str(target_pos), target_a0, target_a1]) + '\n')
                    
                    # next vcf line
                    entry = next_vcf_line(vcf)
                    if entry == 'EOF':
                        continue
                    else:
                        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = entry[:8]
                        CHROM, POS = int(CHROM), int(POS)

                # mpileup position surpasses the target position due to no coverage across samples
                elif (CHROM > target_chr) | ((CHROM == target_chr) & (POS > target_pos)):
                    GT_pseudo = np.repeat('9', n)
                    geno.write(''.join(GT_pseudo) + '\n')
                    snp.write('\t'.join([target_id, str(target_chr), '0', str(target_pos), target_a0, target_a1]) + '\n')
                
                # mpileup position is before target position. Should never happen except at the end of target files
                else: 
                    geno.write(''.join(GT_pseudo) + '\n')
                    snp.write('\t'.join([target_id, str(target_chr), '0', str(target_pos), target_a0, target_a1]) + '\n')


        for s in samples:
            ind.write('\t'.join([s, 'U', 'pop']) + '\n')
            
if __name__== "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--snp', type = str, default = "", help = "biallelic SNP table")
    parser.add_argument('-p', '--pileup', type = str, default = "", help = "bcftools mpileup")
    parser.add_argument('-o', '--output', type = str, default = "", help = "output")   
    args = parser.parse_args()
        
    random_draw(args.pileup, args.snp, args.output)

