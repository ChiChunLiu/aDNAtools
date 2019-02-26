#!/bin/python3

import numpy as np
import argparse
import gzip


base_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

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


def next_target_snp(snps, skip_ambiguous_snp = True):
        s = next(snps)
        target_id, target_chrom, target_pos, target_a0, target_a1 = s.strip().split('\t')
        
        if skip_ambiguous_snp:
            while base_complement[target_a0] == target_a1:
                s = next(snps)
                target_id, target_chrom, target_pos, target_a0, target_a1 = s.strip().split('\t')

        return target_id, target_chrom, target_pos, target_a0, target_a1


def random_draw(pileup, target, output_file, remove_ambiguous_snp = True):

    open_vcf = vcf_open(pileup)
    samples = get_samples(pileup)
    n = len(samples)
    
    with open(output_file + '.geno', 'a') as geno, open(output_file + '.snp', 'a') as snp, open(output_file + '.ind', 'a') as ind:
        with open_vcf(pileup, 'r') as vcf, open(target, 'r') as snps:
            
            # first target snp
            target_id, target_chr, target_pos, target_a0, target_a1 = next_target_snp(snps)

            for line in vcf:
                if not line.startswith(b'#'):
                    entry = line.decode('utf-8').strip().split('\t')
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = entry[:8]

                    if (CHROM == target_chr) & (POS == target_pos):
                        FORMAT = entry[-n-1]
                        records = entry[-n:]

                        GT_pseudo = []

                        match_allele = {target_a0: target_a1,
                                        target_a1: target_a0,
                                        base_complement[target_a0]: base_complement[target_a1],
                                        base_complement[target_a1]: base_complement[target_a0]}

                        matched_alt = match_allele[REF]

                        # locate alternative allele in the target file
                        if matched_alt in ALT.split(','):
                            matched_alt_index = ALT.split(',').index(matched_alt)
                        else:
                            matched_alt_index = 0

                        # check allele swap
                        if (REF == target_a0) | (REF == base_complement[target_a0]):
                            swap = False
                        else:
                            swap = True

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
                        snp.write('\t'.join([target_id, target_chr, '0', target_pos, target_a0, target_a1]) + '\n')
                        # next target snp
                        target_id, target_chr, target_pos, target_a0, target_a1 = next_target_snp(snps)
                        
                    else:
                        GT_pseudo = np.repeat('9', n)
                        geno.write(''.join(GT_pseudo) + '\n')
                        snp.write('\t'.join([target_id, target_chr, '0', target_pos, target_a0, target_a1]) + '\n')
                        

        for s in samples:
            ind.write('\t'.join([s, 'U', 'pop']) + '\n')

            
if __name__== "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--snp', type = str, default = "", help = "biallelic SNP table")
    parser.add_argument('-p', '--pileup', type = str, default = "", help = "bcftools mpileup")
    parser.add_argument('-o', '--output', type = str, default = "", help = "output")   
    args = parser.parse_args()
        
    random_draw(args.pileup, args.snp, args.output)

