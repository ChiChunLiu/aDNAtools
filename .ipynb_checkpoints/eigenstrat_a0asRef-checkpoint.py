#!/bin/python
import os
import argparse
from numpy import where
from pysam import FastaFile
from utils import allele_combination, base_complement

'''
The script convert eigenstrat files a0 alllele to reference 
'''
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type = str, default = "", help = "input prefix")
parser.add_argument('-f', '--reference', type = str, default = "", help = "fasta reference file")
parser.add_argument('-o', '--output', type = str, default = "", help = "output prefix")
args = parser.parse_args()

flip_dict = {'0': '2', '2': '0', '1': '1', '9': '9'}
def flip_geno(geno):
    return [flip_dict[g] for g in geno]    

def find_alt_allele(ref, allele_cb):
    match = [ref == a[0] for a in allele_cb]     
    idx = where(match)[0][0]
    alt = allele_cb[idx][1]
    
    return idx, alt


if os.path.exists(args.output + '.geno'):
    raise Exception('Outputs have already exsited.')

ref_fa = FastaFile(args.reference)

with open(args.input + '.snp', 'r') as snp_input, open(args.input + '.geno', 'r') as geno_input, open(args.output + '.snp', 'a') as snp_output, open(args.output + '.geno', 'a') as geno_output:
    for sline, gline in zip(snp_input, geno_input):
        snpID, chrom, cm, pos, a0, a1 = sline.strip().split()
        fasta_ref = ref_fa.fetch(chrom, int(pos) -1, int(pos))
        allele_cb = allele_combination((a0, a1))
        # if not strand ambiguous
        if allele_cb[0] != allele_cb[3]:
            if fasta_ref in list(base_complement.keys()):
                # find matched alt allele
                idx, alt = find_alt_allele(fasta_ref, allele_cb)              
                snp_output.write('\t'.join([snpID, chrom, cm, pos, fasta_ref, alt]) + '\n')
                # do not flip if matched index = 0 or 1
                if idx in [0,1]:
                    geno_output.write(gline)  
                else:
                    geno_output.write(''.join(flip_geno(list(gline.strip()))) + '\n')
