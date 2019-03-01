#!/bin/python

import pandas as pd
import argparse
from numpy import where

'''
Usage

python fixed2a0a1.py -a eigenfile1 -b eigenfile2 -o output_prefix

where eigenfile2 contains fixed positions denoted as ('0', A) 

-h flag specifies the notation for a fixed allele.
'''


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--input1', type = str, default = "", help = "prefix for anchor")
parser.add_argument('-b', '--input2', type = str, default = "", help = "prefix for target")
parser.add_argument('-o', '--output', type = str, default = "", help = "output prefix")
parser.add_argument('-f', '--fixed_allele', type = str, default = "0", help = "fixed allele")
args = parser.parse_args()


hom = args.fixed_allele

col_name = ['id', 'chrom', 'cm', 'pos', 'a0', 'a1']
snp_df_1 = pd.read_csv(args.input1 + '.snp', header = None, sep = '\t', names = col_name).reset_index()
snp_df_2 = pd.read_csv(args.input2 + '.snp', header = None, sep = '\t', names = col_name).reset_index()
snp_df_merged = snp_df_1.merge(snp_df_2, how = 'inner', on = ['chrom', 'pos'], suffixes = ['_df1', '_df2'])


# finding the right allele pair and remove strand ambiguous variants
base_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

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
    ... [('A', 'A'), ('0', 'C')] is TRUE
    '''
    match = []
    for a in a_list:
        match.append(any(x == y for x, y in zip(a0, a)))
    return match

A0, A1 = [],[]

for index, variant in snp_df_merged.iterrows():    
    if (variant['a0_df2'] != hom) & (variant['a1_df2'] != hom):
        A0.append(variant['a0_df2'])
        A1.append(variant['a1_df2'])
    else:
        # check ambiguity form data 1
        allele_cb = allele_combination((variant['a0_df1'], variant['a1_df1']))
        if allele_cb[0] == allele_cb[1][::-1]:
            A0.append('NA')
            A1.append('NA')
        else:
            # find ma
            alleles = allele_cb[where(any_match((variant['a0_df2'], variant['a1_df2']), allele_cb))[0][0]]
            A0.append(alleles[0])
            A1.append(alleles[1])

snp_df_merged.loc[:, 'A0'] = A0
snp_df_merged.loc[:, 'A1'] = A1

# index for output snps
keep_index = list(snp_df_merged.index_df2)
keep_index.sort(reverse = True)

# write a snp file
snp_df_merged = snp_df_merged[['id_df2', 'chrom', 'cm_df2', 'pos', 'A0', 'A1']].copy()
snp_df_merged.to_csv(args.output + '.snp', sep = '\t', index = False, header = False)

# write a genotype file
with open(args.input2 + '.geno', 'r') as geno_in, open(args.output + '.geno', 'a') as geno_out:
    line = 0
    l = keep_index.pop()
    for g in geno_in:
        if line == l:
            geno_out.write(g)
            if len(keep_index) > 0:
                l = keep_index.pop()
        line += 1

