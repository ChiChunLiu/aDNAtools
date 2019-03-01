#!/usr/bin/python3

from pandas_plink import read_plink
import argparse
import pandas as pd
import numpy as np
import utils as utl

'''
Usage

python aDNAtools/plink2eigenstrat.py -b plink_prefix -o eigenstrat_prefix -p / -P population_label
current implementation reads the whole genotype matrix into memory. It takes ~5G for 300 samples on 
the OmniExpress array (1M variants)

'''


def numpy2eigen_geno(geno):
    gt = geno.compute()
    gt = 2 - gt
    gt[np.isnan(gt)] = 9
    return gt
    

def fam2ind(fam):
    sample = fam[['iid', 'gender', 'fid']].copy()
    sample.columns = ['sample', 'sex', 'population']
    sample.loc[:,'sex'] = [utl.tokenize_sex(s) for s in sample['sex'].values]
    return sample


def bim2snp(bim):
    bim = bim[['snp','chrom','cm','pos','a0','a1']].copy()
    bim.loc[:,'chrom'] = bim.chrom.astype(int)
    bim.loc[:,'a0'] = bim.a0.astype(str)
    bim.loc[:,'a1'] = bim.a1.astype(str)
    bim.columns = ['id', 'chr', 'gpos', 'pos', 'a0', 'a1']
    bim = bim.sort_values(['chr', 'pos'])
    return bim

if __name__== "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input', type = str, default = "", help = "prefix for bed/bim/fam")
    parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for outputs")
    parser.add_argument('-p', '--population', type = str, default = "", help = "one population")
    parser.add_argument('-P', '--ind2pop', type = str, default = "", help = "list of populations")
    args = parser.parse_args()

    (bim, fam, geno) = read_plink(args.input)
   
    snp = bim2snp(bim)
    
    gt = numpy2eigen_geno(geno)
    gt = gt[snp.index,:]
    
    sample = fam2ind(fam)
    if (args.population != '') & (args.ind2pop != ''):
        raise Exception('flag -p and -P cannot be used together.')
    elif args.population is not '':
        sample.loc[:,'population'] = args.population
    elif args.ind2pop is not '':
        pops = []
        with open(args.ind2pop) as pfile:
            for p in pfile:
                p = p.strip()
                pops.append(p)
        sample.loc[:,'population'] = pops

    # write files
    out = utl.eigenstrat(geno = gt,
                         snp = snp,
                         ind = sample)
    out.write_eigenstrat(args.output)

