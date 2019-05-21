#!/bin/python
import argparse
import pandas as pd
import numpy as np


def get_population(file):            
    with open(file, 'r') as f:
        populations = f.readlines()
    populations = [x.strip() for x in populations] 
    return populations
    
def get_kept_index(pops, ind):
    pop_count = {}
    kept_index = []
    with open(ind, 'r') as f:
        for index, line in enumerate(f):
            sample, sex , population = line.strip().split()
            
            if population not int pop_count.keys():
                pop_count[population] = 1
            else pop_count[population] += pop_count[population]
            
            if population in pops:
                kept_index.append(index)
    return kept_index

def get_maf(gline, kept_index):
    gt = list(gline.strip())
    gt = np.array([float(g) for g in gt])
    gt = gt[kept_index]
    gt[gt == 9] = np.nan
    af = np.nanmean(gt) * 0.5
    return np.array([af, 1 - af]).min()

if __name__== "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type = str, default = "", help = "prefix for inputs")
    parser.add_argument('-p', '--population', type = str, default = "", help = "population list")
    parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for outputs")
    parser.add_argument('-e', '--threshold', type = float, default = .05, help = "maf threshold")
    args = parser.parse_args()
    
    
    pops = get_population(args.population)
    kept_index = get_kept_index(pops = pops, args.input + '.ind')  
   
    with open(args.input + '.geno', 'r') as geno, open(args.input + '.snp', 'r') as snp, \
         open(args.output + '.geno', 'a') as geno_out, open(args.input + '.snp', 'a') as snp_out:
            for gline, sline in zip(geno, snp):
                maf = get_maf(gline)
                if gline > eps:
                    geno_out.write(gline)
                    snp_out.write(sline)
