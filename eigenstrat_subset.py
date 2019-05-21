import argparse
import random
import pandas as pd


def read_population_draw_table(file):
    '''specify population and n. e.g UpperMustang 30
    '''
    population_map = {}
    with open(file, 'r') as f:
        for line in f:
            p, n = line.strip().split()
            if p in population_map.keys():
                raise ValueError('sample population subsetted multiple times.')
            else: 
                population_map[p] = n
    return population_map

def read_kept_individual(file):
    '''spefify a list of individual wanting to keep 
    '''
    with open(file, 'r') as f:
        inds = [i.strip() for i in f] 
    return inds


def read_kept_population(file):
    '''spefify a list of population wanting to keep
    '''
    with open(file, 'r') as f:
        pops = [p.strip() for p in f] 
    return pops

def population_draw(clst_df, subset_dict):
    '''draw random samples according to the dictionary
    '''
    pops = list(subset_dict.keys())
    all_dfs = [clst_df[~clst_df['population'].isin(pops)]]
    for p in pops:
        if int(subset_dict[p]) > 0:
            df = clst_df[clst_df['population'] == p]
            df = df.sample(n = int(subset_dict[p]))
            all_dfs.append(df)    

    df = pd.concat(all_dfs)
    df = df.sort_index()
    return df


if __name__== "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type = str, default = "", help = "prefix for inputs")
    parser.add_argument('-r', '--seed', type = str, default = "", help = "random seed")
    parser.add_argument('-T', '--table', type = str, default = "", help = "population subset table. e.g. UpperMustang 30")
    parser.add_argument('-S', '--sample', type = str, default = "", help = "samples to keep")
    parser.add_argument('-P', '--population', type = str, default = "", help = "populations to keep")
    parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for outputs")
    args = parser.parse_args()
    
    
    sub = [args.table, args.sample, args.population]

    clst = pd.read_csv(args.input + '.ind', sep = '\s+', header = None)
    if clst.shape[1] != 3:
        clst = pd.read_csv(args.input + '.ind', sep = '\t', header = None)
    clst.columns = ['sample', 'sex', 'population']
        
    if not ''.join(sub):
        raise ValueError('subset input error')
    
    elif args.table:
        random.seed(args.seed)
        population_map = read_population_draw_table(args.table)
        clst_subset = population_draw(clst, population_map)
        
    elif args.sample:
        inds = read_kept_individual(args.sample)
        clst_subset = clst[clst['sample'].isin(inds)]
        
    else:
        pops = read_kept_population(args.population)
        clst_subset = clst[clst['population'].isin(pops)]
        
    clst_subset.to_csv(args.output + '.ind', sep = '\t', index = False, header = False)
    keep_index = clst_subset.index
    
    with open(args.input + '.geno', 'r') as geno, open(args.output + '.geno', 'a') as geno_out:
        for line in geno:
            gt = list(line.strip())
            gt = [gt[i] for i in keep_index]
            gt = ''.join(gt) + '\n'
            geno_out.write(gt)

    with open(args.input + '.snp', 'r') as snp, open(args.output + '.snp', 'a') as snp_out:
        for line in snp:
            snp_out.write(line)

