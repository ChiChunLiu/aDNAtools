import pandas as pd
import numpy as np
import argparse

'''
interpolate genetic position into an eigenstrat snp file
'''

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--map', type = str, default = "", help = "input prefix")
parser.add_argument('-i', '--input', type = str, default = "", help = "input snp file")
parser.add_argument('-o', '--output', type = str, default = "", help = "output snp file")
args = parser.parse_args()

recomb_map = pd.read_csv(args.map,
                         sep = '\t',
                         header = None,
                         names = ['CHROM', 'POS', 'CM'])


snp_input = pd.read_csv(args.input,
                        sep = '\t',
                        header = None,
                        names = ['ID', 'CHROM', 'CM', 'POS', 'REF', 'ALT'])


def interpolate_genetic_position(recomb_map, chrom, pos):
    '''return interpolated genetic position for a single chromosome
    '''
    recomb_map = recomb_map[recomb_map['CHROM'] == chrom]
    y = np.interp(x = pos, 
                  xp = recomb_map.POS,
                  fp = recomb_map.CM,
                  left = 0,
                  right = np.max(recomb_map.CM))
    return list(y)


all_pos = []
for c in range(1, 23):
    pos = snp_input[snp_input['CHROM'] == c]
    y = interpolate_genetic_position(recomb_map, chrom = c, pos = pos)
    all_pos.extend(y)

snp_output = snp_input.copy()
snp_output['CM'] = all_pos
snp_output[['ID', 'CHROM', 'CM', 'POS', 'REF', 'ALT']].to_csv(
    args.output,
    sep = '\t',
    index= False, 
    header = False)
    
    