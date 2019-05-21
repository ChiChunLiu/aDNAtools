import argparse
import numpy as np
from scipy.sparse import coo_matrix

'''
print fraction of missing entries per individual and per snp
'''

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type = str, default = "", help = "input prefix")
parser.add_argument('-o', '--output', type = str, default = "", help = "output prefix")
args = parser.parse_args()

with open(args.input + '.geno') as f:
    p = sum(1 for line in f if line.strip())
    
with open(args.input + '.geno') as f:
    f.seek(0)
    first_line = f.readline().strip()
    n = len(list(first_line))
    
def find_missing(geno_file, shape):

    row = []
    col = []
    data = []

    with open(geno_file, 'r') as geno:
        for p, gline in enumerate(geno):
            gline = list(gline.strip())
            c_idx = [i for i,g in enumerate(gline) if g == '9']
            if c_idx:
                r_idx = list(np.repeat(p, len(c_idx)))
                d = list(np.repeat(1, len(c_idx)))s
                row.extend(r_idx)
                col.extend(c_idx)
                data.extend(d)
    return coo_matrix((data, (row, col)), shape = shape)


M = find_missing(args.input + '.geno', (p, n))
ind_miss = np.asarray(M.sum(axis = 0))[0] / p
snp_miss = np.asarray(M.sum(axis = 1))[:,0] / n

np.savetxt(args.output + '.snp.miss', snp_miss, fmt='%.4f')
np.savetxt(args.output + '.ind.miss', ind_miss, fmt='%.4f')
