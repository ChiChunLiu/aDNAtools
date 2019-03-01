import os
import pandas as pd
import numpy as np
import itertools as it

'''

Eigenstrat I/O utilities

'''

class eigenstrat(object):
    """
    """
    def __init__(self, geno=None, snp=None, ind=None):

        # p x n genotype matrix
        self.geno = geno

        # DataFrame of SNP level information
        self.snp = snp
        
        # list of individual ids
        self.ind = ind

    def __geno2string(self, x):
        
        #source: https://stackoverflow.com/questions/2721521/
        #fastest-way-to-generate-delimited-string-from-1d-numpy-array/13861407
        #generate an array with strings
        x_arrstr = np.char.mod('%i', x)
        #combine to a string
        x_str = ''.join(x_arrstr)
        return x_str

    def __write_eigenstrat_geno(self, path):
        '''
        argument:
        ----------
        geno: ndarray
            genotype numpy array
        path: string
            file path
        '''
        if os.path.exists(path):
            raise Exception('file already exists!')
        # To Do: check 0/1/2/9
        for p in range(self.geno.shape[0]):
            with open(path, 'a') as f:
                f.write(self.__geno2string(self.geno[p,:]) + '\n')

    def write_eigenstrat(self, prefix):
        '''
        write prefix.geno, prefix.snp, prefix.ind

        argument:
        ---------
        geno: ndarray
            genotype numpy array
        prefix: string
            prefix for the outputs
        '''
        self.__write_eigenstrat_geno(path = prefix + '.geno')

        if list(self.snp.columns.values) != ['id', 'chr', 'gpos', 'pos', 'a0', 'a1']:
            raise Exception('please reformat your snp file')
        self.snp.to_csv(prefix + '.snp', sep='\t', index = False, header = False)

        if list(self.ind.columns.values) != ['sample', 'sex', 'population']:
            raise Exception('please reformat your ind file')
        self.ind.to_csv(prefix + '.ind', sep='\t', index = False, header = False)

        
'''

miscellaneous utility function

'''

def numpy2allel(x):
    '''
    convert 2d numpy array into scikit allele
    compatible 3d numpy array
    2 -> [1, 1]
    0 -> [0, 0]
    1 -> [0, 1]
    '''
    x1 = np.copy(x)
    x1[x1 == 1] = 0
    x1[x1 == 2] = 1
    
    x2 = np.copy(x)
    x2[x2 == 2] = 1
    
    return np.stack((x1,x2), axis = -1)

'''tokenization
'''

chrom_digit = {'X':23, 'Y':24, 'chrM': 90, 'MT':90, 'M': 90, 'XY':91, '0': 0}
for c in list(range(1,23)):
    chrom_digit[str(c)] = int(c)

sex_token = {'Female': 'F', 'Male': 'M',
             '1': 'M', '2': 'F'}

def digitize_chrom(chrom, chrom_digitize_dict = chrom_digit):
    '''mapping chromosome from {'1'..'22','X','Y','MT'} to {1~24,90}
    '''
    return chrom_digitize_dict[chrom]

def tokenize_sex(sex, token_dict = sex_token):
    if sex in sex_token.keys():
        return sex_token[sex]
    else:
        return 'U'
    
def space2underscore(x):
    '''
    return
    ------
        string: str
        input with underscores inplace of spaces
    '''
    return x.replace(" ", "_")
