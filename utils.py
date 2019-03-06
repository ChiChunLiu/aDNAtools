import os
import pandas as pd
import numpy as np
import itertools as it

'''
Complementary utilities to ancpy from the NovembreLab
'''

class eigenstrat(object):
    """Eigenstrat I/O utilities
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
    

class vcf(object):
    """
    Simple VCF I/O utilities
    """
    def __init__(self, geno=None, snp_df=None, inds=None):

        # p x n genotype matrix
        self.geno = geno
        # DataFrame of SNP level information
        self.snp_df = snp_df
        # list of individual ids
        self.inds = inds
        self.gt2string = {0: "0/0", 1: "0/1", 2: "1/1", np.nan: "./."}

    def __write_rec(snp_info, gt):
        gt = [self.gt2string[g] for g in gt]
        rec = '\t'.join(map(str, snp_info)) + '\t'
        rec += '\t'.join([".", ".", ".", "GT"]) + '\t'
        rec += '\t'.join(gt)
        vcf_print(rec)

    def __print_header(samples):
        self.__vprint("##fileformat=VCFv4.2", overwrite = True)
        self.__vprint("##Source=ancpy util extension")
        self.__vprint("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        s = '\t'.join(samples)
        self.__vprint("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + s)

    def __vprint(x, overwrite = False):
        if overwrite:
            print(x, file = open(prefix, "w"))
        else:
            print(x, file = open(prefix, "a"))  

    def write_vcf(self, prefix):
        
        self.__print_header(self.inds)
        for index, snp in self.snp_df.iterrows():
            snp_info = snp[["CHROM", "POS", "ID", "A1", "A2"]].tolist()
            snp_info = [str(i) for i in snp_info]
            gt = self.geno[index,:]
            self.__write_rec(snp_info, gt)


#class plink(geno, snp, ind):
# TBD   


        
        
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
            [('A', 'A'), ('0', 'C')] is TRUE
    '''
    match = []
    for a in a_list:
        match.append(any(x == y for x, y in zip(a0, a)))
    return match
