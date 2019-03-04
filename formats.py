import sys
import os
import itertools as it
import pandas as pd
import numpy as np
#import pysam

# from ancpy project on the Novembre Lab
# To be removed after the repo stablizes
# https://github.com/NovembreLab/ancpy

class Genotypes(object):
    """A class for storing genotype data which at its core 
    consists of a dense p x n genotype matrix, a dataframe 
    of metadata associated which each snp, and a list of 
    individual ids. We also provide computation options for 
    setting up dataset for commonly using analyses

    Arguments
    ---------
    Y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    Attributes
    ----------
    Y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int 
        number of snps (features)

    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, Y=None, snp_df=None, inds=None, clst=None, population_df=None, build='hg19'):

        # p x n genotype matrix
        self.Y = Y

        # DataFrame of SNP level information
        self.snp_df = snp_df
        
        # list of individual ids
        self.inds = inds
        
        # individual information
        self.clst = clst

        # individual information
        self.population_df = population_df
        
        # genome build
        self.build = build
        
        # below is adapted from 
        # https://github.com/bulik/ldsc/blob/master/munge_sumstats.py

        # complement of each base
        self.base_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

        # bases 
        self.bases = self.base_complement.keys()

        # dict of allele combinations that are strand ambiguous
        self.strand_ambiguous_alleles = {"".join(x): x[0] == self.base_complement[x[1]]
                                         for x in it.product(self.bases, self.bases)
                                         if x[0] != x[1]} 

    def remove_chrom_prefix(self):
        """
        Remove the chromsome prefix from the snp_df
        this helps for making datasets consistent before 
        merging
        """
        self.snp_df["CHROM"] = self.snp_df["CHROM"].apply(lambda x: int(x[3:]))

    def keep_autosomes(self):
        """Keep snps only present on the autosomes
        """
        # filter out snps not on the autosome
        idx = np.array((self.snp_df["CHROM"] <= 22).index)
        
        return(idx)

    def remove_rare_common(self, eps):
        """
        Removes SNPs that are too rare or 
        too common given some threshold eps

        Arguments
        ---------
        eps : float
            frequency threshold eps
        """
        # compute allele frequencies at each SNP
        f = np.nansum(self.Y, axis=1) / (2 * self.n)

        # get the indices of SNPs that dont pass the thresh
        idx = np.where((f >= eps) & (f <= (1. - eps)))[0]
        
        return(idx)
    
    def filter_bad_snps(self, eps):
        """
        Removes SNPs with too many missing
        entries
        
        Arguments
        ---------
        eps : float
            frequency threshold eps
        """
        # compute allele frequencies at each SNP
        f = np.isnan(self.Y).sum(axis = 1) / self.Y.shape[1]

        # get the indices of SNPs that dont pass the thresh
        idx = np.where(f <= eps)[0]

        return(idx)

    def remove_strand_ambiguous_alleles(self):
        """
        Removes any snps that have alleles that
        are strand ambiguous. 

        TODO: speed this up its currently a bit slow
        also add removal of incosistent alleles
        """
        # lambda function for applying to each snp
        strand_unamb_func = lambda row: not self.strand_ambiguous_alleles[row["A1"] + row["A2"]]

        # snps that have unambigous strand
        idx = self.snp_df.apply(strand_unamb_func, axis=1)

        return(idx)
    
    def binarize(self, scheme = '02'):
        """
        Convert genotype data to binary
        by mapping 2->1 (2), 0->0 and randomly 
        selecting a 0 or 1 (2) for heterozygotes.
        This emulates the read sampled data found
        in ancient DNA studies
        """
        # randomly sample hets
        p_het = np.sum(self.Y[self.Y == 1.])
        if scheme == '02':
            if p_het > 0:
                samps = np.random.binomial(1, .5, int(p_het)).astype(np.float32) 
                np.place(self.Y, self.y == 1., samps)
        
        if scheme == '01':
            if p_het > 0:
                samps = np.random.binomial(1, .5, int(p_het)).astype(np.float32) 
                np.place(self.Y, self.y == 1., samps)
        # set 2s to 1s
        self.Y[self.Y == 2.] = 1.
        

    def normalize(self, scale, impute = True):
        """
        Mean center the data so each snp has mean 0. 
        If scale true then scale each SNP to have var 1. 
        If impute is true it sets missing sites to 0. i.e. 
        the approx mean. Note for all computations we 
        ignore nans.

        Arguments
        ---------
        scale : bool
            scale the data to have variance 1

        impute : bool
            predict missing sites with the expected mean 0
        """
        if scale:
            # compute mean genotype for every SNP
            mu = np.nanmean(self.Y, axis=1)

            # compute std dev for every SNP
            std = np.nanstd(self.Y, axis=1)

            # center and scale
            Zt = (self.Y.T - mu) / std
            self.Z = Zt.T
        else:
            # compute mean genotype for every SNP
            mu = np.nanmean(self.Y, axis=1)

            # center
            Zt = self.Y.T - mu
            self.Z = Zt.T

        if impute:
            # set nans to zero
            self.Z = np.nan_to_num(self.Z)
            
    def orient_alleles(self, ref_df):
        """
        
        intersect snp_df and ref_df
        find snps that need to be flipped
        flip the snps in the genotpye matrix
        
        Do we need to make this a class? Orienter?
        
        TODO
        
        """
        pass

class AncestryMap(Genotypes):
    """A class for the eigenstrat genotype format
    which consists of a geno, snp, and ind file

    Arguments
    ---------
    prefix : str
        prefix of  path to data files

    Attributes
    ----------
    geno_path : str
        path to eigenstrat geno file

    ind_path : str
        path to eigenstrat ind file

    snp_path : str
        path to eigenstart snp file
    """
    def __init__(self, prefix):
    
        # inherit attributes from Genotypes
        super().__init__()
    
        # prefix to paths
        self.prefix = prefix
        
        # path to eigenstart geno file
        self.geno_path = "{}.geno".format(self.prefix)
        
        # path to eigenstrat ind file
        self.ind_path = "{}.ind".format(self.prefix)
        
        # path to eigenstrat snp file
        self.snp_path = "{}.snp".format(self.prefix)

        # path to eigenstart geno file
        self.packed_geno_path = "{}.packedancestrymapgeno".format(self.prefix)
        
    def _get_snp_dataframe(self):
        """Gets snp level dataframe stored as 
        pandas DataFrame
        """
        # read into pandas df
        self.snp_df = pd.read_csv(self.snp_path, header=None,
                                    sep = '\t')
        
        # check if the dataframe has 6 columns
        if self.snp_df.shape[1] < 6:
             self.snp_df = pd.read_csv(self.snp_path, header=None,
                                    sep = '\s+')
        if self.snp_df.shape[1] < 6:
            raise ValueError("{}.snp must have 6 cols".format(self.snp_path))
        
        # name columns
        self.snp_df.columns = ["SNP", "CHROM", "CM_POS", "POS", "A1", "A2"]
        
        # number of snps
        self.p = self.snp_df.shape[0]

    def _get_inds(self):
        """Get list of individual ids
        """
        ind_df = pd.read_table(self.ind_path, header=None,
                               sep = '\t')
        
        # check if the dataframe has 6 columns
        if ind_df.shape[1] < 3:
            ind_df = pd.read_table(self.ind_path, header=None,
                               sep = '\s+')
        if ind_df.shape[1] < 3:
            raise ValueError("{}.ind must have 3 cols".format(self.ind_path))

        ind_df.columns = ["IND", "SEX", "CLST"]
        
        self.clst = ind_df
        self.inds = ind_df["IND"].tolist()
        
        # number of inds
        self.n = len(self.inds)


class PackedAncestryMap(AncestryMap):
    """Class for packed ancestry map eigenstrat format. 
    The packed format's geno is binary so it requires a 
    couple steps processing before being loaded into a 
    numpy array. We modify some code from
    
    https://github.com/mathii/pyEigenstrat
    
    Arguments
    ---------
    prefix : str
        prefix of  path to data files
    """
    def __init__(self, prefix):

        # inhert attributes from AncestryMap
        super().__init__(prefix)

        # get the snp dataframe
        super()._get_snp_dataframe()

        # get list of individuals
        super()._get_inds()

        # get the genotype matrix
        self._get_genotype_matrix()

    def _get_genotype_matrix(self):
        """Gets the genotype matrix stored as a
        numpy array
        """
        rlen = max(48, int(np.ceil(self.n * 2 / 8)))

        # read in binary 
        self.Y = np.fromfile(self.packed_geno_path, dtype='uint8')[rlen:]
        self.Y.shape = (self.p, rlen)

        # unpack 
        self.Y = np.unpackbits(self.Y, axis=1)[:, :(2 * self.n)]
        self.Y = 2 * self.Y[:, ::2] + self.Y[:, 1::2]

        # convert to float to allow missing data stored as nan
        self.Y = self.Y.astype(np.float)

        # set missing data to nan
        self.Y[self.Y == 3] = np.nan

class UnpackedAncestryMap(AncestryMap):
    """Class for unpacked ancestry map eigenstrat format. The unpacked format's
    geno is not binary so it doesn't require the same processing before 
    being loaded into a numpy array. 
    
    Arguments
    ---------
    prefix : str
        prefix of  path to data files
    """
    def __init__(self, prefix):

        # inhert attributes from AncestryMap
        super().__init__(prefix)

        # get the snp dataframe
        super()._get_snp_dataframe()

        # get list of individuals
        super()._get_inds()
        
        # get the genotype matrix
        self._get_genotype_matrix()

    def _get_genotype_matrix(self):
        """Gets the genotype matrix stored as a
        numpy array
        """
        # read the geno file
        with open(self.geno_path, "r") as f:
            matstr = f.read().replace("\n", "")

        # convert to long array
        ys = np.fromstring(matstr, dtype=np.int8) - 48

        # reshape to p x n matrix
        self.Y = ys.reshape(self.p, self.n).astype(np.float)

        # replace 9s with nans
        self.Y[self.Y == 9] = np.nan

