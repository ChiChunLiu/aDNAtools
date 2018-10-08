# requires numpy, pandas, argparse, scikit-allel
import numpy as np
import pandas as pd
import argparse, os, allel

def random_draw(snp_list, pileup):
    """
    Draw a random allele at every given position from reads
    
    Parameters
    ----------
    snp_list : array_like
        table with 5 columns: id, chr, pos, ref, alt
    pileup :
        Sser should create from:
        bcftools mpileup --ignore-RG -B -q30 -Q30 \
        -T pos.tsv.gz \
        -f hs37d5.fa \
        -b bam_list \ 
        -a AD,DP | \
        bcftools norm -m-both -Oz -o out.vcf.gz
        
    Returns
    -------
    pseudo: ndarray(dtype=float, ndim=2)
        allele randomly drawn
    ad: ndarray(dtype=float, ndim=2)
        allele depth
    snp_out: DataFrame
        pandas data frame with 5 columns: 
        id, chr, pos, ref, alt
    samples: list
         samples
    """
    # read in mpileup file
    callset = allel.read_vcf(pileup, fields = '*')
    ad = callset['calldata/AD'][:,:,0:2]

    # parse mpileup info with scikit-allel and add one index column
    var_info = pd.DataFrame({'chr':callset['variants/CHROM'], 'pos': callset['variants/POS'], 
                  'ref': callset['variants/REF'], 'alt': callset['variants/ALT'][:,0]})
    var_info['index1'] = var_info.index
    
    sites_in = pd.read_table(snp_list, names = ['var_id','chr','pos','ref','alt_input'])
    sites_in['chr'] = sites_in['chr'].astype(str)

    # joining tables (these are all dealing with pileup format)
    new_df = sites_in.merge(var_info, how='left', on = ['chr','pos','ref'])
    # mark duplicated positions and keep either rows with consistent alt allele or no duplicate
    new_df['dup'] = new_df.duplicated(subset = ['chr','pos'], keep = False)
    new_df_ddup = new_df[(new_df.dup == False) | (new_df.alt == new_df.alt_input)]
    
    # above wrongly removes sites with reads of other alleles (machine error or multiallelic)
    # deal with that with the below. Suboptimal, to be changed.
    new_df['bad'] = False
    idx = new_df[~new_df['var_id'].isin(new_df_ddup['var_id'])].index
    new_df.loc[idx,'bad'] = True
    snp_df = new_df[(new_df.dup == False) | (new_df.alt == new_df.alt_input) | ((new_df.bad == True) & (new_df.alt == '<*>'))]
    snp_df.reset_index(drop = True, inplace = True)

    # original and new index of SNPs with reads
    i = snp_df[snp_df['index1'].isnull() == False].index1.astype(int)
    j = snp_df[snp_df['index1'].isnull() == False].index
    ad = ad[i,:,:]
    ad_tot = ad.sum(axis = 2)
    
    # random draw of alleles
    pseudo = -1 * np.ones(ad_tot.shape)
    r, c = np.where(ad_tot > 0)
    pseudo[r,c] = np.random.binomial(1, ad[:,:,1][r,c]/ad_tot[r,c])
    
    # outputs
    pseudo_all = -1 * np.ones((snp_df.shape[0], pseudo.shape[1]))
    pseudo_all[j,:] = pseudo
    ad_all = np.zeros((snp_df.shape[0], pseudo.shape[1], 2))
    ad_all[j,:,:] = ad
    ad_all = ad_all.astype(int)
    
    snp_out = snp_df[['var_id', 'chr', 'pos', 'ref', 'alt_input']]
    snp_out.columns = ['id', 'chr', 'pos', 'ref', 'alt']
    
    return pseudo_all, ad_all, snp_out, callset['samples']


def array2gt(x):
    if x == 0:
        return "0/0"
    elif x == 1:
        return "1/1"
    elif x == -1:
        return "./."
    
def vcf_gt(x):
    return list(map(array2gt, x))

def vcf_ad(x):
    return list(map(lambda x,y: str(x) + ',' + str(y), x[:, 0], x[:, 1]))

def write_vcf_rec(snp_info, gt, ad):
    rec = '\t'.join(map(str, snp_info)) + '\t'
    rec += '\t'.join([".", ".", ".", "GT:AD"]) + '\t'
    ad = vcf_ad(ad)
    gt = vcf_gt(gt)
    rec += '\t'.join(list(map(lambda x,y: x + ':' + y, gt, ad)))
    vcf_print(rec)

   
def vcf_header(sample_list):
    vcf_print("##fileformat=VCFv4.2", overwrite = True)
    vcf_print("##Source=random_draw.py")
    vcf_print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    vcf_print("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth\">")
    s = '\t'.join(sample_list)
    vcf_print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + s)

def vcf_print(x, overwrite = False):
    if overwrite:
        print(x, file = open(args.output, "w"))
    else:
        print(x, file = open(args.output, "a"))
  
if __name__== "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--snp', type = str, default = "", help = "biallelic SNP table")
    parser.add_argument('-p', '--pileup', type = str, default = "", help = "bcftools mpileup")
    parser.add_argument('-o', '--output', type = str, default = "", help = "output")   
    args = parser.parse_args()
        
    # stdout vcf 
    gts, ads, snps, samples = random_draw(args.snp, args.pileup)
    vcf_header(samples)
    for index, rec in snps.iterrows():
        snp = rec[['chr','pos','id','ref','alt']]
        gt = gts[index,:]
        ad = ads[index,:,:]
        write_vcf_rec(snp, gt, ad)
