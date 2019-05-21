# aDNAtools
## Scripts for analyzing ancient DNA

### File Convertion
1. plink2eigenstrat_new.py
2. eigenstrat2plink.py
3. vcf2eigenstrat_wgs.py
  - convert vcf to eigenstrat files with non-present positions filled as fixed as REF

### Eigenstrat Utils
1. eigenstrat_fixed2a0a1.py (for eigenstrat files)
  - fill in another allele on fixed positions based on an anchor file. It also Removes strand ambiguous variants.

3. eigenstrat_subset.py
  - subset eigenstrat file with population:number table/population/sample ID

4. eigenstrat_QC.py
  - compute fraction of missing entries per individual and per snp

5. eigenstrat_a0asRef.py
  - set a0 as reference

6. eigenstrat_CM_filler.py
  - fill in genetic positions

7. eigenstrat_pca_preprocessing.py
  - produce various filters (e.g. MAF) for SMARTPCA

### Variant Calling
1. pseudo_haploid_pulldown.py
  - pull down bcftools mpileup files into eigenstrat files

### Python utilities
1. util/format
  - forked from Novembre Lab github
3. geo_util
  - manipulate geographic information

TO DO: 
- plink2eigenstrat_new: test
- pulldown option: disable strand ambiguous SNPs removal (useful when allele0 are always on the forward strand)
- pulldown option: complete I/O. Especially into vcf and packedAncestry formats
- pulldown option: support sex and population labels
- window-f3: window-based f3 for selection scan
- util: merge eigenstrat files based on position
- vcf2eigenstrat: create an .ind file
