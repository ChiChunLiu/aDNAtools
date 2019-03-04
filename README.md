# aDNAtools
## Scripts for analyzing ancient DNA

### File Convertion
1. plink2eigenstrat.py
2. eigenstrat2plink.py
3. fixed2a0a1.py
  - Fill in another allele on fixed positions based on an anchor file. Remove strand ambiguous variants.

### Variant Calling
2. pseudo_haploid_pulldown.py

TO DO: 
- plink2eigenstrat.py: chunk chromosome to make it memory efficient
- pulldown option: disable strand ambiguous SNPs removal (useful when allele0 are always on the forward strand)
- pulldown option: complete I/O. Especially into vcf and packedAncestry formats
- pulldown option: support sex and population labels
- window-f3: window-based f3 for selection scan
