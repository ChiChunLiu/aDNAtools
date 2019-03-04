# aDNAtools
## Scripts for analyzing ancient DNA

### File Convertion
1. plink2eigenstrat_new.py
2. eigenstrat2plink.py

### Variant Calling
2. pseudo_haploid_pulldown.py

TO DO: 
- plink2eigenstrat_new: test
- pulldown option: disable strand ambiguous SNPs removal (useful when allele0 are always on the forward strand)
- pulldown option: complete I/O. Especially into vcf and packedAncestry formats
- pulldown option: support sex and population labels
- window-f3: window-based f3 for selection scan

# Miscellaneous
3. fixed2a0a1.py (for ancestryMap files)
  - Fill in another allele on fixed positions based on an anchor file. It also Removes strand ambiguous variants.
