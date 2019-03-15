# aDNAtools
## Scripts for analyzing ancient DNA

### File Convertion
1. plink2eigenstrat_new.py
2. eigenstrat2plink.py

### Variant Calling
1. pseudo_haploid_pulldown.py

TO DO: 
- plink2eigenstrat_new: test
- pulldown option: disable strand ambiguous SNPs removal (useful when allele0 are always on the forward strand)
- pulldown option: complete I/O. Especially into vcf and packedAncestry formats
- pulldown option: support sex and population labels
- window-f3: window-based f3 for selection scan
- util: merge eigenstrat files based on position
- vcf2eigenstrat: create an .ind file

### Miscellaneous
1. fixed2a0a1.py (for eigenstrat files)
  - Fill in another allele on fixed positions based on an anchor file. It also Removes strand ambiguous variants.

2. vcf2eigenstrat_wgs.py (for vcf files subset into array position)
  - fill not variants not presenting in vcf files in as homozygous as reference. Removes inconsistent variants. 
