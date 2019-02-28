# aDNAtools
## Scripts for analyzing ancient DNA
1. (old) mpileup2vcf.py: random pseudo-haploid calls from bcftools mpileup.
2. (latest) pseudo_haploid_pulldown.py: random pseudo-haploid calls from bcftools mpileup.

TO DO: 
- pulldown option: disable strand ambiguous SNPs removal (useful when allele0 are always on the forward strand)
- pulldown option: complete I/O. Especially into vcf and packgedAncestry formats
- pulldown option: support sex and population labels
- window-f3: window-based f3 for selection scan
