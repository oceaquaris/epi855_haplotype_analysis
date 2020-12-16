To replicate results:

To create haplotype data:
1) Un-bgzip src/data/g2f_hybrid_biallelic_maf03_prune99.vcf.gz
2) Run src/haplotype.py
3) Run src/haplotype_data_process.R

To perform GWAS and haplotype association:
1) Un-gzip src/rMVP_in/*.gz
2) Run src/haplotype_analysis.R