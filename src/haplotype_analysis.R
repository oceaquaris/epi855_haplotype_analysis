#!/usr/bin/env Rscript
setwd("C:/Users/rs14/Documents/classes/epi855_statistical_genetics/project/src/")

# load libraries
library(rMVP)

# define file names we'll be using
fname_geno_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.geno.tsv"
fname_haplo_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.haplo.tsv"
fname_htaxa_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.htaxa.tsv"
fname_hpos_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.hpos.tsv"
fname_hmap_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.hmap.tsv"
fname_gmap_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.gmap.tsv"
fname_ind_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.ind.tsv"
fname_kin_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.kin.tsv"
fname_feff_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.feff.tsv"
fname_pheno_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.pheno.tsv"

# define prefixes we'll be using
prefix_geno <- "genoAnalysis"
prefix_haplo <- "haploAnalysis"

##################################################################
########################## prepare data ##########################
##################################################################

# prepare genotype analysis data
MVP.Data(
    fileNum = fname_geno_rMVP,
    fileMap = fname_gmap_rMVP,
    filePhe = fname_pheno_rMVP,
    fileInd = fname_ind_rMVP,
    fileKin = fname_kin_rMVP,
    filePC  = fname_feff_rMVP,
    out = paste0("rMVP_bin/", prefix_geno),
    sep.num = '\t',
    sep.map = '\t',
    sep.phe = '\t',
    sep.kin = '\t',
    sep.pc  = '\t',
    SNP.impute = NULL,
    verbose = TRUE
)

# prepare haplotype analysis data
MVP.Data(
    fileNum = fname_haplo_rMVP,
    fileMap = fname_hmap_rMVP,
    filePhe = fname_pheno_rMVP,
    fileInd = fname_ind_rMVP,
    fileKin = fname_kin_rMVP,
    filePC  = fname_feff_rMVP,
    out = paste0("rMVP_bin/", prefix_haplo),
    sep.num = '\t',
    sep.map = '\t',
    sep.phe = '\t',
    sep.kin = '\t',
    sep.pc  = '\t',
    SNP.impute = NULL,
    verbose = TRUE
)

##################################################################
########################## analyze data ##########################
##################################################################

################
# Perform GWAS #
################
# load genotype data
geno <- attach.big.matrix(paste0("rMVP_bin/", prefix_geno, ".geno.desc"))

# load covariates
covar <- attach.big.matrix(paste0("rMVP_bin/", prefix_geno, ".pc.desc"))

# load map information
map <- read.table(paste0("rMVP_bin/", prefix_geno, ".geno.map"), header = FALSE)

# load phenotype date
pheno <- read.table(paste0("rMVP_bin/", prefix_geno, ".phe"), header = TRUE)

# load kinship matrix
kinship <- attach.big.matrix(paste0("rMVP_bin/", prefix_geno, ".kin.desc"))

# custom covariate numbers for each model
covars = list(
    list(covar[,c(1,4)], covar[,4]), # Anthesis
    list(covar[,c(1,4)], covar[,4]), # Silking
    list(covar[,c(1,4)], covar[,4]), # Plant height
    list(covar[,4],      covar[,4])  # Ear height
)

# conduct GWAS
for(i in 2:ncol(pheno)) {
    mvp <- MVP(
        phe = pheno[,c(1,i)],
        geno = geno,
        map = map,
        K = kinship,
        CV.MLM = covars[[i-1]][[1]],
        CV.FarmCPU = covars[[i-1]][[2]],
        vc.method = "EMMA",
        method = c("MLM","FarmCPU"),
        outpath = paste0("rMVP_out/", prefix_geno, "/")
    )
}

#######################################
# Perform Haplotype Association Study #
#######################################
# load haplotype data
geno <- attach.big.matrix(paste0("rMVP_bin/", prefix_haplo, ".geno.desc"))

# load covariates
covar <- attach.big.matrix(paste0("rMVP_bin/", prefix_haplo, ".pc.desc"))

# load map information
map <- read.table(paste0("rMVP_bin/", prefix_haplo, ".geno.map"), header = FALSE)

# load phenotype date
pheno <- read.table(paste0("rMVP_bin/", prefix_haplo, ".phe"), header = TRUE)

# load kinship matrix
kinship <- attach.big.matrix(paste0("rMVP_bin/", prefix_haplo, ".kin.desc"))

# custom covariate numbers for each model
covars = list(
    list(covar[,4],      covar[,4]), # Anthesis
    list(covar[,4],      covar[,4]), # Silking
    list(covar[,4],      covar[,4]), # Plant height
    list(covar[,c(1,4)], covar[,4])  # Ear height
)

# conduct GWAS
for(i in 2:ncol(pheno)) {
    mvp <- MVP(
        phe = pheno[,c(1,i)],
        geno = geno,
        map = map,
        K = kinship,
        CV.MLM = covars[[i-1]][[1]],
        CV.FarmCPU = covars[[i-1]][[2]],
        vc.method = "EMMA",
        method = c("MLM","FarmCPU"),
        outpath = paste0("rMVP_out/", prefix_haplo, "/")
    )
}
