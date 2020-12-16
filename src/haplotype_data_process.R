#!/usr/bin/env Rscript
setwd("C:/Users/rs14/Documents/classes/epi855_statistical_genetics/project/src/")

# load libraries
library(data.table)
library(vcfR)

# define file names we'll be using
fname_geno <- "data/g2f_hybrid_biallelic_maf03_prune99.vcf"
fname_haplo <- "data/hmat.txt"
fname_hposgz <- "data/hpos.txt"
fname_htaxagz <- "data/htaxa.txt"

fname_geno_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.geno.tsv"
fname_haplo_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.haplo.tsv"
fname_htaxa_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.htaxa.tsv"
fname_hpos_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.hpos.tsv"
fname_hmap_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.hmap.tsv"
fname_gmap_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.gmap.tsv"
fname_ind_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.ind.tsv"
fname_kin_rMVP <- "rMVP_in/g2f_hybrid_biallelic_maf03_prune99.kin.tsv"

fname_geno_Rdata <- "Rdata/g2f_hybrid_biallelic_maf03_prune99.geno.mat.Rdata"
fname_haplo_Rdata <- "Rdata/g2f_hybrid_biallelic_maf03_prune99.haplo.mat.Rdata"
fname_htaxa_Rdata <- "Rdata/g2f_hybrid_biallelic_maf03_prune99.htaxa.df.Rdata"
fname_hpos_Rdata <- "Rdata/g2f_hybrid_biallelic_maf03_prune99.hpos.df.Rdata"
fname_kin_Rdata <- "Rdata/g2f_hybrid_biallelic_maf03_prune99.kin.mat.Rdata"


##################################################################
########################### read files ###########################
##################################################################

#################
# read vcf file #
#################
# read file definition
vcf <- read.vcfR(fname_geno, verbose=FALSE)

# using as.numeric = TRUE results in undefined behavior due to / in file.
# we must read fields manually and use a hash table to extract values
# extract strings: p x n = nmarker x ntaxa
geno_str <- extract.gt(vcf, as.numeric = FALSE)

# create empty matrix: n x p = ntaxa x nmarker
geno <- matrix(0, nrow = ncol(geno_str), ncol = nrow(geno_str))
colnames(geno) <- rownames(geno_str)
rownames(geno) <- colnames(geno_str)

# create hash table
hashtab <- c(0,1,1,2)
names(hashtab) <- c("0/0","0/1","1/0","1/1")

# populate fields with genotypes
for(j in 1:ncol(geno_str)) {
    geno[j,] <- hashtab[geno_str[,j]]
}

# save matrix for easy reading
write.table(t(geno), file = fname_geno_rMVP, sep = '\t', row.names = FALSE, col.names = FALSE)
save(geno, file=fname_geno_Rdata)

gmap_df <- data.frame(SNP_ID = getID(vcf), Chromosome = getCHROM(vcf), Position = getPOS(vcf))
write.table(gmap_df, file = fname_gmap_rMVP, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

#########################
# read haplotype matrix #
#########################
# read file using fread for speed
haplo <- as.matrix(fread(fname_haplo, header = FALSE))

# read taxa names
htaxa <- read.table(fname_htaxagz, header = FALSE)

# read haplotype positions
hpos <- read.table(fname_hposgz, header = FALSE)

# get positions
v <- hpos[,2]
hposgrp <- vector(mode = "integer", length = length(v))
g = 1
hposgrp[1] <- g
for(i in 2:length(v)) {
    g = ifelse(v[i-1] == v[i], g + 1, 1)
    hposgrp[i] <- g
}

hapID <- paste0("H", hpos[,1], "_", hpos[,2], "_", hposgrp)

# name column and row names
rownames(haplo) <- htaxa[,1]
colnames(haplo) <- hapID

# save matrix for easy reading
write.table(t(haplo), file = fname_haplo_rMVP, sep = '\t', row.names = FALSE, col.names = FALSE)
save(haplo, file = fname_haplo_Rdata)

# save map file for rMVP
map_df <- data.frame(SNP_ID = as.character(hapID), Chromosome = hpos[,1], Position = hpos[,2])
write.table(map_df, file = fname_hmap_rMVP, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# save individual name file
ind_df <- data.frame(Name = htaxa[,1])
write.table(ind_df, file = fname_ind_rMVP, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# calculate IBS kinship matrix
invnloci = 1.0 / ncol(geno)                                 # calculate inverse number of loci
X = geno - 1                                                # convert 0,1,2 to -1,0,1
kinship <- 0.5 * ((invnloci * tcrossprod(X, X)) + 1.0)      # calculate kinship

# save kinship to tile
write.table(kinship, file = fname_kin_rMVP, sep = '\t', row.names = FALSE, col.names = FALSE)
save(kinship, file = fname_kin_Rdata)

# save data frams for easy reading
save(htaxa, file=fname_htaxa_Rdata)
save(hpos, file=fname_hpos_Rdata)
