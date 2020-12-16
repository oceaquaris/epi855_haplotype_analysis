#!/usr/bin/env python3
import cyvcf2
import numpy

def extract_vcf(fname):
    """
    Extract genotype information.
    Parameters
    ----------
    fname : str
        File name
    """
    vcf = cyvcf2.VCF(fname) # make VCF iterator
    taxa = vcf.samples      # extract taxa names from vcf header (list)
    geno_dict = {}          # make empty dict to store extracted values
    ref_dict = {}           # reference state
    alt_dict = {}           # alternative state
    pos_dict = {}           # physical position
    name_dict = {}          # marker name
    for variant in vcf:     # iterate through VCF file and accumulate variants
        chr = variant.CHROM # get CHROM
        if chr not in geno_dict.keys():     # add CHROM key if needed
            geno_dict[chr] = []
            ref_dict[chr] = []
            alt_dict[chr] = []
            pos_dict[chr] = []
            name_dict[chr] = []
        # append entries
        geno_dict[chr].append(variant.genotype.array()[:,0:2].T)
        ref_dict[chr].append(variant.REF)
        alt_dict[chr].append(variant.ALT[0])
        pos_dict[chr].append(variant.POS)
        name_dict[chr].append(variant.ID)
    # convert all geno_dict entries into numpy arrays
    for chr in geno_dict.keys():
        geno_dict[chr] = numpy.stack(geno_dict[chr], 2)
        pos_dict[chr] = numpy.int64(pos_dict[chr])
    print("read file %s" % fname)
    return geno_dict, ref_dict, alt_dict, pos_dict, name_dict, taxa

def group_haplotype(geno_dict, pos_dict, wsize):
    haplo_dict = {}
    state_dict = {}
    hpos_dict = {}
    hpos_start_dict = {}
    hpos_stop_dict = {}
    for chr in geno_dict.keys():
        if chr not in haplo_dict.keys():
            haplo_dict[chr] = []
            state_dict[chr] = []
            hpos_dict[chr] = []
            hpos_start_dict[chr] = []
            hpos_stop_dict[chr] = []
        nmkr = geno_dict[chr].shape[2]
        for st,sp in zip(range(0,nmkr-wsize), range(wsize,nmkr)):
            x = geno_dict[chr][:,:,st:sp].reshape(-1,wsize) # get block
            uhap, hcode = numpy.unique(x, return_inverse = True, axis = 0)
            haplo_dict[chr].append(hcode.reshape(2,-1,1))
            state_dict[chr].append(uhap)
            hpos_dict[chr].append(numpy.mean(pos_dict[chr][st:sp]))
            hpos_start_dict[chr].append(numpy.min(pos_dict[chr][st:sp]))
            hpos_stop_dict[chr].append(numpy.max(pos_dict[chr][st:sp]))
        print("group_haplotype: CHROM %s done." % chr)
    for chr in haplo_dict.keys():
        haplo_dict[chr] = numpy.concatenate(haplo_dict[chr], 2)
    return haplo_dict, state_dict, hpos_dict, hpos_start_dict, hpos_stop_dict

def make_haplotype_matrix(haplo_dict, state_dict, hpos_dict, hpos_start_dict, hpos_stop_dict):
    nhap = sum(len(e) for k in state_dict.keys() for e in state_dict[k])
    nind = haplo_dict[list(haplo_dict.keys())[0]].shape[1]
    hmat = numpy.zeros((nind, nhap), dtype = 'int16')
    chrom = []
    state = []
    pos = []
    pos_start = []
    pos_stop = []
    k = 0
    for chr in haplo_dict.keys():
        mat = haplo_dict[chr]
        sd = state_dict[chr]
        pd = hpos_dict[chr]
        ad = hpos_start_dict[chr]
        bd = hpos_stop_dict[chr]
        ncol = mat.shape[2]
        for j in range(ncol):
            for i in range(len(sd[j])):
                state.append(sd[j][i])
                pos.append(pd[j])
                pos_start.append(ad[j])
                pos_stop.append(bd[j])
                chrom.append(chr)
            for i in range(nind):
                for p in range(2):
                    hmat[i,k + mat[p,i,j]] += 1
            k += len(sd[j])
        print("make_haplotype_matrix: CHROM %s done." % chr)
    return hmat, state, chrom, pos, pos_start, pos_stop

if __name__ == "__main__":
    geno_dict, ref_dict, alt_dict, pos_dict, name_dict, taxa = extract_vcf("data/g2f_hybrid_biallelic_maf03_prune99.vcf")
    haplo_dict, state_dict, hpos_dict, hpos_start_dict, hpos_stop_dict = group_haplotype(geno_dict, pos_dict, 10)
    hmat, state, chrom, pos, pos_start, pos_stop = make_haplotype_matrix(haplo_dict, state_dict, hpos_dict, hpos_start_dict, hpos_stop_dict)
    numpy.savetxt("data/hmat.txt", hmat, fmt='%i', delimiter='\t')
    with open('data/htaxa.txt', 'w') as f:
        for e in taxa:
            f.write("%s\n" % e)
    with open('data/hstate.txt', 'w') as f:
        for e in state:
            f.write("%s\n" % e)
    with open('data/hpos.txt', 'w') as f:
        for e,ee,eee,eeee in zip(chrom, pos, pos_start, pos_stop):
            f.write("%s\t%s\t%s\t%s\n" % (e,ee,eee,eeee))
