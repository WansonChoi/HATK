#-*- coding: utf-8 -*-

import os, sys, re
import numpy as np

def getDictHLAs(_dict_seq):

    with open(_dict_seq, 'r') as f_dict_seq:
        for line in f_dict_seq:
            l = line.split()
            yield l[0] # HLA allele (ex. A*01:01)


def getHLAgenes(_dict_seq) -> np.ndarray:
    return np.unique([allele.split('*')[0] for allele in getDictHLAs(_dict_seq)]) # is both `HLA_avail` and `HLA_target`.

def getMaxNfield(_dict_seq) -> int:
    return np.max([len(allele.split('*')[1].split(':')) for allele in getDictHLAs(_dict_seq)])

def checkSameHLAset(_dict_AA_seq, _dict_SNPS_seqs) -> bool:
    HLA_target_AA = getHLAgenes(_dict_AA_seq)
    print(HLA_target_AA)
    HLA_target_SNPS = getHLAgenes(_dict_SNPS_seqs)
    print(HLA_target_SNPS)

    return np.array_equal(HLA_target_AA, HLA_target_SNPS)


if __name__ == '__main__':

    dict_AA_seq = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/IMGT3470_dict_hg18/HLA_DICTIONARY_AA.hg18.imgt3470.txt"
    dict_AA_map = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/IMGT3470_dict_hg18/HLA_DICTIONARY_AA.hg18.imgt3470.map"
    dict_SNPS_seq = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/IMGT3470_dict_hg18/HLA_DICTIONARY_SNPS.hg18.imgt3470.txt"
    dict_SNPS_map = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/IMGT3470_dict_hg18/HLA_DICTIONARY_SNPS.hg18.imgt3470.map"

    tmp1 = getHLAgenes(dict_AA_seq)
    print(tmp1)

    tmp2 = getMaxNfield(dict_AA_seq)
    print(tmp2)

    tmp3 = checkSameHLAset(dict_AA_seq, dict_SNPS_seq)
    print(tmp3)
    pass