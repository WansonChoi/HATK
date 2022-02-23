#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join, exists
import numpy as np
import pandas as pd

from src.util import Exists, printDF, printDict


def getAvailableHLAs(_imgt_dir, _p_data, _hg):

    """
    """

    ## (1) HLA available in IMGT database with exact 3 files(prot, nuc, gen).
    """
    In v3320, for example,
    HLA-DPA2 doesn't have 'DPA2_prot.txt'.
    HLA-DPB2 doesn't have 'DPB2_prot.txt'.
    HLA-DRB3 doesn't have 'DRB3-{prot,nuc}.txt'.
    HLA-DRB4 doesn't have 'DRB4-{prot,nuc}.txt'.
    HLA-H doesn't have 'H_prot.txt'.
    HLA-J doesn't have 'J_prot.txt'.
    HLA-K doesn't have 'K_prot.txt'.
    HLA-L doesn't have 'L_prot.txt'.
    HLA-P doesn't have 'P-{prot,nuc}.txt'.
    HLA-T doesn't have 'T_prot.txt'.
    HLA-V doesn't have 'V_prot.txt'.
    HLA-W doesn't have 'W_prot.txt'.
    HLA-Y doesn't have 'Y_prot.txt'.
    
    """
    l_files = os.listdir(_imgt_dir+'/alignments')
    # print(l_files)

    HLA_all = np.unique([file.split('_')[0] for file in l_files])
    # Only HLAs of which the exact 3 files are available.

    dict_files_prot = getTargetProtFiles(HLA_all, _imgt_dir, "prot") # Mainly because of DRB1. (DRB_prot.txt)
    dict_files_nuc = getTargetProtFiles(HLA_all, _imgt_dir, "nuc")
    dict_files_gen = getTargetProtFiles(HLA_all, _imgt_dir, "gen")

    # printDict(dict_files_prot)
    # printDict(dict_files_nuc)
    # printDict(dict_files_gen)

    arr_HLA_avail_files = \
        [hla for hla in HLA_all if Exists(dict_files_prot[hla]) and Exists(dict_files_nuc[hla]) and Exists(dict_files_gen[hla])]

    arr_HLA_avail_files = np.array(arr_HLA_avail_files)
    # print(arr_HLA_avail_files)


    ## (2) HLA available in Base position information.
    df_HLA_EXON1_start_bp = \
        pd.read_csv(join(_p_data, "HLA_EXON1_START_CODON_POSITIONS_hg{}.txt".format(_hg)),
                    sep='\s+', header=None, dtype=str) \
        .dropna()
    # printDF("df_HLA_EXON1_start_bp", df_HLA_EXON1_start_bp)

    arr_HLA_avail_BP = df_HLA_EXON1_start_bp.iloc[:, 0].values
    # print(arr_HLA_avail_BP)

    return arr_HLA_avail_files, arr_HLA_avail_BP


def getTargetHLAs(_HLA_req, _HLA_avail_files, _HLA_avail_BP):

    """

    Conditions to be a target gene.
    1. Its prot, nuc, and gen files must be in `_imgt_dir`
    2. Its BP and direction information must be given in 'HLA_EXON1_START_CODON_POSITIONS_hgxx.txt' in data folder.
        (cf) The # of available HLA genes are different as to hg version.
    3. It must be in the requested HLA genes.

    """


    f1 = np.isin(_HLA_req, _HLA_avail_files)
    arr_excluded1 = np.array(_HLA_req)[~f1]
    # print(arr_excluded1)

    f2 = np.isin(_HLA_req, _HLA_avail_BP)
    arr_excluded2 = np.array(_HLA_req)[~f2]
    # print(arr_excluded2)

    arr_excluded_all = np.union1d(arr_excluded1, arr_excluded2)
    arr_target = np.setdiff1d(_HLA_req, arr_excluded_all)

    return arr_target, arr_excluded1, arr_excluded2


def getTargetProtFiles(_HLA_target, _imgt_dir, _type):

    _imgt_alignment_dir = _imgt_dir+'/alignments'

    l_files = os.listdir(_imgt_alignment_dir)

    l_files_type = [file for file in l_files if file.endswith('_{type}.txt'.format(type=_type))]
    # print(l_files_prot)


    dict_files_type_target = {hla: None for hla in _HLA_target}

    for hla in _HLA_target:
        filename1 = "{hla}_{type}.txt".format(hla=hla, type=_type)

        if filename1 in l_files_type:
            dict_files_type_target[hla] = join(_imgt_alignment_dir, filename1)
        else:
            dict_files_type_target[hla] = None

        if hla == 'DRB1':# Hard exception handling for HLA-DRB1
            hla2 = re.sub(r'\d$', '', hla)
            filename2 = "{hla}_{type}.txt".format(hla=hla2, type=_type) # mainly due to 'DRB_{prot,nuc}.txt' files.

            if filename2 in l_files_type:
                dict_files_type_target[hla] = join(_imgt_alignment_dir, filename2)


    # for k, v in dict_files_type_target.items():
    #     print("{}: {}".format(k, v))

    return dict_files_type_target



if __name__ == '__main__':

    _imgt_dir = "/home/wansonchoi/sf_VirtualBox_Share/HATK/example/IMGTHLA3320"
    _p_data = "/home/wansonchoi/sf_VirtualBox_Share/HATK/IMGT2Seq/data"
    _hg = "18"

    arr_HLA_avail = getAvailableHLAs(_imgt_dir, _p_data, _hg)
    arr_HLA_target = \
        getTargetHLAs(("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "E", "F", "G", "H", "J"),
                      *arr_HLA_avail)

    # target_files_prot = getTargetProtFiles(arr_HLA_target, _imgt_dir, "prot")
    # target_files_gen = getTargetProtFiles(arr_HLA_target, _imgt_dir, "gen")
    # target_files_nuc = getTargetProtFiles(arr_HLA_target, _imgt_dir, "nuc")


