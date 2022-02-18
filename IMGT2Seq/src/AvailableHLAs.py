#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import numpy as np
import pandas as pd

from src.util import printDF


def getAvailableHLAs(_imgt_dir, _p_data, _hg):

    ## (1) HLA available in sequence files.
    l_files = os.listdir(_imgt_dir+'/alignments')
    # print(l_files)

    arr_HLA_avail_files = np.unique([file.split('_')[0] for file in l_files])
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
    arr_target1 = np.array(_HLA_req)[f1]
    arr_excluded1 = np.array(_HLA_req)[~f1]
    # print(arr_target1)
    # print(arr_excluded1)

    f2 = np.isin(arr_target1, _HLA_avail_BP)
    arr_target2 = arr_target1[f2] # 2 is final output so far.
    arr_excluded2 = arr_target1[~f2]
    # print(arr_target2)
    # print(arr_excluded2)

    return arr_target2, arr_excluded1, arr_excluded2


def getTargetProtFiles(_HLA_target, _imgt_dir, _type):

    _imgt_alignment_dir = _imgt_dir+'/alignments'

    l_files = os.listdir(_imgt_alignment_dir)

    l_files_type = [file for file in l_files if file.endswith('_{type}.txt'.format(type=_type))]
    # print(l_files_prot)


    dict_files_type_target = {hla: None for hla in _HLA_target}

    for hla in _HLA_target:
        hla2 = re.sub(r'\d$', '', hla)
        filename1 = "{hla}_{type}.txt".format(hla=hla, type=_type)
        filename2 = "{hla}_{type}.txt".format(hla=hla2, type=_type) # mainly due to 'DRB_{prot,nuc}.txt' files.

        if filename1 in l_files_type:
            dict_files_type_target[hla] = join(_imgt_alignment_dir, filename1)
        elif filename2 in l_files_type:
            dict_files_type_target[hla] = join(_imgt_alignment_dir, filename2)
        else:
            dict_files_type_target[hla] = None

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


