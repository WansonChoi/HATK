#-*- coding: utf-8 -*-
import os, sys, re
from os.path import join, basename, dirname
import numpy as np
import pandas as pd


def checkBIM(_bim_HATK):

    """
    1. Count the number of each type of markers. (SNP, HLA, AA, SNPS, INS)
    2. Acquire Unique HLA gene set('HLA_target'). (ex. ["A", "B", "C", ..., "DRB1"])
    """

    df_bim = pd.read_csv(_bim_HATK, sep='\s+', header=None, dtype=str)
    # print("df_bim:\n{}\n".format(df_bim))


    ## df_bim => HLA + AA + SNPS + INS + (intergenic) SNP.

    # HLA
    """
    (old): "HLA_DQB1_02", "HLA_DQB1_0201", "HLA_A_0101"
    (new; HATK): "HLA_F*01:01", "HLA_A*01:04N", "HLA_A*01:100", "HLA_A*01:123N", "HLA_DOA*01:01"
    """
    p_HLA = re.compile(r'^HLA_([A-Z]+\d?)')
    f_HLA = np.array([bool(p_HLA.match(label)) for label in df_bim.iloc[:, 1].values])
    df_bim_HLA = df_bim.loc[f_HLA, :]
    # print("df_bim_HLA:\n{}\n".format(df_bim_HLA))

    # AA
    """
    (old): "AA_A_246_30020067", "AA_A_-15_30018338_L", "AA_DQB1_203_32737170_V"
    (new): "AA_F_-9_29723501_exon1", "AA_DRB1_233_32580271_exon5_R", 
    """
    p_AA = re.compile(r'^AA_([A-Z]+\d?)')
    f_AA = np.array([bool(p_AA.match(label)) for label in df_bim.iloc[:, 1].values])
    df_bim_AA = df_bim.loc[f_AA, :]
    # print("df_bim_AA:\n{}\n".format(df_bim_AA))

    # SNPS
    """
    (old): SNP_DQB1_32737107, SNP_DQB1_32740579_G, SNP_DRB1_32657335
    (new): SNPS_DQA1_22_32637480_exon1, SNPS_DPA1_3181_33070390_intron1
    
    (cf): "SNP_A-2064274" <- just a (intergenic) SNP.
    """
    p_SNPS = re.compile(r'^SNPS?_([A-Z]+\d?)_')
    f_SNPS = np.array([bool(p_SNPS.match(label)) for label in df_bim.iloc[:, 1].values])
    df_bim_SNPS = df_bim.loc[f_SNPS, :]
    # print("df_bim_SNPS:\n{}\n".format(df_bim_SNPS))

    # INS
    """
    (old): INS_DQB1_226x227_32736002_PQGPPPAG
    (new): INS_SNPS_DQA1_3126x3127_32640584, INS_AA_C_300x301_31270016
    """
    p_INS = re.compile(r'^INS_(AA_|SNPS_)?([A-Z]+\d?)')
    f_INS = np.array([bool(p_INS.match(label)) for label in df_bim.iloc[:, 1].values])
    df_bim_INS = df_bim.loc[f_INS, :]
    # print("df_bim_INS:\n{}\n".format(df_bim_INS))

    # SNP
    f_SNP = ~(f_HLA | f_AA | f_SNPS | f_INS)
    df_bim_SNP = df_bim.loc[f_SNP, :]
    # print("df_bim_SNP:\n{}\n".format(df_bim_SNP))


    ## Unique HLA genes.
    get_HLAgenes = np.vectorize(lambda x, p, idx: p.match(x).group(idx))

    arr_HLA = list(get_HLAgenes(df_bim_HLA.iloc[:, 1].values, p_HLA, 1)) if df_bim_HLA.shape[0] > 0 else []
    arr_AA = list(get_HLAgenes(df_bim_AA.iloc[:, 1].values, p_AA, 1)) if df_bim_AA.shape[0] > 0 else []
    arr_SNPS = list(get_HLAgenes(df_bim_SNPS.iloc[:, 1].values, p_SNPS, 1)) if df_bim_SNPS.shape[0] > 0 else []
    arr_INS = list(get_HLAgenes(df_bim_INS.iloc[:, 1].values, p_INS, 2)) if df_bim_INS.shape[0] > 0 else []

    l_HLA_target = list(np.unique(arr_HLA + arr_AA + arr_SNPS + arr_INS))
    # print(l_HLA_target)


    return df_bim.shape[0], df_bim_SNP.shape[0], \
           df_bim_HLA.shape[0], df_bim_AA.shape[0], df_bim_SNPS.shape[0], df_bim_INS.shape[0], l_HLA_target



if __name__ == '__main__':

    bim_HATK = "/Users/wansonchoi/Git_Projects/HATK/tests/20221007_bMG/wtccc_filtered_58C_RA.hatk.300+300.hg18.chr6.29-34mb.with38.bim"

    # old version
    # bim_HATK = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.bim"


    checkBIM(bim_HATK)