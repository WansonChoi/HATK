#-*- coding: utf-8 -*-
import os, sys, re
from os.path import join, basename, dirname
import numpy as np
import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.PLINK import Genotype
from src.util import Exists

class bMarker(Genotype):

    def __init__(self, _file_prefix):

        ### Main Variables ###
        self.FRQ = _file_prefix+'.FRQ.frq' if Exists(_file_prefix+'.FRQ.frq') else None

        # SNP marker
        self.m_SNP = -1

        # HLA marker
        self.m_HLA = -1
        self.m_AA = -1
        self.m_SNPS = -1
        self.m_INS = -1


        ### Main Actions ###
        super().__init__(_file_prefix)
        self.checkBIM()


    def checkBIM(self):

        df_bim = pd.read_csv(self.bim, sep='\s+', header=None, dtype=str)
        # print("df_bim:\n{}\n".format(df_bim))

        ## HLA markers

        # HLA
        p_HLA = re.compile(r'^HLA_')
        f_HLA = np.array([bool(p_HLA.match(label)) for label in df_bim.iloc[:, 1].values])
        df_bim_HLA = df_bim.loc[f_HLA, :]
        # print("df_bim_HLA:\n{}\n".format(df_bim_HLA))

        # AA
        p_AA = re.compile(r'^AA_')
        f_AA = np.array([bool(p_AA.match(label)) for label in df_bim.iloc[:, 1].values])
        df_bim_AA = df_bim.loc[f_AA, :]
        # print("df_bim_AA:\n{}\n".format(df_bim_AA))

        # SNPS
        p_SNPS = re.compile(r'^SNPS_')
        f_SNPS = np.array([bool(p_SNPS.match(label)) for label in df_bim.iloc[:, 1].values])
        df_bim_SNPS = df_bim.loc[f_SNPS, :]
        # print("df_bim_SNPS:\n{}\n".format(df_bim_SNPS))

        # INS
        p_INS = re.compile(r'^INS_')
        f_INS = np.array([bool(p_INS.match(label)) for label in df_bim.iloc[:, 1].values])
        df_bim_INS = df_bim.loc[f_INS, :]
        # print("df_bim_INS:\n{}\n".format(df_bim_INS))

        # SNP
        f_SNP = ~(f_HLA | f_AA | f_SNPS | f_INS)
        df_bim_SNP = df_bim.loc[f_SNP, :]
        # print("df_bim_SNP:\n{}\n".format(df_bim_SNP))


        self.M_markers = df_bim.shape[0]
        self.m_SNP = df_bim_SNP.shape[0]
        self.m_HLA = df_bim_HLA.shape[0]
        self.m_AA = df_bim_AA.shape[0]
        self.m_SNPS = df_bim_SNPS.shape[0]
        self.m_INS = df_bim_INS.shape[0]


    def __repr__(self):

        str_BED = \
            "- BED: {}\n".format(self.bed)
        str_BIM = \
            "- BIM: {}\n".format(self.bim)
        str_FAM = \
            "- FAM: {}\n".format(self.fam)
        str_FRQ = \
            "- FRQ: {}\n".format(self.FRQ)

        str_N_samples = \
            "- # of samples: {}\n".format(self.N_samples)
        str_M_markers = \
            "- # of Total markers: {}\n".format(self.M_markers)

        str_m_SNP = \
            "- # of SNPs: {}\n".format(self.m_SNP)
        str_m_HLA = \
            "- # of HLA markers: {}\n".format(self.m_HLA)
        str_m_AA = \
            "- # of Amino acid markers: {}\n".format(self.m_AA)
        str_m_SNPS = \
            "- # of Intragenic DNA markers: {}\n".format(self.m_SNPS)
        str_m_INS = \
            "- # of Insertion markers: {}\n".format(self.m_INS)

        str_hasSexInfo = \
            "- has Sex Info?: {}\n".format(self.f_hasSexInfo)
        str_hasPheInfo = \
            "- has Phenotype Info?: {}\n".format(self.f_hasPheInfo)


        str_summary = ''.join([
            str_BED, str_BIM, str_FAM, str_FRQ,
            str_N_samples, str_M_markers, str_m_SNP, str_m_HLA, str_m_AA, str_m_SNPS, str_m_INS,
            str_hasSexInfo, str_hasPheInfo
        ]).rstrip('\n')

        return str_summary




if __name__ == '__main__':

    bmarker = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_bMG_20220913/wtccc_filtered_58C_RA.hatk.300+300.hg18.chr6.29-34mb"
    r = bMarker(bmarker)
    print(r)
    pass