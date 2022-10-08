#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists
import pandas as pd
import numpy as np

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.util import Exists, checkFile, hasHeader, getHeader

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class GT(object):
    """
    A CONTAINER class to manage PLINK genotype files(*.{bed,bim,fam}).
    """
    def __init__(self, _file_prefix):

        ### Main Variables ###
        self.bed = checkFile(_file_prefix+'.bed', std_ERROR+"Given PLINK BED file('{}') can't be found.".format(_file_prefix+'.bed'))
        self.bim = checkFile(_file_prefix+'.bim', std_ERROR+"Given PLINK BIM file('{}') can't be found.".format(_file_prefix+'.bim'))
        self.fam = checkFile(_file_prefix+'.fam', std_ERROR+"Given PLINK FAM file('{}') can't be found.".format(_file_prefix+'.fam'))

        self.N_samples, \
        self.f_hasSexInfo, \
        self.f_hasPheInfo = checkFAM(self.fam)

        self.M_markers = checkBIM(self.bim)


    def __repr__(self):

        str_BED = \
            "- BED: {}\n".format(self.bed)
        str_BIM = \
            "- BIM: {}\n".format(self.bim)
        str_FAM = \
            "- FAM: {}\n".format(self.fam)

        str_N_samples = \
            "- # of samples: {}\n".format(self.N_samples)
        str_M_markers = \
            "- # of markers: {}\n".format(self.M_markers)

        str_hasSexInfo = \
            "- has Sex Info?: {}\n".format(self.f_hasSexInfo)
        str_hasPheInfo = \
            "- has Phenotype Info?: {}\n".format(self.f_hasPheInfo)

        str_summary = ''.join([
            str_BED, str_BIM, str_FAM,
            str_N_samples, str_M_markers,
            str_hasSexInfo, str_hasPheInfo
        ]).rstrip('\n')

        return str_summary



def checkFAM(_fam):

    df_fam = pd.read_csv(_fam, sep='\s+', header=None, dtype=str)
    # print("df_fam:\n{}\n".format(df_fam))

    arr_sex = df_fam.iloc[:, 4].values
    f_hasSexInfo = not np.all( np.logical_or((arr_sex == '-9'), (arr_sex == '0')) )
    # print(f_hasSexInfo)

    arr_phe = df_fam.iloc[:, 5].values
    f_hasPheInfo = not np.all( np.logical_or((arr_phe == '-9'), (arr_phe == '0')) )
    # print(f_hasPheInfo)

    return df_fam.shape[0], f_hasSexInfo, f_hasPheInfo


def checkBIM(_bim):

    df_bim = pd.read_csv(_bim, sep='\s+', header=None, dtype=str)

    return df_bim.shape[0]



if __name__ == '__main__':

    ### Genotype class

    ## Mac

    # _gt = '/Users/wansonchoi/Git_Projects/HATK/example/wtccc_filtered_58C_RA.hatk.300+300.hg18.chr6.29-34mb'
    _gt = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged"


    gt = GT(_gt)
    print(gt)
