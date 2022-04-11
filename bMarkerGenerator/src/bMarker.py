#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists
import numpy as np
import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.PLINK import Genotype

class bMarker(Genotype):

    def __init__(self, _file_prefix):

        super().__init__(_file_prefix)


    def checkBIM(self):

        df_bim = pd.read_csv(self.bim, sep='\s+', header=None, dtype=str)

        ## HLA markers

        # HLA
        p_HLA = re.compile(r'^HLA_')
        # f_HLA =

        # AA
        p_AA = re.compile(r'^AA_')

        # SNPS
        p_SNPS = re.compile(r'^SNPS_')

        # INS
        p_INS = re.compile(r'^INS_')





if __name__ == '__main__':

    bmarker = "/Users/wansonchoi/Git_Projects/HATK/example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18"
    r = bMarker(bmarker)
    print(r)
    pass