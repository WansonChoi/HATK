# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np
import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error
from src.util import Exists

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


class HPED(object):

    def __init__(self, _hped):
        ### Main Variables ###
        self.hped = None
        self.df_hped = None
        self.df_hped_Left = None; self.df_hped_Right = None
        self.HLA_avail = None

        self.f_hasHeader = None
        self.f_hasPrefix = None # All elements(HLA types) have the gene prefix(ex. 'A*').
        self.f_hasFieldSep = None


        ### Main Actions ###
        self.setHPED(_hped)
        self.diagnoseHPED()

        # print(self.__repr__())


    def setHPED(self, _hped):
        if Exists(_hped):
            self.hped = _hped

            self.f_hasHeader = hasHeader(self.hped)

            self.df_hped = pd.read_csv(self.hped, sep='\s+', header=(0 if self.f_hasHeader else None), dtype=str)
            self.df_hped_Left = self.df_hped.iloc[:, :6]
            self.df_hped_Right = self.df_hped.iloc[:, 6:]
        else:
            raise HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME +
                "The given HPED file('{}') can't be found. ".format(self.hped)
            )


    def diagnoseHPED(self):
        """
        # Any allele that has the field separator?
        # Any allele that has the HLA gene prefix?
        # Any G/P-group allele?

        # Guess the field format of the given hped file.
        """

        if self.f_hasHeader:
            ## Dimension check(== 6 + N*2)
            if self.df_hped_Right.shape[1] % 2:
                raise HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "# of columns of the given HPED is {col_all}(6 + {col_HLA}). It must be even number." \
                    .format(col_all=self.df_hped.shape[1], col_HLA=self.df_hped_Right.shape[1])
                )

            ## `HLA_avail`
            self.HLA_avail = np.unique([re.sub(r'_\d$', '', item) for item in self.df_hped_Right.columns])

        else:
            ## Dimension check(== 6 + 8*2  == 22)
            if self.df_hped.shape[1] != 22:
                raise HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "# of columns of a HPED file without a header must be 22(6 + 8*2). The given HPED file('{}') has {}." \
                    .format(self.hped, self.df_hped.shape[1])
                )

            ## `HLA_avail`
            self.HLA_avail = np.array(["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]) # fixed 8 major HLA genes.


        ## `f_hasPrefix`
        # p_GenePrefix = re.compile(r'^\w+\*')
        # func = lambda x :
        # arr =


    def __repr__(self):

        str_hped = \
            "- HPED file: {}\n".format(self.hped)

        str_hasHeader = \
            "- has Header?: {}\n".format(self.f_hasHeader)

        str_HLA = \
            "- HLA available: {}\n".format(list(self.HLA_avail))

        str_summary = ''.join([
            str_hped,
            str_hasHeader,
            str_HLA
        ]).rstrip('\n')

        return str_summary


    def __bool__(self):
        return Exists(self.hped)


def hasHeader(_hped):
    p_Header = re.compile(r'^FID\s+IID')
    with open(_hped, 'r') as f_hped:
        line_1st = f_hped.readline()
        return bool(p_Header.match(line_1st))



if __name__ == '__main__':

    """
    Wrapper class for the HPED file('*.hped').
    
    """


    hped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220221/Dummy.N50.hped"
    H = HPED(hped)
    print(H)