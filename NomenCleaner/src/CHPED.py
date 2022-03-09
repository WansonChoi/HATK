# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np
import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error
from src.util import Exists
from NomenCleaner.src.HPED import HPED, hasHeader

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


class CHPED(object):

    def __init__(self, _chped, _Nfield=None):
        ### Main Variables ###
        self.chped = None
        self.Nfield = _Nfield

        self.df_chped = None
        self.df_chped_Left = None; self.df_chped_Right = None

        self.HLA_avail = None
        # No `self.HLA_target`.

        self.f_hasHeader = None
        # self.f_hasPrefix = None # All elements(HLA types) have the gene prefix(ex. 'A*').
        # self.f_hasFieldSep = None
        self.isAllCHPEDallele = None

        ### Main Actions ###
        self.setCHPED(_chped)
        self.checkCHPED()


    def setCHPED(self, _chped):
        if Exists(_chped):
            self.chped = _chped
            self.f_hasHeader = hasHeader(_chped)
        else:
            raise HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME +
                "Given CHPED file('{}') can't be found.".format(_chped)
            )


    def checkCHPED(self):
        """
        ## Dimension check
        ## 'HLA_avail' (Header or not?)
        ## isCHPEDallele? (GenePrefix, FieldSep)
        ## Guess N-field(i.e. get the maximum N-field.)
        """
        self.df_chped = pd.read_csv(self.chped, sep='\s+', header=(0 if self.f_hasHeader else None), dtype=str)
        self.df_chped_Left = self.df_chped.iloc[:, :6]
        self.df_chped_Right = self.df_chped.iloc[:, 6:]

        ## Dimension check(== 6 + N*2)
        if self.df_chped_Right.shape[1] % 2:
            raise HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME +
                "# of columns of the given CHPED is {col_all}(6 + {col_HLA}). It must be even number." \
                .format(col_all=self.df_chped.shape[1], col_HLA=self.df_chped_Right.shape[1])
            )

        ## `HLA_avail`
        if self.f_hasHeader:
            self.HLA_avail = np.unique([re.sub(r'_\d$', '', item) for item in self.df_chped_Right.columns])
        else:
            self.HLA_avail = getAvailableHLA(self.df_chped_Right)

        ## isCHPEDallele?
        self.isAllCHPEDallele = np.all(np.vectorize(isCHPEDallele)(self.df_chped_Right.values))

        ## Guess N-field
        self.Nfield = np.max(np.vectorize(guessNfield)(self.df_chped_Right.values))



    def __repr__(self):

        str_hped = \
            "- CHPED file: {}\n".format(self.chped)

        str_hasHeader = \
            "- has Header?: {}\n".format(self.f_hasHeader)

        str_HLA = \
            "- HLA available: {}\n".format(list(self.HLA_avail))

        str_isAllCHPEDallele = \
            "- Are all alleles are CHPED alleles?: {}\n".format(self.isAllCHPEDallele)

        str_Nfield = \
            "- Field format: {}-field\n".format(self.Nfield)

        str_summary = ''.join([
            str_hped,
            str_hasHeader,
            str_HLA,
            str_isAllCHPEDallele,
            str_Nfield
        ]).rstrip('\n')

        return str_summary


    def __bool__(self):
        return Exists(self.chped)


def getAvailableHLA(_df_chped_Right):

    l_RETURN = [] # HLA_avail

    N_col = int(_df_chped_Right.shape[1]/2)

    for i in np.arange(N_col):
        idx1 = 2*i
        idx2 = idx1 + 1

        arr_temp = np.setdiff1d(np.union1d(_df_chped_Right.iloc[:, idx1], _df_chped_Right.iloc[:, idx2]), ['0'])
        # print(arr_temp)

        func = np.vectorize(lambda x : x.split('*')[0])
        arr_HLA = np.unique(func(arr_temp))

        if len(arr_HLA) > 1:
            print(std_WARNING_MAIN_PROCESS_NAME +
                  "There are more than two HLA genes('{}') for {}-th and {}-th columns" \
                  .format(arr_HLA, idx1, idx2))
        else:
            l_RETURN.append(arr_HLA[0])

    return l_RETURN


isCHPEDallele = lambda x : bool(re.match(r'\w+\*\d{2,3}(:\d{2,3})*[A-Z]?', x)) if x != '0' else True
guessNfield = lambda x : len(x.split(':')) if x != '0' else 0



if __name__ == '__main__':

    """
    Wrapper class for the CHPED file('*.chped').

    """

    # Header
    _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220308_BKreport/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.chped"

    # w/o Header
    # _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220308_BKreport/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.subset.8major.chped"

    CH = CHPED(_chped)
    print(CH)

    # df_chped = pd.read_csv(_chped, sep='\s+', header=None, dtype=str)
    # HLA_avail = getAvailableHLA(df_chped.iloc[:, 6:])
    # print(HLA_avail)