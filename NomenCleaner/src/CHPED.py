# -*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np
import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error
from src.util import Exists, checkFile, FieldFormat2Label
from NomenCleaner.src.HPED import HPED, hasHeader

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


class CHPED(object):

    def __init__(self, _chped):

        ### Main Variables ###
        self.chped = checkFile(_chped, std_ERROR+"Given CHPED file('{}') can't be found.".format(_chped))
        self.f_hasHeader = hasHeader(self.chped)
        self.HLA_avail = None
        self.field_format = None

        # DataFrame
        self.df_chped = None
        self.df_chped_Left = None
        self.df_chped_Right = None

        # Flags
        # self.f_hasPrefix = None # All elements(HLA types) have the gene prefix(ex. 'A*').
        # self.f_hasFieldSep = None
        self.isAllCHPEDallele = None


        ### Main Actions ###
        self.checkCHPED()


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
                std_ERROR +
                "CHPED must have even number of columns(6 + 2*(# of HLA genes)).\n" 
                "Meanwhile, Given CHPED file('{}') has {} columns." \
                .format(self.chped, self.df_chped.shape[1])
            )

        ## `HLA_avail`
        if self.f_hasHeader:
            self.HLA_avail = list(np.unique([re.sub(r'_\d$', '', item) for item in self.df_chped_Right.columns]))
        else:
            self.HLA_avail = getAvailableHLA(self.df_chped_Right)

        ## isCHPEDallele?
        self.isAllCHPEDallele = np.all(np.vectorize(isCHPEDallele)(self.df_chped_Right.values))

        ## Guess N-field
        self.field_format = np.max(np.vectorize(guessFieldFormat)(self.df_chped_Right.values))


    def __repr__(self):

        str_chped = \
            "- CHPED file: {}\n".format(self.chped)

        str_hasHeader = \
            "- has Header?: {}\n".format(self.f_hasHeader)

        str_HLA = \
            "- HLA:\n" \
            "   (requested): {}\n".format(list(self.HLA_avail))

        str_isAllCHPEDallele = \
            "- Are all alleles CHPED allele?: {}\n".format(self.isAllCHPEDallele)

        str_Nfield = \
            "- Field format: {}\n".format(FieldFormat2Label(self.field_format))

        str_summary = ''.join([
            str_chped,
            str_hasHeader,
            str_HLA,
            str_isAllCHPEDallele,
            str_Nfield
        ]).rstrip('\n')

        return str_summary


    def __bool__(self):
        return Exists(self.chped)


    def subsetCHPED(self, _out_prefix, _HLA_ToExtract:list, _HLA_ToExclude=()):

        _HLA_ToExtract = list(np.setdiff1d(_HLA_ToExtract, _HLA_ToExclude))

        if not self.f_hasHeader:
            print(std_ERROR + "Given CHPED('{}') can't be subsetted to the HLA genes('{}') "
                              "because the CHPED doesn't have the header.".format(self.chped, _HLA_ToExtract))
            return -1

        _out = _out_prefix + '.chped'

        l_target_columns = \
            ['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe'] + \
            ['{}_{}'.format(hla, i) for hla in _HLA_ToExtract for i in np.arange(1, 3)]

        df_chped_subset = self.df_chped.loc[:, l_target_columns]
        # print(df_chped_subset)

        df_chped_subset.to_csv(_out, sep='\t', header=True, index=False)
        return CHPED(_out)



def getAvailableHLA(_df_chped_Right) -> list:

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
            print(std_WARNING +
                  "There are more than two HLA genes('{}') for {}-th and {}-th columns" \
                  .format(arr_HLA, idx1, idx2))
        else:
            l_RETURN.append(arr_HLA[0])

    return l_RETURN


isCHPEDallele = lambda x : bool(re.match(r'\w+\*\d{2,3}(:\d{2,3})*[A-Z]?', x)) if x != '0' else True
guessFieldFormat = lambda x : len(x.split(':')) if x != '0' else 0



if __name__ == '__main__':

    """
    Wrapper class for the CHPED file('*.chped').

    """

    # Header
    # _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/Dummy.N50.imgt3460.3major.header.chped"
    # _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.8major.header.chped"
    _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.chped"

    # no Header
    # _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.8major.noheader.chped"
    # _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/Dummy.N50.imgt3460.3major.noheader.chped"
    # _chped = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.noheader.chped"

    CH = CHPED(_chped)
    print(CH)

    CH2 = CH.subsetCHPED("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220412/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset",
                         _HLA_ToExtract=['A', 'B', 'C'], _HLA_ToExclude=['DRB1'])
    print(CH2)

    # df_chped = pd.read_csv(_chped, sep='\s+', header=None, dtype=str)
    # HLA_avail = getAvailableHLA(df_chped.iloc[:, 6:])
    # print(HLA_avail)