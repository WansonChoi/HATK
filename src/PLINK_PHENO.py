#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists
import pandas as pd
import numpy as np
from src.util import asSingleList

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.util import Exists, checkFile, hasHeader, getHeader

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class PHENO(object):
    """
    A Container class to manage PLINK phenotype files.
    """

    def __init__(self, _file):

        ### Main Variables ###
        self.phe = checkFile(_file, std_ERROR + "Given PLINK Phenotype file('{}') can't be found.".format(_file))

        self.f_hasHeader = hasHeader(self.phe)
        if not self.f_hasHeader:
            RaiseError(
                HATK_InputPreparation_Error,
                std_ERROR + "Given Phenotype('{}') doesn't have header.".format(_file)
            )

        self.pheno_name_avail, \
        self.N_samples = checkPhe(self.f_hasHeader, self.phe)

        self.trait_types = checkPheDtypes(self.phe, self.pheno_name_avail)



    def __repr__(self):
        str_file = \
            "- Phenotype file: {}\n".format(self.phe)
        str_pheno_name_avail = \
            "- Phenotypes available: {}\n".format(self.pheno_name_avail)
        str_pheno_dtype_target = \
            "- Phenotypes dtype: {}\n".format(self.trait_types)
        str_hasHeader = \
            "- has Header?: {}\n".format(self.f_hasHeader)
        str_N_samples = \
            "- # of samples: {}\n".format(self.N_samples)

        str_summary = \
            ''.join([str_file, str_pheno_name_avail, str_pheno_dtype_target,
                     str_hasHeader, str_N_samples]).rstrip('\n')

        return str_summary



def checkPheDtypes(_phe, _pheno_name_avail:list):

    df_phe = pd.read_csv(_phe, sep='\s+', header=0) # given dtype as it is

    # NA, Int, Float
    isNA = np.vectorize(lambda x : (x == -9) or (x == 0) or (x == -1))
    isInt = np.vectorize(lambda x: isinstance(x, int))
    isFloat = np.vectorize(lambda x: isinstance(x, float))


    l_RETURN = [] # dtypes of available traits.
    for ph in _pheno_name_avail:
        arr_phe = df_phe.loc[:, ph].values

        arr_isNA = isNA(arr_phe)
        arr_isInt = isInt(arr_phe)
        arr_isFloat = isFloat(arr_phe)

        prop_Int = float(sum(arr_isInt)) / len(arr_isInt)
        prop_Float = float(sum(arr_isFloat)) / len(arr_isFloat)

        if np.all(arr_isNA):
            print(std_WARNING + "Given phenotype('{}') values are all NAs(-9, 0, -1).".format(ph))
            l_RETURN.append("NA")
        elif prop_Int > 0.8 and prop_Float < 0.5:
            l_RETURN.append("Binary")
        elif prop_Int < 0.5 and prop_Float > 0.8:
            l_RETURN.append("Continuous")
        else:
            l_RETURN.append("???")

    return l_RETURN



def checkPhe(_f_hasHeader, _phe):

    pheno_name_avail, N_samples = [], 0

    if _f_hasHeader:
        df_phe = pd.read_csv(_phe, sep='\s+', header=0, dtype=str)
        pheno_name_avail = df_phe.columns.to_list()[2:]
        N_samples = df_phe.shape[0]

    else:
        df_phe = pd.read_csv(_phe, sep='\s+', header=None, dtype=str)
        if df_phe.shape[1] > 3: print(std_WARNING + "There are more than 2 phenotypes without header.")
        N_samples = df_phe.shape[0]

    return pheno_name_avail, N_samples


# def writePheSubset(_phe, _out):
#     if self.f_hasHeader:
#         self.DF_pheno.loc[:, self.pheno_name_target].to_csv(_out, sep='\t', header=True, index=False)
#         return _out
#     else:
#         print(std_WARNING + "No header line to subset.")
#         return -1


# def isPheBinary2(_arr_phe, _pheno_name, _pheno_file):
#
#     # All NA values.
#     isNA = np.vectorize(lambda x: (x == -9) or (x == 0) or (x == -1))
#     arr_isNA = isNA(_arr_phe.values)
#     if np.all(arr_isNA):
#         print(std_WARNING + "The given phenotype('{}') values are all NAs(-9, 0, -1).".format(_pheno_name))
#
#     # Mixed ints and floats.
#     isInt = np.vectorize(lambda x: isinstance(x, int))
#     isFloat = np.vectorize(lambda x: isinstance(x, float))
#
#     arr_isInt = isInt(_arr_phe.values)
#     arr_isFloat = isFloat(_arr_phe.values)
#
#     prop_Int = float(sum(arr_isInt)) / len(arr_isInt)
#     prop_Float = float(sum(arr_isFloat)) / len(arr_isFloat)
#
#     if prop_Int > 0.8 and prop_Float < 0.5:
#         return True
#     elif prop_Int < 0.5 and prop_Float > 0.8:
#         print("sadf")
#         return False
#     else:
#         # print("prop_Int: {}(%)".format(prop_Int * 100))
#         # print("prop_Float: {}(%)".format(prop_Float * 100))
#
#         print(
#             std_WARNING.lstrip("\n") +
#             "HATK can't decide whether the target phenotype('{}') is binary(Case-Control) or continuous(Quantitative).\n"
#             "Please check that phenotype column in the given phenotype file('{}')." \
#             .format(_pheno_name, _pheno_file)
#         )
#         print("He")
#         return -1
#
#     """
#     Somthing wrong...
#     """



if __name__ == '__main__':

    ### Phenotype class

    ## Linux
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe')
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD'])
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD', 'All_CD'])
    # print(pt)

    ## Mac
    _pt = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.phe"


    pt = PHENO(_pt)
    print(pt)