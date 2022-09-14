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



class COVAR(object):
    """
    A wrapper class to manage PLINK covariate files.
    """
    def __init__(self, _file, _covar_name:list=()):
        ### Main Variables ###
        self.file = checkFile(_file, std_ERROR + "Given PLINK Covariate file('{}') can't be found.".format(_file))

        self.covar_name_req = _covar_name
        self.covar_name_avail = []
        self.covar_name_target = []

        self.f_hasHeader = hasHeader(self.file)

        self.N_covar_name_target = 0
        self.DF_covar = None
        self.N_samples = None

        ### Main Actions ###
        # Investigate given COVAR file. (depending on header) - `self.DF_covar`, `self.covar_name_avail`, `self.N_samples`.
        if self.f_hasHeader:
            self.DF_covar = pd.read_csv(self.file, sep='\s+', header=0, dtype=str)
            self.covar_name_avail = self.DF_covar.columns.to_list()[2:]
            self.N_samples = self.DF_covar.shape[0]
        else:
            self.DF_covar = pd.read_csv(self.file, sep='\s+', header=None, dtype=str)
            if self.DF_covar.shape[1] > 3: print(std_WARNING + "There are more than 2 covariates without header.")
            self.N_samples = self.DF_covar.shape[0]

        # get Target Phenotype(s). - `self.covar_name_target`, `self.N_covar_name_target`
        if len(self.covar_name_req) > 0:
            self.covar_name_target = list(np.intersect1d(self.covar_name_req, self.covar_name_avail))
            self.N_covar_name_target = len(self.covar_name_target)

            # Excluded ones.
            l_excluded = list(np.setdiff1d(self.covar_name_req, self.covar_name_avail))
            if len(l_excluded) > 0:
                print(std_WARNING + "Next covariates are NOT IN given covariate file: {}".format(l_excluded))


    def __repr__(self):
        # str_main = "< PLINK Covariate file summary >\n"
        str_file = \
            "- Covariate file: {}\n".format(self.file)
        str_pheno_name_req = \
            "- Covariates requested: {}\n".format(self.covar_name_req)
        str_pheno_name_avail = \
            "- Covariates available: {}\n".format(self.covar_name_avail)
        str_pheno_name_target = \
            "- Covariates target: {}\n".format(self.covar_name_target)
        str_hasHeader = \
            "- has Header?: {}\n".format(self.f_hasHeader)
        str_N_samples = \
            "- # of samples: {}\n".format(self.N_samples)


        str_summary = ''.join([str_file,
                               str_pheno_name_req, str_pheno_name_avail, str_pheno_name_target,
                               str_hasHeader, str_N_samples]).rstrip('\n')
        return str_summary


    def __bool__(self): return Exists(self.file)


    def writeSubset(self, _out):
        if self.f_hasHeader:
            self.DF_covar.loc[:, self.covar_name_target].to_csv(_out, sep='\t', header=True, index=False)
            return _out
        else:
            print(std_WARNING + "No header line to subset.")
            return -1



if __name__ == '__main__':

    ### Covariate class
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar')
    # print(cv)
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar', ['Immunochip2'])
    # print(cv)
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar', ['Immunochip2', 'GWAS'])
    # print(cv)
    cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar', ['Immunochip2', 'GWAS', 'asdf'])
    print(cv)
