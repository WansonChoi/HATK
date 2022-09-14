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



class PHENO(object):
    """
    A wrapper class to manage PLINK phenotype files.
    """
    def __init__(self, _file, _pheno_name:list=()):

        ### Main Variables ###
        self.file = checkFile(_file, std_ERROR + "Given PLINK Phenotype file('{}') can't be found.".format(_file))
        self.pheno_name_req = _pheno_name
        self.pheno_name_avail = []
        self.pheno_name_target = []

        self.f_hasHeader = hasHeader(self.file)

        self.N_pheno_name_target = 0
        self.DF_pheno = None
        self.N_samples = None


        ### Main Actions ###
        # Investigate given PHENO file. (depending on header) - `self.DF_pheno`, `self.pheno_name_avail`, `self.N_samples`.
        if self.f_hasHeader:
            self.DF_pheno = pd.read_csv(self.file, sep='\s+', header=0, dtype=str)
            self.pheno_name_avail = self.DF_pheno.columns.to_list()[2:]
            self.N_samples = self.DF_pheno.shape[0]
        else:
            self.DF_pheno = pd.read_csv(self.file, sep='\s+', header=None, dtype=str)
            if self.DF_pheno.shape[1] > 3: print(std_WARNING + "There are more than 2 phenotypes without header.")
            self.N_samples = self.DF_pheno.shape[0]

        # get Target Phenotype(s). - `self.pheno_name_target`, `self.N_pheno_name_target`
        if len(self.pheno_name_req) > 0:
            self.pheno_name_target = list(np.intersect1d(self.pheno_name_req, self.pheno_name_avail))
            self.N_pheno_name_target = len(self.pheno_name_target)

            # Excluded ones.
            l_excluded = list(np.setdiff1d(self.pheno_name_req, self.pheno_name_avail))
            if len(l_excluded) > 0:
                print(std_WARNING + "Next phenotypes are NOT IN given phenotype file: {}".format(l_excluded))


    def __repr__(self):
        # str_main = "< PLINK Phenotype file summary >\n"
        str_file = \
            "- Phenotype file: {}\n".format(self.file)
        str_pheno_name_req = \
            "- Phenotypes requested: {}\n".format(self.pheno_name_req)
        str_pheno_name_avail = \
            "- Phenotypes available: {}\n".format(self.pheno_name_avail)
        str_pheno_name_target = \
            "- Phenotypes target: {}\n".format(self.pheno_name_target)
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
            self.DF_pheno.loc[:, self.pheno_name_target].to_csv(_out, sep='\t', header=True, index=False)
            return _out
        else:
            print(std_WARNING + "No header line to subset.")
            return -1



if __name__ == '__main__':

    ### Phenotype class
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe')
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD'])
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD', 'All_CD'])
    # print(pt)
    pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD', 'All_CD', 'asdf'])
    print(pt)