#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists
import pandas as pd
import numpy as np

import src.HATK_Error as HATK_Error
from src.util import Exists, checkFile, hasHeader, getHeader

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class Genotype(object):
    """
    A wrapper class to manage PLINK genotype files(*.{bed,bim,fam}).
    """
    def __init__(self, _file_prefix):

        ### Main Variables ###
        self.file_gt = -1
        self.bed = -1
        self.bim = -1
        self.fam = -1

        self.N_samples = -1
        self.M_markers = -1

        self.f_hasSexInfo = -1
        self.f_hasPheInfo = -1

        ### Main Actions ###
        self.setGenotype(_file_prefix)
        self.checkFAM()
        self.checkBIM()

    # ========== ========== #

    def setGenotype(self, _file_prefix):

        # (1) '*.bed'
        if exists(_file_prefix+'.bed'):
            self.bed = _file_prefix+'.bed'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            std_ERROR +
            "The '*.bed' file('{}') can't be found.".format(_file_prefix+'.bed'))

        # (2) '*.bim'
        if exists(_file_prefix+'.bim'):
            self.bim = _file_prefix+'.bim'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            std_ERROR +
            "The '*.bim' file('{}') can't be found.".format(_file_prefix + '.bim'))

        # (3) '*.fam'
        if exists(_file_prefix+'.fam'):
            self.fam = _file_prefix+'.fam'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            std_ERROR +
            "The '*.fam' file('{}') can't be found.".format(_file_prefix + '.fam'))

        self.file_gt = _file_prefix


    def checkFAM(self):

        df_fam = pd.read_csv(self.fam, sep='\s+', header=None, dtype=str)
        # print("df_fam:\n{}\n".format(df_fam))

        arr_sex = df_fam.iloc[:, 4].values
        f_hasSexInfo = not np.all( np.logical_or((arr_sex == '-9'), (arr_sex == '0')) )
        # print(f_hasSexInfo)

        arr_phe = df_fam.iloc[:, 5].values
        f_hasPheInfo = not np.all( np.logical_or((arr_phe == '-9'), (arr_phe == '0')) )
        # print(f_hasPheInfo)

        ### Set summary values.
        self.N_samples = df_fam.shape[0]
        self.f_hasSexInfo = f_hasSexInfo
        self.f_hasPheInfo = f_hasPheInfo


    def checkBIM(self):

        df_bim = pd.read_csv(self.bim, sep='\s+', header=None, dtype=str)

        ### Set summary values.
        self.M_markers = df_bim.shape[0]


    def __repr__(self):

        str_files = "bed: {}\nbim: {}\nfam:{}\n".format(self.bed, self.bim, self.fam)
        str_shape = "N(Samples): {}\nM(Markers): {}\n".format(self.N_samples, self.M_markers)
        str_Info = "hasSexInfo: {}\nhasPheInfo: {}\n".format(self.f_hasSexInfo, self.f_hasPheInfo)

        str_summary = ''.join([str_files, str_shape, str_Info]).rstrip('\n')

        return str_summary



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



class CONDITION(object):

    def __init__(self, _condition:str=None, _condition_list:str=None, _l_markers:list=()):
        ### Main Variables ###
        self.condition = _condition # comma-separated string (ex. "AA_A_11_exon1,rs1234")
        self.condition_list = _condition_list # file path.

        self.l_condition = []
        self.N_condition = 0


        ### Main Actions ###
        # Condition itme check.
        if self.condition:
            self.l_condition = self.condition.split(',')
        else:
            with open(self.condition_list, 'r') as f_condlist:
                for line in f_condlist:
                    self.l_condition.append(line.rstrip('\n'))
        self.N_condition = len(self.l_condition)

        # All given conditions are in target marekrs?
        if len(_l_markers) > 0:
            l_excluded = list(np.setdiff1d(self.l_condition, _l_markers)) # supposed to [].
            if len(l_excluded) > 0:
                print(std_WARNING + "Some conditions are not in the target marker set: {}".format(l_excluded))


    def __bool__(self): return bool(self.condition) != bool(self.condition_list) # xor
    def isGivenAsFile(self): return (not self.condition) and self.condition_list

    def __repr__(self):
        # str_main = "< PLINK Condition file summary >\n"
        str_file = \
            "- Condition file: {}\n".format(self.condition_list) if self.isGivenAsFile() else \
            "- Condition string: {}\n".format(self.condition)
        str_conditions = \
            "- Conditions: {}\n".format(self.l_condition)
        str_N_conditions = \
            "- # of Conditions: {}\n".format(self.N_condition)

        str_summary = ''.join([str_file, str_conditions, str_N_conditions]).rstrip('\n')
        return str_summary



# ==================== ==================== ==================== #

if __name__ == '__main__':

    ### Genotype class
    # gt = Genotype('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged')
    # print(gt)

    ### Phenotype class
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe')
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD'])
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD', 'All_CD'])
    # print(pt)
    # pt = PHENO('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe', ['Dis_CD', 'All_CD', 'asdf'])
    # print(pt)

    ### Covariate class
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar')
    # print(cv)
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar', ['Immunochip2'])
    # print(cv)
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar', ['Immunochip2', 'GWAS'])
    # print(cv)
    # cv = COVAR('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar', ['Immunochip2', 'GWAS', 'asdf'])
    # print(cv)

    ### Condition class
    cond = CONDITION(_condition_list="/home/wansonchoi/sf_VirtualBox_Share/UC-CD-HLA/analysis/02-conditional-5/ToCond_CD.txt")
    print(cond)

    cond = CONDITION(_condition="AA_DRB1_37_32660037_S,AA_DRB1_37_32660037_N")
    print(cond)


    pass