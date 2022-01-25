#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists
import pandas as pd
import numpy as np
import src.HATK_Error as HATK_Error

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)

print_MAIN = lambda x : print(std_MAIN_PROCESS_NAME + x)
print_ERROR = lambda x : print(std_ERROR_MAIN_PROCESS_NAME + x)
print_WARNING = lambda x : print(std_WARNING_MAIN_PROCESS_NAME + x)



class PLINK(object):
    """
    Wrapper class to manage input PLINK files for association test.

    Subclasses to be composed. (has-a relationship)
    - Genotype
    - Phenetype
    - Covariate
    - AssociationTest

    """
    def __init__(self):

        Phenotype('asdf')

        pass



class Genotype(object):
    """
    A wrapper class to manage PLINK genotype files(*.{bed,bim,fam}).
    """
    def __init__(self, _file_prefix):

        ### Main Variables ###
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
            "The '*.bed' file('{}') can't be found.".format(_file_prefix+'.bed'), print_ERROR)

        # (2) '*.bim'
        if exists(_file_prefix+'.bim'):
            self.bim = _file_prefix+'.bim'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            "The '*.bim' file('{}') can't be found.".format(_file_prefix + '.bim'), print_ERROR)

        # (3) '*.fam'
        if exists(_file_prefix+'.fam'):
            self.fam = _file_prefix+'.fam'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            "The '*.fam' file('{}') can't be found.".format(_file_prefix + '.fam'), print_ERROR)


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

        str_summary = ''.join([str_files, str_shape, str_Info])

        return str_summary



class Phenotype(object):
    """
    A wrapper class to manage PLINK phenotype files.
    """
    def __init__(self, _file):

        ### Main Variables ###
        self.file_phe = -1

        self.f_hasHeader = -1
        self.N_phenotypes = -1
        self.l_phenotypes = -1

        self.N_samples = -1

        ### Main Actions ###
        self.setPhenotype(_file)
        self.checkPhenotype()


    def setPhenotype(self, _file):
        if exists(_file):
            self.file_phe = _file
        else: raise HATK_Error.HATK_InputPreparation_Error(
            "The given Phenotype file('{}') can't be found.".format(_file))


    def checkPhenotype(self):

        def hasHeader():

            p_Header = re.compile(r'^FID\s+IID\s+')
            p_Line = re.compile(r'^(\S+)\s+(\S+)\s+(-?\d+)')

            with open(self.file_phe, 'r') as f_phe:
                line_1st = f_phe.readline()

                m_Header = p_Header.match(line_1st)
                # print(m_Header)
                m_Line = p_Line.match(line_1st)
                # print(m_Line)

            return bool(m_Header) or (not bool(m_Line))

        f_hasHeader = hasHeader()

        if f_hasHeader:
            df_phe = pd.read_csv(self.file_phe, sep='\s+', header=0, dtype=str)
            # print("df_phe:\n{}\n".format(df_phe))

            l_phenotypes = df_phe.columns.tolist()[2:]

            ### Set summary values.
            self.f_hasHeader = f_hasHeader
            self.N_phenotypes = len(l_phenotypes)
            self.l_phenotypes = l_phenotypes

            self.N_samples = df_phe.shape[0]

        else:
            raise HATK_Error.HATK_InputPreparation_Error(
                "The phenotype file('{}') must have the header line starting with 'FID', 'IID', "
                "and phenotype name(s).".format(self.file_phe))


    def __repr__(self):
        str_phe = "Phenotype file: {}\n".format(self.file_phe)
        str_info = "N(samples): {}\nhasHeader: {}\nPhenotypes: {} ({} item(s))\n" \
                    .format(self.N_samples, self.f_hasHeader, self.l_phenotypes, self.N_phenotypes)

        str_summary = ''.join([str_phe, str_info])
        return str_summary



class Covariate(object):
    """
    A wrapper class to manage PLINK covariate files.
    """
    def __init__(self, _file):
        ### Main Variables ###
        self.file_covar = -1

        self.f_hasHeader = -1
        self.N_covariates = -1
        self.l_covariates = -1

        self.N_samples = -1

        ### Main Actions ###
        self.setCovariates(_file)
        self.checkCovariates()


    def setCovariates(self, _file):
        if exists(_file):
            self.file_covar = _file
        else: raise HATK_Error.HATK_InputPreparation_Error(
            "The given Covariate file('{}') can't be found.".format(_file))


    def checkCovariates(self):

        def hasHeader():

            p_Header = re.compile(r'^FID\s+IID\s+')

            with open(self.file_covar, 'r') as f_phe:
                line_1st = f_phe.readline()

                m_Header = p_Header.match(line_1st)
                # print(m_Header)

            return bool(m_Header)

        f_hasHeader = hasHeader()

        if f_hasHeader:
            df_covar = pd.read_csv(self.file_covar, sep='\s+', header=0)

            l_covariates = df_covar.columns.tolist()[2:]

            ### Set summary values.
            self.f_hasHeader = f_hasHeader
            self.N_covariates = len(l_covariates)
            self.l_covariates = l_covariates

            self.N_samples = df_covar.shape[0]

        else:
            raise HATK_Error.HATK_InputPreparation_Error(
                "The covariate file('{}') must have the header line starting with 'FID', 'IID', "
                "and covariate name(s).".format(self.file_covar))


    def __repr__(self):
        str_file_covar = "Covariate file: {}\n".format(self.file_covar)
        str_info = "N(samples): {}\nhasHeader: {}\nCovariates: {} ({} item(s))\n" \
                    .format(self.N_samples, self.f_hasHeader, self.l_covariates, self.N_covariates)

        str_summary = ''.join([str_file_covar, str_info])
        return str_summary



class AssociationTest(object):
    """
    A wrapper class to manage PLINK association test result.
    """
    def __init__(self):
        pass

# ==================== ==================== ==================== #

if __name__ == '__main__':

    ### Genotype class
    # gt = Genotype('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged')
    # print(gt)

    ### Phenotype class
    # pt = Phenotype('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe')
    # print(pt)

    ### Covariate class
    cv = Covariate('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar')
    print(cv)

    pass