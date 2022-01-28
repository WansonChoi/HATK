#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists
import pandas as pd
import numpy as np
import src.HATK_Error as HATK_Error

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)



class PLINK(object):
    """
    Wrapper class to manage input PLINK files for association test study.

    Subclasses to be composed. (has-a relationship)
    - Genotype
    - Phenetype
    - Covariate
    - AssociationTest

    """
    def __init__(self, _file_prefix_forAll=-1, _file_GT=-1, _file_PT=-1, _file_CV=-1):
        ### Main Variables ###
        self.GT = -1
        self.PT = -1
        self.CV = -1

        ### Main Actions ###
        self.setPlinkFiles(_file_prefix_forAll, _file_GT, _file_PT, _file_CV)
        # self.checkPlinkFiles()


    def setPlinkFiles(self, _file_prefix_forAll, _file_GT, _file_PT, _file_CV):
        if _file_prefix_forAll != -1:
            _file_PT = _file_prefix_forAll+'.phe' if exists(_file_prefix_forAll+'.phe') else \
                        _file_prefix_forAll +'.pheno' if exists(_file_prefix_forAll+'.pheno') else -1
            _file_CV = _file_prefix_forAll+'.cov' if exists(_file_prefix_forAll+'.cov') else \
                        _file_prefix_forAll + '.covar' if exists(_file_prefix_forAll+'.covar') else -1

        # set instances
        self.GT = Genotype(_file_GT) if _file_GT != -1 else -1
        self.PT = Phenotype(_file_PT) if _file_PT != -1 else -1
        self.CV = Covariate(_file_CV) if _file_CV != -1 else -1


    def checkPlinkFiles(self):
        # ex) Sample number in the files are consistent or not.
        pass


    def __repr__(self):
        str_GT = "=====< GENOTYPE >=====\n{}\n".format(self.GT)
        str_PT = "=====< PHENOTYPE >=====\n{}\n".format(self.PT)
        str_CV = "=====< COVARIATE >=====\n{}\n".format(self.CV)
        str_endline = "=======================\n"

        str_summary = ''.join([str_GT, str_PT, str_CV, str_endline]).rstrip('\n')
        return str_summary



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
            std_ERROR_MAIN_PROCESS_NAME +
            "The '*.bed' file('{}') can't be found.".format(_file_prefix+'.bed'))

        # (2) '*.bim'
        if exists(_file_prefix+'.bim'):
            self.bim = _file_prefix+'.bim'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "The '*.bim' file('{}') can't be found.".format(_file_prefix + '.bim'))

        # (3) '*.fam'
        if exists(_file_prefix+'.fam'):
            self.fam = _file_prefix+'.fam'
        else: raise HATK_Error.HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
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
            std_ERROR_MAIN_PROCESS_NAME +
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
                std_ERROR_MAIN_PROCESS_NAME +
                "The phenotype file('{}') must have the header line starting with 'FID', 'IID', "
                "and phenotype name(s).".format(self.file_phe))


    def __repr__(self):
        str_phe = "Phenotype file: {}\n".format(self.file_phe)
        str_info = "N(samples): {}\nhasHeader: {}\nPhenotypes: {} ({} item(s))\n" \
                    .format(self.N_samples, self.f_hasHeader, self.l_phenotypes, self.N_phenotypes)

        str_summary = ''.join([str_phe, str_info]).rstrip('\n')
        return str_summary


    def __bool__(self):
        return self.file_phe != -1



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
            std_ERROR_MAIN_PROCESS_NAME +
            "The given Covariate file('{}') can't be found.".format(_file))


    def checkCovariates(self):

        def hasHeader():

            p_Header = re.compile(r'^FID\s+IID\s+')

            with open(self.file_covar, 'r') as f_covar:
                line_1st = f_covar.readline()

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
                std_ERROR_MAIN_PROCESS_NAME +
                "The covariate file('{}') must have the header line starting with 'FID', 'IID', "
                "and covariate name(s).".format(self.file_covar))


    def __repr__(self):
        str_file_covar = "Covariate file: {}\n".format(self.file_covar)
        str_info = "N(samples): {}\nhasHeader: {}\nCovariates: {} ({} item(s))\n" \
                    .format(self.N_samples, self.f_hasHeader, self.l_covariates, self.N_covariates)

        str_summary = ''.join([str_file_covar, str_info]).rstrip('\n')
        return str_summary


    def __bool__(self):
        return self.file_covar != -1



class AssociationTest(object):
    """
    A wrapper class to manage PLINK association test result.
    """
    def __init__(self, _file):
        ### Main Variables ###
        self.assoc = -1

        ### Main Actions ###
        self.setAssoc(_file)


    def setAssoc(self, _file):
        if exists(_file):
            self.assoc = _file
        else: raise HATK_Error.HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "The given association test result can't be found('{}').".format(_file))


    def __repr__(self):
        str_assoc = "Association Test Result: '{}'\n".format(self.assoc)

        str_summary = ''.join([str_assoc])
        return str_summary


# ==================== ==================== ==================== #

if __name__ == '__main__':

    ### Genotype class
    # gt = Genotype('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged')
    # print(gt)

    ### Phenotype class
    # pt = Phenotype('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe')
    # print(pt)

    ### Covariate class
    # cv = Covariate('/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar')
    # print(cv)

    ### AssociationTest class
    # assoc = AssociationTest('/home/wansonchoi/sf_VirtualBox_Share/UC-CD-HLA/analysis/01-association/All_CD.assoc.logistic')
    # print(assoc)

    ### PLINK class
    # myplink = PLINK(_file_GT='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged',
    #                 _file_PT='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe',
    #                 _file_CV='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar')
    # print(myplink)

    pass