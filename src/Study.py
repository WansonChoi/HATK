#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.PLINK_GT import GT
from src.PLINK_PHENO import PHENO, getTargetPheDtype
from src.PLINK_COVAR import COVAR
from src.PLINK_CONDITION import CONDITION
from src.PLINK_Bash import getPLINKexec
from src.util import findExec

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class Study(object):

    """
    A Container class for Association test.

    """

    ### External Software ###
    plink = getPLINKexec()
    Rscript = findExec("Rscript", std_ERROR+"'Rscript' command can't be found. Please install R.")


    def __init__(self, _out_prefix, _hg, _bfile, _pheno, _pheno_name, _pheno_name_dtype=None, _covar=None, _covar_name=None, _condition_list=None,
                 _f_save_intermediates=False):

        ### Main Variables ###
        self.hg = str(_hg)
        self.f_save_intermediates = _f_save_intermediates

        # (1) Genotype (required)
        self.GT = GT(_bfile)

        # (2) Phenotype (required)
        self.PHENO = PHENO(_pheno)
        self.pheno_name = _pheno_name # (***) must be A SINGLE phenotype name. (1 trait - 1 study instance)
        if self.pheno_name not in self.PHENO.pheno_name_avail: # Major exception handling.
            RaiseError(HATK_InputPreparation_Error,
                       std_ERROR + "Requested phenotype name('{}') is NOT IN the given phenotype file('{}')." \
                       .format(self.pheno_name, self.PHENO.pheno_name_avail))
        self.pheno_name_dtype = _pheno_name_dtype if _pheno_name_dtype else \
                                getTargetPheDtype(self.PHENO.pheno_name_avail, self.PHENO.trait_types, self.pheno_name)

        # (3) Covariate (optional)
        self.COVAR = COVAR(_covar, _covar_name.split(',') if isinstance(_covar_name, str) else [])
        # (4) Condition (optional)
        self.CONDITION = CONDITION(None, _condition_list)

        self.out_prefix = _out_prefix


    def __repr__(self):
        str_hg = \
            "Human Genome Build: hg{}\n".format(self.hg)
        str_out_prefix = \
            "Output prefix: {}\n".format(self.out_prefix)
        str_GT = \
            "=====< GENOTYPE >=====\n{}\n".format(self.GT)
        str_PT = \
            "=====< PHENOTYPE >=====\n{}\n".format(self.PHENO)
        str_pheno_name_target = \
            "\n- Target Phenotype: {}\n".format(self.pheno_name)
        str_pheno_dtype_target = \
            "- Target Phenotype Dtype: {}\n".format(self.pheno_name_dtype)
        str_CV = \
            "=====< COVARIATE >=====\n{}\n".format(self.COVAR) if self.COVAR else ""
        str_covar_name = \
            "\n- Target Covariate(s): {}\n".format(self.COVAR.covar_name_target) if self.COVAR else ""
        str_condition = \
            "=====< CONDITION >=====\n{}\n".format(self.CONDITION) if self.CONDITION else ""

        str_external_soft = \
            "=====< External Software >=====\n" \
            "- PLINK: {}\n" \
            "- R: {}".format(self.plink, self.Rscript)

        str_summary = \
            ''.join([str_hg, str_out_prefix, str_GT, str_PT, str_pheno_name_target, str_pheno_dtype_target,
                     str_CV, str_covar_name, str_condition,
                     str_external_soft]).rstrip('\n')
        return str_summary



if __name__ == '__main__':

    ## LINUX
    # OUT = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_Study_20220914/wtccc_filtered_58C_RA.hatk.300+300.hg18.chr6.29-34mb.STUDY"
    # GT = Genotype("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_bMG_20220913/wtccc_filtered_58C_RA.hatk.300+300.hg18.chr6.29-34mb")
    # print(GT)
    # PT = PHENO("/home/wansonchoi/sf_VirtualBox_Share/HATK/example/wtccc_filtered_58C_RA.hatk.300+300.phe")
    # print(PT)

    # OUT = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_Study_20220914/merged.STUDY"
    # GT = GT("/home/wansonchoi/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged")
    # PT = PHENO("/home/wansonchoi/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe", ["Dis_CD"])

    ## Mac
    # OUT = "/Users/wansonchoi/Git_Projects/HATK/tests/20221007_bMG/merged.STUDY"
    # _bfile = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged"
    # PT = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.phe"
    # pheno_name = "Dis_CD"

    OUT = "/Users/wansonchoi/Git_Projects/HATK/tests/20221007_bMG/merged.STUDY"
    _bfile = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged"
    PT = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.phe"
    pheno_name = "Dis_UC"
    CV = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.covar"
    covar_name = "Immunochip2,GWAS"
    hg = 18

    r = Study(OUT, hg, _bfile, PT, pheno_name, CV, covar_name)
    print(r)

    pass



#     def plotManhattan(self, _out_prefix):
#
#         # Temporarily interfaced with this way.
#         m_plot = Manhattan([self.assoc], self.assoc, '18', _p_src='../HLA_Manhattan/src',
#                            _p_data='../HLA_Manhattan/data')
#         return m_plot


        # def doStepwiseConditional(self):
        #     # 'doAssociationTest()' and 'plotManhattan()' will be used, here.
        #     pass