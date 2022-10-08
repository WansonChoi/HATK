#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join

import src.HATK_Error as HATK_Error
from src.PLINK_GT import GT
from src.PLINK_PHENO import PHENO, getTargetPheDtype
from src.PLINK_COVAR import COVAR
from src.PLINK_CONDITION import CONDITION

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class Study(object):

    """
    A Container class for Association test.

    """

    # def __init__(self, _out_prefix, _GT:GT, _pheno_name, _PT:PHENO=None, _CV:COVAR=None, _CD:CONDITION=None):
    def __init__(self, _out_prefix, _bfile, _pheno, _pheno_name, _covar=None, _covar_name=None, _condition_list=None):

        ### Main Variables ###
        self.GT = GT(_bfile)

        self.PHENO = PHENO(_pheno)
        self.pheno_name = _pheno_name # (***) must be A SINGLE phenotype name. (1 trait - 1 study instance)
        self.pheno_name_dtype = getTargetPheDtype(self.PHENO.pheno_name_avail, self.PHENO.trait_types, self.pheno_name)

        self.COVAR = COVAR(_covar, _covar_name.split(',')) if isinstance(_covar, str) else None
        self.covar_name = _covar_name.split(',')

        self.CONDITION = CONDITION(None, _condition_list) if isinstance(_condition_list, str) else None

        self.out_prefix = _out_prefix


    def __repr__(self):
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
            "\n- Target Covariate(s): {}\n".format(self.covar_name) if self.COVAR else ""
        str_condition = \
            "=====< CONDITION >=====\n{}\n".format(self.CONDITION) if self.CONDITION else ""

        str_summary = \
            ''.join([str_GT, str_PT, str_pheno_name_target, str_pheno_dtype_target,
                     str_CV, str_covar_name, str_condition]).rstrip('\n')
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

    r = Study(OUT, _bfile, PT, pheno_name, CV, covar_name)
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