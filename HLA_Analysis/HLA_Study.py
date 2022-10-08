#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.util import findExec
from src.Study import Study
from bMarkerGenerator.bMarker import bMarker

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class HLA_Study(Study):

    """
    A container class for HLA fine-mapping analysis.
    """

    ### External Software ###
    beagle = findExec("beagle", std_ERROR+"'beagle' command can't be found. Please install it.")
    java = findExec("java", std_ERROR+"'java' command can't be found. Please install it.")

    beagle2vcf = findExec("beagle2vcf.jar",
        std_ERROR+"'beagle2vcf.jar' module can't be found.\n"
                  "Please (1) Download that jar file from 'https://faculty.washington.edu/browning/beagle_utilities/beagle2vcf.jar' AND \n"
                  "(2) Locate it in the 'dependency/' folder.")
    vcf2beagle = findExec("vcf2beagle.jar",
        std_ERROR+"'vcf2beagle.jar' module can't be found.\n"
                  "Please (1) Download that jar file from 'https://faculty.washington.edu/browning/beagle_utilities/vcf2beagle.jar' AND \n"
                  "(2) Locate it in the 'dependency/' folder.")
    linkage2beagle = findExec("linkage2beagle.jar",
        std_ERROR+"'linkage2beagle.jar' module can't be found.\n"
                  "Please (1) Download that jar file from 'https://faculty.washington.edu/browning/beagle_utilities/linkage2beagle.jar' AND \n"
                  "(2) Locate it in the 'dependency/' folder.")


    def __init__(self, _out_prefix, _bfile, _pheno, _pheno_name, _covar=None, _covar_name=None, _condition_list=None):

        ### Init. of super class ###
        super().__init__(_out_prefix, _bfile, _pheno, _pheno_name, _covar, _covar_name, _condition_list)


        ### Main Variables ###
        self.bmarker = bMarker(_bfile) # binary markers for HLA fine-mapping.
        self.bgl_phased = None

        # results
        self.assoc = None
        self.omnibus = None

        self.manhattan = None
        self.Heatmap = None



    def __repr__(self):
        str_GT = \
            "=====< Binary Markers(Genotype) >=====\n{}\n".format(self.bmarker)
        str_PT = \
            "=====< PHENOTYPE >=====\n{}\n".format(self.PHENO)
        str_pheno_name = \
            "\n- Target Phenotype: {}\n".format(self.pheno_name)
        str_pheno_dtype_target = \
            "- Target Phenotype Dtype: {}\n".format(self.pheno_name_dtype)
        str_CV = \
            "=====< COVARIATE >=====\n{}\n".format(self.COVAR) if self.COVAR else ""
        str_covar_name = \
            "\n- Target Covariate(s): {}\n".format(self.covar_name) if self.COVAR else ""
        str_condition = \
            "=====< CONDITION >=====\n{}\n".format(self.CONDITION) if self.CONDITION else ""

        str_external_soft = \
            "=====< External Software >=====\n" \
            "- PLINK: {}\n" \
            "- R: {}\n" \
            "- java: {}\n" \
            "- Beagle: {}\n" \
            "- beagle2vcf: {}\n" \
            "- vcf2beagle: {}\n" \
            "- linkage2beagle: {}\n" \
            .format(self.plink, self.Rscript, self.java, self.beagle, self.beagle2vcf, self.vcf2beagle,
                    self.linkage2beagle)

        str_summary = \
            ''.join([str_GT, str_PT, str_pheno_name, str_pheno_dtype_target,
                     str_CV, str_covar_name, str_condition,
                     str_external_soft]).rstrip('\n')
        return str_summary



if __name__ == '__main__':

    ## Linux

    ## Mac
    OUT = "/Users/wansonchoi/Git_Projects/HATK/tests/20221007_bMG/merged.STUDY"
    _bfile = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged"
    PT = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.phe"
    pheno_name = "Dis_CD"
    CV = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.covar"
    covar_name = "Immunochip2,GWAS"

    r = HLA_Study(OUT, _bfile, PT, pheno_name, CV, covar_name)
    print(r)


    pass