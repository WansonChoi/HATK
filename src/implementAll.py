#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join
import time

from IMGT2Seq.__main__ import HATK_IMGT2Seq
from NomenCleaner.__main__ import HATK_NomenCleaner
from bMarkerGenerator.__main__ import HATK_bMarkerGenertor
from HLA_Analysis.HLA_Study import HLA_Study

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)

class implementAll(object):
    """
    An action class for whole implementation.
    - IMGT2Seq
    - NomenCleaner
    - bMarkerGenerator
    - Association test
    - Phasing
    - Ominbus test
    - Manhattan for assoc
    - Heatmap for assoc
    - Manhattan for ominbus test.

    HLA gene set will be set as the HLA genes in hped files.
    """

    def __init__(self, _out_prefix, _hg, _imgt, _imgt_dir, _hped, _bfile, _pheno, _pheno_name, _pheno_name_dtype=None,
                _covar=None, _covar_name=None, _condition=None, _condition_list=None,
                _F_one:bool=False, _F_two:bool=False, _F_three:bool=False, _F_four:bool=False, _F_Ggroup:bool=False, _F_Pgroup:bool=False,
                _MultiP=1, _java_mem='1G', _nthreads=1,
                _f_save_intermediates=False):

        _out_dir = dirname(_out_prefix)

        ### =====< IMGT2Seq >=====
        self.imgt_output = \
            HATK_IMGT2Seq(_imgt, _hg, _out_prefix, _imgt_dir,
                          _F_one, _F_two, _F_three, _F_four, _F_Ggroup, _F_Pgroup,
                          _MultiP=_MultiP, _f_save_intermediates=_f_save_intermediates)

        # variables to be used next.
        HAT = self.imgt_output.IMGT2Seq_Output.HAT # => NC
        dict_AA = self.imgt_output.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_AA_seq.rstrip('.txt')
        dict_SNPS = self.imgt_output.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_SNPS_seq.rstrip('.txt')

        # self.imgt_output = "imgt"
        # HAT = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATKv2_wholeImple_20221012/HLA_ALLELE_TABLE.imgt3470.hat"
        # dict_AA = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATKv2_wholeImple_20221012/HLA_DICTIONARY_AA.hg18.imgt3470"
        # dict_SNPS = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATKv2_wholeImple_20221012/HLA_DICTIONARY_SNPS.hg18.imgt3470"


        ### =====< NomenCleaner >=====
        self.NC_output = \
            HATK_NomenCleaner(
                _hped, HAT, _out_prefix,
                _F_one, _F_two, _F_three, _F_four, _F_Ggroup, _F_Pgroup,
            )

        chped = self.NC_output.chped


        ### =====< bMarkerGenerator >=====
        self.bMG_output = \
            HATK_bMarkerGenertor(
                chped, _out_prefix, _hg, dict_AA, dict_SNPS, _bfile=_bfile,
                _f_save_intermediates=_f_save_intermediates
            )

        bMARKER_prefix = self.bMG_output.bmarker.file_prefix


        self.hStudy = HLA_Study(_out_prefix, _hg, bMARKER_prefix, _pheno, _pheno_name, _pheno_name_dtype,
                                _covar, _covar_name, _condition_list,
                                _f_save_intermediates, _imgt_out_dir=_out_dir, _java_mem=_java_mem, _nthreads=_nthreads)

        ### =====< Associaion test >=====
        self.hStudy.doRegression()

        ### =====< Manhattan plot - Assoc >=====
        self.hStudy.doManhattanPlot_assoc()

        ### =====< Heatmap - Assoc >=====
        self.hStudy.doHeatmapPlot()

        ### =====< Phasing >=====
        self.hStudy.doPhasing()

        ### =====< Omnibus test >=====
        self.hStudy.doOmnibusTest()

        ### =====< Manhattan plot - Omnibus >=====
        self.hStudy.doManhattanPlot_omnibus()




    def __repr__(self):
        str_imgt_output = \
            "=====< IMGT2Seq result >=====\n{}\n".format(self.imgt_output)

        str_NC_output = \
            "=====< NomenCleaner result >=====\n{}\n".format(self.NC_output)

        str_bMG_output = \
            "=====< bMarkerGenerator result >=====\n{}\n".format(self.bMG_output)

        str_hStudy_output = \
            "=====< bMarkerGenerator result >=====\n{}\n".format(self.hStudy)


        str_summary = ''.join([
            str_imgt_output, str_NC_output, str_bMG_output

        ]).rstrip("\n")
        return str_summary