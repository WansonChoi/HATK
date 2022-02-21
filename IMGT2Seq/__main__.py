#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import argparse, textwrap
import numpy as np

from IMGT2Seq.IMGT2Seq import IMGT2Seq
from IMGT2Seq.src.AvailableHLAs import getAvailableHLAs, getTargetHLAs
from src.HATK_Error import RaiseError, HATK_InputPreparation_Error
from src.util import Exists

std_MAIN_PROCESS_NAME = "\n[IMGT2Seq]: "
std_ERROR_MAIN_PROCESS_NAME = "\n[IMGT2Seq::ERROR]: "
std_WARNING_MAIN_PROCESS_NAME = "\n[IMGT2Seq::WARNING]: "


class HATK_IMGT2Seq(object):

    def __init__(self, _imgt, _hg, _out,
                 _F_one, _F_two, _F_three, _F_four, _F_Ggroup, _F_Pgroup, _imgt_dir,
                 _no_Indel=False, _MultiP=False, _save_intermediates=False, _no_prime=True, _p_data='./IMGT2Seq/data',
                 _HLA_req=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

        """
        Either (1) generate a new Dictionary files or (2) find existing dictionaries.

        (2022.02.15.) Compare the requested HLAs vs. Available HLAs.
        """

        ### Main Variables ###
        self.imgt = _imgt
        self.hg = _hg
        self.out = _out

        self.out_dir = self.out if isdir(self.out) else dirname(self.out)
        self.version_label = "hg{}.imgt{}".format(self.hg, self.imgt)

        self.imgt_dir = _imgt_dir if Exists(_imgt_dir, isdir) \
            else RaiseError(HATK_InputPreparation_Error, "The given directory for the IMGT database('{}') isn't a directory.".format(_imgt_dir))
        self.data_dir = _p_data if Exists(_p_data, isdir) \
            else RaiseError(HATK_InputPreparation_Error, "The given data directory('{}') isn't a directory.".format(_p_data))

        # Requested HLAs
        self.HLA_req = _HLA_req
        self.HLA_avail = getAvailableHLAs(self.imgt_dir, self.data_dir, self.hg)
        self.HLA_target, self.HLA_excluded1, self.HLA_excluded2 = getTargetHLAs(self.HLA_req, *self.HLA_avail)

        # Output Field format
        self.F_one = _F_one
        self.F_two = _F_two
        self.F_three = _F_three
        self.F_four = _F_four
        self.F_Ggroup = _F_Ggroup
        self.F_Pgroup = _F_Pgroup
        self.which_format = None

        # Output results
        self.__dict_AA__ = None # Prefix
        self.__dict_AA_seq__ = None
        self.__dict_AA_map__ = None
        self.__dict_SNPS__ = None # Prefix
        self.__dict_SNPS_seq__ = None
        self.__dict_SNPS_map__ = None
        self.__HAT__ = None
        self.__d_MapTable__ = None

        # Optional arguments
        self.no_indel = _no_Indel
        self.no_prime = _no_prime
        self.multiprocess = _MultiP
        self.save_intermediates = _save_intermediates

        # flags
        self.f_hasPreviousResult = self.findExistingResult()


        ### Main Actions ###
        if self.f_hasPreviousResult:
            print(std_WARNING_MAIN_PROCESS_NAME.lstrip('\n')
                  + "Using the previous result in the output directory('{}')\n".format(self.out_dir))

            print(self.__repr__())

        else:

            if self.F_one:
                self.which_format = 1
            elif self.F_two:
                self.which_format = 2
            elif self.F_three:
                self.which_format = 3
            elif self.F_four:
                self.which_format = 4
            elif self.F_Ggroup:
                self.which_format = 5
            elif self.F_Pgroup:
                self.which_format = 6
            else:
                self.which_format = 2
                print(std_WARNING_MAIN_PROCESS_NAME.lstrip('\n') +
                      "Which output field format to use wans't given. IMGT2Seq 'overrides' it as '2-field'. "
                      "Please check '--1field', '--2field', ..., '--Ggroup' and '--Pgroup' arguments again.")

            # print(self.__repr__()) # debug
            self.performIMGT2Seq()
            self.findExistingResult()
            print(self.__repr__())


    def findExistingResult(self):
        """
        Check whether some of the output exists in the output directory folder.
        When the last output of IMGT2Seq exists in the `outdir` directory, Use that instead of generating a new one redundantly.

        """

        # (1) dict_AA
        t_dict_AA = join(self.out_dir, "HLA_DICTIONARY_AA.{}".format(self.version_label))

        self.__dict_AA_seq__ = t_dict_AA+".txt" if exists(t_dict_AA+".txt") else None
        self.__dict_AA_map__ = t_dict_AA+".map" if exists(t_dict_AA+".map") else None
        self.__dict_AA__ = t_dict_AA if self.__dict_AA_seq__ and self.__dict_AA_map__ else None

        # (2) dict_SNPS
        t_dict_SNPS = join(self.out_dir, "HLA_DICTIONARY_SNPS.{}".format(self.version_label))

        self.__dict_SNPS_seq__ = t_dict_SNPS+".txt" if exists(t_dict_SNPS+".txt") else None
        self.__dict_SNPS_map__ = t_dict_SNPS+".map" if exists(t_dict_SNPS+".map") else None
        self.__dict_SNPS__ = t_dict_SNPS if self.__dict_SNPS_seq__ and self.__dict_SNPS_map__ else None

        # (3) hat
        t_hat = join(self.out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(self.imgt))
        self.__HAT__ = t_hat if exists(t_hat) else None

        # (4) maptables
        self.__d_MapTable__ = {hla: None for hla in self.HLA_target}

        for hla in self.HLA_target:
            t_maptable = join(self.out_dir, "HLA_MAPTABLE_{}.{}.txt".format(hla, self.version_label))
            self.__d_MapTable__[hla] = t_maptable if exists(t_maptable) else None


        return self.__bool__()


    def performIMGT2Seq(self):
        # self.__dict_AA__, self.__dict_SNPS__, self.__HAT__, self.__d_MapTable__ = \
        #     IMGT2Seq(self.imgt, self.hg, self.out, _imgt_dir=self.imgt_dir, _no_Indel=self.no_indel,
        #              _MultiP=self.multiprocess, _save_intermediates=self.save_intermediates,
        #              _p_data="IMGT2Seq/data", __Nfield_OUTPUT_FORMAT=self.which_format)

        temp = \
            IMGT2Seq(self.imgt, self.hg, self.out, _imgt_dir=self.imgt_dir, _no_Indel=self.no_indel,
                     _MultiP=self.multiprocess, _save_intermediates=self.save_intermediates, _p_data="IMGT2Seq/data",
                     __Nfield_OUTPUT_FORMAT=self.which_format, _HLA_target=self.HLA_target)


    def __repr__(self):

        str_header = \
                "< IMGT2Sequence Summary (Using existing results) >\n" if self.f_hasPreviousResult else \
                "< IMGT2Sequence Summary (Newly generated.) >\n"

        str_imgt_version = \
            "- IMGT version: v{}\n".format(self.imgt)

        str_hg = \
            "- Human Genome version: hg{}\n".format(self.hg)

        str_out_dir = \
            "- Output directory: {}\n".format(self.out_dir)

        str_HLA = \
            "- HLA:\n" \
            "   (requested): {}\n" \
            "   (excluded1): {}   (No 3 '*_prot.txt', '*_nuc.txt', and '*_gen.txt' files in '{}')\n" \
            "   (excluded2): {}   (No exon1 start BP information(in '{}'))\n" \
            "   (target): {}\n" \
            .format(self.HLA_req,
                    self.HLA_excluded1, self.imgt_dir,
                    self.HLA_excluded2, join(self.data_dir, "HLA_EXON1_START_CODON_POSITIONS_hg{}.txt".format(self.hg)),
                    list(self.HLA_target))

        str_Nfield = \
            "" if self.f_hasPreviousResult else \
            "- Requested Output format field: {}\n".format(
                "1-field" if self.which_format == 1 else \
                "2-field" if self.which_format == 2 else \
                "3-field" if self.which_format == 3 else \
                "4-field" if self.which_format == 4 else \
                "G-group" if self.which_format == 5 else \
                "P-group" if self.which_format == 6 else None
            )

        str_dict_AA = \
            "- Dictionary(Amino acids)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.__dict_AA_seq__, self.__dict_AA_map__)

        str_dict_SNPS = \
            "- Dictionary(Intragenic SNPs)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.__dict_SNPS_seq__, self.__dict_SNPS_map__)

        str_HLA_Allele_Table = \
            "- HLA Allele Table:\n" \
            "   {}\n".format(self.__HAT__)

        str_maptable = \
            "- HLA Maptables for Heatmap plot: \n" \
            + ''.join([
                "   {hla}   : {maptable}\n".format(hla=hla, maptable=self.__d_MapTable__[hla]) for hla in self.HLA_target
            ])

        str_summary = ''.join([
            str_header,
            str_imgt_version,
            str_hg,
            str_out_dir,
            str_HLA,
            str_Nfield,
            str_dict_AA,
            str_dict_SNPS,
            str_HLA_Allele_Table,
            str_maptable
        ]).rstrip('\n')

        return str_summary


    def __bool__(self):
        return Exists(self.__dict_AA_seq__) and Exists(self.__dict_AA_map__) and \
               Exists(self.__dict_SNPS_seq__) and Exists(self.__dict_SNPS_map__) and \
               Exists(self.__HAT__) and \
               (len(self.__d_MapTable__) == len(self.HLA_target) and np.all([Exists(x) for x in self.__d_MapTable__.values()]))


    def getResult(self):
        return self.__dict_AA__, self.__dict_SNPS__, self.__HAT__, self.__d_MapTable__





if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='IMGT2Seq',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        IMGT2Sequence



    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
    parser.add_argument("--hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)
    parser.add_argument("--out", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser.add_argument("--imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                        required=True)

    parser.add_argument("--imgt-dir", help="\nIn case User just want to specify the directory of IMGT data folder.\n\n",
                        required=True)


    parser.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n",
                        default=['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')

    parser.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    parser.add_argument("--multiprocess", help="\nSetting off parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')


    # Output format selection
    format_selection = parser.add_mutually_exclusive_group()
    format_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n", action="store_true", dest="F_one")
    format_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n", action="store_true", dest="F_two")
    format_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n", action="store_true", dest="F_three")
    format_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields.\n\n", action="store_true", dest="F_four")
    format_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n", action="store_true", dest="F_Ggroup")
    format_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n", action="store_true", dest="F_Pgroup")



    ##### < for Test > #####

    # hg18 / imgt3320
    # args = parser.parse_args(["--hg", "18",
    #                           "--out", out,
    #                           "--imgt", "3320",
    #                           "--imgt-dir", imgt_dir,
    #                           '--save-intermediates'])

    # hg18 / imgt3320
    # args = parser.parse_args(["--hg", "18",
    #                           "--out", out,
    #                           "--imgt", "3320"])



    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    ##### < Main function Execution. > #####


    HATK_IMGT2Seq(args.imgt, args.hg, args.out, args.F_one, args.F_two, args.F_three, args.F_four, args.F_Ggroup,
                  args.F_Pgroup, _imgt_dir=args.imgt_dir, _no_Indel=args.no_indel, _MultiP=args.multiprocess,
                  _save_intermediates=args.save_intermediates, _HLA_req=args.HLA)