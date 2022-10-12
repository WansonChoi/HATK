#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import argparse, textwrap
import numpy as np

from IMGT2Seq.IMGT2Seq import IMGT2Seq
from IMGT2Seq.src.AvailableHLAs import getAvailableHLAs, getTargetHLAs
from src.HATK_Error import RaiseError, HATK_InputPreparation_Error
from src.util import Exists, FieldFormat2Label, which_format
from IMGT2Seq.src.IMGT2Seq_Output import IMGT2Seq_Output

std_MAIN = "\n[IMGT2Seq]: "
std_ERROR = "\n[IMGT2Seq::ERROR]: "
std_WARNING = "\n[IMGT2Seq::WARNING]: "


class HATK_IMGT2Seq(object):

    def __init__(self, _imgt, _hg, _out, _imgt_dir,
                 _F_one:bool, _F_two:bool, _F_three:bool, _F_four:bool, _F_Ggroup:bool, _F_Pgroup:bool,
                 _no_Indel=False, _MultiP=False, _f_save_intermediates=False, _no_prime=True, _p_data='./IMGT2Seq/data',
                 _HLA_req=None):

        """
        Either (1) generate a new Dictionary files or (2) find existing dictionaries.

        (2022.02.15.) Compare the requested HLAs vs. Available HLAs.
        """

        ### Main Variables ###
        self.imgt = _imgt
        self.hg = _hg
        self.out = _out

        self.out_dir = self.out if isdir(self.out) else dirname(self.out)
        # self.version_label = "hg{}.imgt{}".format(self.hg, self.imgt)

        self.imgt_dir = _imgt_dir if Exists(_imgt_dir, isdir) \
            else RaiseError(HATK_InputPreparation_Error, "Given directory for the IMGT database('{}') isn't a directory.".format(_imgt_dir))
        self.data_dir = _p_data if Exists(_p_data, isdir) \
            else RaiseError(HATK_InputPreparation_Error, "Given data directory('{}') isn't a directory.".format(_p_data))

        # Requested HLAs
        self.HLA_req = _HLA_req
        self.HLA_avail = getAvailableHLAs(self.imgt_dir, self.data_dir, self.hg)
        if _HLA_req:
            self.HLA_target, self.HLA_excluded1, self.HLA_excluded2 = \
                getTargetHLAs(self.HLA_req, *self.HLA_avail)
        else:
            # If no HLA_req, then take all available.
            self.HLA_target = list(np.intersect1d(*self.HLA_avail))
            self.HLA_excluded1 = list(np.setdiff1d(self.HLA_avail[0], self.HLA_target))
            self.HLA_excluded2 = list(np.setdiff1d(self.HLA_avail[1], self.HLA_target))


        # Output Field format
        self.which_format = which_format(_F_one, _F_two, _F_three, _F_four, _F_Ggroup, _F_Pgroup)

        # Optional arguments
        self.no_indel = _no_Indel
        self.no_prime = _no_prime
        self.multiprocess = _MultiP
        self.f_save_intermediates = _f_save_intermediates


        ### Main Actions ###
        self.IMGT2Seq_Output = \
            IMGT2Seq(self.imgt, self.hg, self.out_dir, _imgt_dir=self.imgt_dir, _no_Indel=self.no_indel,
                     _MultiP=self.multiprocess, _f_save_intermediates=self.f_save_intermediates, _p_data="IMGT2Seq/data",
                     __Nfield_OUTPUT_FORMAT=self.which_format, _HLA_target=self.HLA_target)
        # print(self.__repr__()) # debug


    def __repr__(self):

        # str_header = \
        #         "< IMGT2Sequence Summary (Using existing results) >\n" if self.f_hasPreviousResult else \
        #         "< IMGT2Sequence Summary (Newly generated.) >\n"

        str_imgt_version = \
            "- IMGT version: v{}\n".format(self.imgt)

        str_hg = \
            "- Human Genome version: hg{}\n".format(self.hg)

        str_out_dir = \
            "- Output directory: {}\n".format(self.out_dir)

        # str_Nfield = \
        #     "" if self.f_hasPreviousResult else \
        #     "- Requested Output Field format: {}\n".format(FieldFormat2Label(self.which_format))
        str_Nfield = \
            "- Requested Output Field format: {}\n".format(FieldFormat2Label(self.which_format))

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

        str_dict_AA = \
            "- Dictionary(Amino acids)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_AA_seq,
                    self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_AA_map)
            # .format(self.__dict_AA_seq__, self.__dict_AA_map__)

        str_dict_SNPS = \
            "- Dictionary(Intragenic SNPs)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_SNPS_seq,
                    self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_SNPS_map)
            # .format(self.__dict_SNPS_seq__, self.__dict_SNPS_map__)

        str_HLA_Allele_Table = \
            "- HLA Allele Table:\n" \
            "   {}\n".format(self.IMGT2Seq_Output.HAT)
            # "   {}\n".format(self.__HAT__)

        str_maptable = \
            "- HLA Maptables for Heatmap plot: \n" \
            + ''.join([
                "   {hla}   : {maptable}\n" \
                    .format(hla=hla,
                            maptable=self.IMGT2Seq_Output.HLA_MAPTABLE[hla] if hla in self.IMGT2Seq_Output.HLA_MAPTABLE else None) \
                for hla in self.HLA_target
                    # .format(hla=hla, maptable=self.__d_MapTable__[hla]) for hla in self.HLA_target
            ])

        str_summary = ''.join([
            # str_header,
            str_imgt_version,
            str_hg,
            str_out_dir,
            str_Nfield,
            str_HLA,
            str_dict_AA,
            str_dict_SNPS,
            str_HLA_Allele_Table,
            str_maptable
        ]).rstrip('\n')

        return str_summary


    def __bool__(self):
        return bool(self.IMGT2Seq_Output)


# class HATK_IMGT2Seq(object):
#     def __init__(self, _imgt, _hg, _out, _imgt_dir,
#                  _F_one:bool, _F_two:bool, _F_three:bool, _F_four:bool, _F_Ggroup:bool, _F_Pgroup:bool,
#                  _no_Indel=False, _MultiP=False, _f_save_intermediates=False, _no_prime=True, _p_data='./IMGT2Seq/data',
#                  _HLA_req=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):
#
#         """
#         Either (1) generate a new Dictionary files or (2) find existing dictionaries.
#
#         (2022.02.15.) Compare the requested HLAs vs. Available HLAs.
#         """
#
#         ### Main Variables ###
#         self.imgt = _imgt
#         self.hg = _hg
#         self.out = _out
#
#         self.out_dir = self.out if isdir(self.out) else dirname(self.out)
#         self.version_label = "hg{}.imgt{}".format(self.hg, self.imgt)
#
#         self.imgt_dir = _imgt_dir if Exists(_imgt_dir, isdir) \
#             else RaiseError(HATK_InputPreparation_Error, "Given directory for the IMGT database('{}') isn't a directory.".format(_imgt_dir))
#         self.data_dir = _p_data if Exists(_p_data, isdir) \
#             else RaiseError(HATK_InputPreparation_Error, "Given data directory('{}') isn't a directory.".format(_p_data))
#
#         # Requested HLAs
#         self.HLA_req = _HLA_req
#         self.HLA_avail = getAvailableHLAs(self.imgt_dir, self.data_dir, self.hg)
#         self.HLA_target, self.HLA_excluded1, self.HLA_excluded2 = getTargetHLAs(self.HLA_req, *self.HLA_avail)
#
#         # Output Field format
#         self.which_format = which_format(_F_one, _F_two, _F_three, _F_four, _F_Ggroup, _F_Pgroup)
#
#         # Output result
#         self.IMGT2Seq_Output = IMGT2Seq_Output(self.out_dir, self.HLA_target)
#         self.f_hasPreviousResult = bool(self.IMGT2Seq_Output)
#
#         # Optional arguments
#         self.no_indel = _no_Indel
#         self.no_prime = _no_prime
#         self.multiprocess = _MultiP
#         self.f_save_intermediates = _f_save_intermediates
#
#
#         ### Main Actions ###
#         print(self.__repr__()) # debug
#         print(self.IMGT2Seq_Output)
#         print(bool(self.IMGT2Seq_Output))
#
#         if self.f_hasPreviousResult:
#             print(std_MAIN.lstrip('\n')
#                   + "Using the previous result in the output directory('{}')\n".format(self.out_dir))
#
#         else:
#             self.IMGT2Seq_Output = \
#                 IMGT2Seq(self.imgt, self.hg, self.out, _imgt_dir=self.imgt_dir, _no_Indel=self.no_indel,
#                          _MultiP=self.multiprocess, _f_save_intermediates=self.f_save_intermediates, _p_data="IMGT2Seq/data",
#                          __Nfield_OUTPUT_FORMAT=self.which_format, _HLA_target=self.HLA_target)
#             # pass
#
#
#     def __repr__(self):
#
#         str_header = \
#                 "< IMGT2Sequence Summary (Using existing results) >\n" if self.f_hasPreviousResult else \
#                 "< IMGT2Sequence Summary (Newly generated.) >\n"
#
#         str_imgt_version = \
#             "- IMGT version: v{}\n".format(self.imgt)
#
#         str_hg = \
#             "- Human Genome version: hg{}\n".format(self.hg)
#
#         str_out_dir = \
#             "- Output directory: {}\n".format(self.out_dir)
#
#         str_Nfield = \
#             "" if self.f_hasPreviousResult else \
#             "- Requested Output Field format: {}\n".format(FieldFormat2Label(self.which_format))
#
#         str_HLA = \
#             "- HLA:\n" \
#             "   (requested): {}\n" \
#             "   (excluded1): {}   (No 3 '*_prot.txt', '*_nuc.txt', and '*_gen.txt' files in '{}')\n" \
#             "   (excluded2): {}   (No exon1 start BP information(in '{}'))\n" \
#             "   (target): {}\n" \
#             .format(self.HLA_req,
#                     self.HLA_excluded1, self.imgt_dir,
#                     self.HLA_excluded2, join(self.data_dir, "HLA_EXON1_START_CODON_POSITIONS_hg{}.txt".format(self.hg)),
#                     list(self.HLA_target))
#
#         str_dict_AA = \
#             "- Dictionary(Amino acids)\n" \
#             "   (seq): {}\n" \
#             "   (map): {}\n" \
#             .format(self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_AA_seq,
#                     self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_AA_map)
#             # .format(self.__dict_AA_seq__, self.__dict_AA_map__)
#
#         str_dict_SNPS = \
#             "- Dictionary(Intragenic SNPs)\n" \
#             "   (seq): {}\n" \
#             "   (map): {}\n" \
#             .format(self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_SNPS_seq,
#                     self.IMGT2Seq_Output.HLA_DICTIONARY.HLA_dict_SNPS_map)
#             # .format(self.__dict_SNPS_seq__, self.__dict_SNPS_map__)
#
#         str_HLA_Allele_Table = \
#             "- HLA Allele Table:\n" \
#             "   {}\n".format(self.IMGT2Seq_Output.HAT)
#             # "   {}\n".format(self.__HAT__)
#
#         str_maptable = \
#             "- HLA Maptables for Heatmap plot: \n" \
#             + ''.join([
#                 "   {hla}   : {maptable}\n" \
#                     .format(hla=hla,
#                             maptable=self.IMGT2Seq_Output.HLA_MAPTABLE[hla] if hla in self.IMGT2Seq_Output.HLA_MAPTABLE else None) \
#                 for hla in self.HLA_target
#                     # .format(hla=hla, maptable=self.__d_MapTable__[hla]) for hla in self.HLA_target
#             ])
#
#         str_summary = ''.join([
#             str_header,
#             str_imgt_version,
#             str_hg,
#             str_out_dir,
#             str_Nfield,
#             str_HLA,
#             str_dict_AA,
#             str_dict_SNPS,
#             str_HLA_Allele_Table,
#             str_maptable
#         ]).rstrip('\n')
#
#         return str_summary
#
#
#     def __bool__(self):
#         return bool(self.IMGT2Seq_Output)







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
    # parser.add_argument("--multiprocess", help="\nSetting off parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    parser.add_argument("--multiprocess", help="\nSetting off parallel multiprocessing.\n\n", type=int, nargs='?', default=1, const=8)
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

    print(sys.argv[0:])
    args = parser.parse_args()
    print(args)


    ##### < Main function Execution. > #####

    HATK_IMGT2Seq(args.imgt, args.hg, args.out, args.imgt_dir, args.F_one, args.F_two, args.F_three, args.F_four, args.F_Ggroup,
                  args.F_Pgroup, _no_Indel=args.no_indel, _MultiP=args.multiprocess,
                  _f_save_intermediates=args.save_intermediates, _HLA_req=args.HLA)