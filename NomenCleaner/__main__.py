#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join
import argparse, textwrap
import numpy as np
import pandas as pd

from NomenCleaner.NomenCleaner import NomenCleaner
from NomenCleaner.HPED import HPED
from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.util import Exists

std_MAIN_PROCESS_NAME = "\n[NomenCleaner]: "
std_ERROR_MAIN_PROCESS_NAME = "\n[NomenCleaner::ERROR]: "
std_WARNING_MAIN_PROCESS_NAME = "\n[NomenCleaner::WARNING]: "


class HATK_NomenCleaner(object):

    def __init__(self, _hped, _hat, _out, _F_one, _F_two, _F_three, _F_four, _F_Ggroup, _F_Pgroup,
                 _f_NoGenePrefix=False, _f_leave_NotFound=False,
                 _HLA_req:tuple=None):

        """

        """
        ### Main Variables ###
        self.HPED = HPED(_hped)
        self.hped_temp = None # Filtered hped.

        self.hat = _hat if Exists(_hat) \
            else RaiseError(HATK_InputPreparation_Error, std_ERROR_MAIN_PROCESS_NAME + "Given HAT file('{}') can't be found.".format(_hat))
        self.imgt = re.search('imgt(\d+)', basename(self.hat)).group(1) if bool(re.search('imgt(\d+)', basename(self.hat))) else None


        self.out = _out
        self.out_prefix = self.out.rstrip('.chped') if self.out.endswith('.chped') else self.out
        self.out_prefix = self.out_prefix + (".imgt{}".format(self.imgt) if self.imgt else "")

        self.HLA_req = _HLA_req if _HLA_req else self.HPED.HLA_avail # Use HLAs in HPED file if user didn't request any HLA genes.

        self.HLA_avail_hped = self.HPED.HLA_avail
        self.HLA_avail_hat = np.unique(np.genfromtxt(self.hat, delimiter='\t', dtype=str, usecols=[0], skip_header=1))
        self.HLA_excluded_hped = np.setdiff1d(self.HLA_req, self.HLA_avail_hped)
        self.HLA_excluded_hat = np.setdiff1d(self.HLA_req, self.HLA_avail_hat)

        self.HLA_target = np.setdiff1d(self.HLA_req, np.union1d(self.HLA_excluded_hped, self.HLA_excluded_hat)) if self.HPED.f_hasHeader else \
                            ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")


        # Output Field format
        self.F_one = _F_one
        self.F_two = _F_two
        self.F_three = _F_three
        self.F_four = _F_four
        self.F_Ggroup = _F_Ggroup
        self.F_Pgroup = _F_Pgroup
        self.which_format = 1 if self.F_one else \
                            2 if self.F_two else \
                            3 if self.F_three else \
                            4 if self.F_four else \
                            5 if self.F_Ggroup else \
                            6 if self.F_Pgroup else \
                            RaiseError(HATK_InputPreparation_Error, "Wrong Output Field Format.")

        # Optional Flags
        self.f_NoGenePrefix = _f_NoGenePrefix
        self.f_leave_NotFound = _f_leave_NotFound

        # Output results
        self.chped = None
        self.chped_log = None

        # to Remove
        self.l_toRemove = []


        ### Main Actions ###
        # print(self.__repr__())
        self.genTempHPED()

        self.chped, self.chped_log = \
            NomenCleaner(self.hped_temp, self.hat, self.out_prefix, self.which_format,
                         _f_NoGenePrefix=self.f_NoGenePrefix, _f_leave_NotFound=self.f_leave_NotFound,
                         _HLA_target=self.HLA_target)

        print(self.__repr__())


    def genTempHPED(self):

        if self.HPED.f_hasHeader:

            # a temp HPED file with `self.HLA_target`.
            # this temp HPED file will be used by NomenCleaner().

            l_temp = np.array([(hla+'_1', hla+'_2') for hla in self.HLA_target]).flatten()
            df_temp = pd.concat([
                self.HPED.df_hped_Left,
                self.HPED.df_hped_Right.loc[:, l_temp]
            ], axis=1)
            # print("Temp HPED to be used by NomenCleaner():\n{}\n".format(df_temp))

            self.hped_temp = self.out_prefix + '.temp.hped'
            df_temp.to_csv(self.hped_temp, sep='\t', header=False, index=False)

            # Things to remove
            self.l_toRemove.append(self.hped_temp)

        else:
            self.hped_temp = self.HPED.hped


    def __repr__(self):
        str_HPED = \
            "- HPED file: {}\n".format(self.HPED.hped)

        str_hasHeader = \
            "- HPED has Header?: {}\n".format(self.HPED.f_hasHeader)

        str_HAT = \
            "- HAT file: {}\n".format(self.hat)

        str_imgt_version = \
            "- IMGT version: v{}\n".format(self.imgt)

        str_Nfield = \
            "- Requested Output Field format: {}\n".format(
                "1-field" if self.which_format == 1 else \
                "2-field" if self.which_format == 2 else \
                "3-field" if self.which_format == 3 else \
                "4-field" if self.which_format == 4 else \
                "G-group" if self.which_format == 5 else \
                "P-group" if self.which_format == 6 else None
            )

        str_out_dir = \
            "- Output prefix: {}(.chped)\n".format(self.out_prefix)

        str_HLA = \
            "- HLA:\n" \
            "   (requested): {}\n" \
            "   (available - hped): {}\n" \
            "   (available - hat): {}\n" \
            "   (excluded - Not in hped): {}\n" \
            "   (excluded - Not in hat): {}\n" \
            "   (target): {}\n" \
            .format(list(self.HLA_req),
                    list(self.HLA_avail_hped), list(self.HLA_avail_hat),
                    list(self.HLA_excluded_hped), list(self.HLA_excluded_hat),
                    list(self.HLA_target))

        str_output = \
            "- CHPED output: {}\n" \
            "- CHPED log: {}\n" \
            .format(self.chped, self.chped_log)


        str_summary = ''.join([
            str_HPED,
            str_hasHeader,
            str_HAT,
            str_imgt_version,
            str_Nfield,
            str_out_dir,
            str_HLA,
            str_output
        ]).rstrip('\n')

        return str_summary


    def __bool__(self): return Exists(self.chped)


    def __del__(self):
        for item in self.l_toRemove:
            if Exists(item): os.remove(item)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='NomenCleaner',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        < NomenCleaner.py >

        - Transforms *.hped file to *.chped file.
        - *.hat file must be given as input file.



    #########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--hped", help="\nHLA Type Data which hasn't gone through NomenCleaner yet.\n\n",
                          dest="hped", required=True)

    parser.add_argument("--hat", help="\nHLA Allele Table file(*.hat).\n\n", required=True)
    # parser.add_argument("--imgt", help="\nSpecifying the IMGT-HLA version.\n\n", required=True)

    parser.add_argument("--out", help="\nOutput file prefix.\n\n", required=True)
    parser.add_argument("--leave-NotFound",
                        help="\nLeaving HLA alleles which can't be found in given *.hat file(Novel or Erroneous allele) intact.\n\n",
                        action='store_true')

    parser.add_argument("--HLA", help="\nHLA genes to process.\n\n", nargs='+')

    format_selection = parser.add_mutually_exclusive_group(required=True)
    format_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n", action="store_true", dest="F_one")
    format_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n", action="store_true", dest="F_two")
    format_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n", action="store_true", dest="F_three")
    format_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields.\n\n", action="store_true", dest="F_four")
    format_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n", action="store_true", dest="F_Ggroup")
    format_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n", action="store_true", dest="F_Pgroup")

    # Flag to remove HLA gene caption.
    parser.add_argument("--NoGenePrefix", help="\nMake converted HLA alleles NOT have HLA gene prefix(ex. \"A*\").\n\n",
                        action='store_true')

    ##### <for Test> #####


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    HATK_NomenCleaner(args.hped, args.hat, args.out,
                      args.F_one, args.F_two, args.F_three, args.F_four, args.F_Ggroup, args.F_Pgroup,
                      _f_NoGenePrefix=args.NoGenePrefix, _f_leave_NotFound=args.leave_NotFound,
                      _HLA_req=args.HLA)