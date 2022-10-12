# -*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join, isdir
from shutil import which
import argparse, textwrap

from bMarkerGenerator.bMarkerGenerator import bMarkerGenerator
from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.util import Exists, FieldFormat2Label, findExec
from src.PLINK_GT import GT
from NomenCleaner.src.CHPED import CHPED
from IMGT2Seq.src.IMGT2Seq_Output import HLA_DICTIONARY

std_MAIN = "\n[bMarkerGenerator]: "
std_ERROR = "\n[bMarkerGenerator::ERROR]: "
std_WARNING = "\n[bMarkerGenerator::WARNING]: "


class HATK_bMarkerGenertor(object):

    # External software
    plink = findExec("plink", std_ERROR+"'plink' command can't be found. Please install PLINK.")

    def __init__(self, _chped, _out_prefix, _hg, _dictionary_AA, _dictionary_SNPS, _bfile=None,
                 _f_save_intermediates=False):

        """

        """

        ### Main Variables ###
        self.out_prefix = _out_prefix
        self.out_dir = dirname(self.out_prefix)
        self.hg = _hg

        # PLINK Genotype
        self.bfile = GT(_bfile) if bool(_bfile) else None

        # CHPED
        self.CHPED = CHPED(_chped)

        # IMGT2Seq HLA dictionary
        self.HLA_DICTIONARY = \
            HLA_DICTIONARY(_dictionary_AA + '.txt', _dictionary_AA + '.map',
                           _dictionary_SNPS + '.txt', _dictionary_SNPS + '.map',
                           _HLA_req=self.CHPED.HLA_avail)
        """
        Intersection between HLA gene sets from CHPED and HLA_DICTIONARY is needed for bMG.
        
        This intersection is acquired in the above 'HLA_DICTIONARY' class with the '_HLA_req' argument.
        
        In the end, `HLA_DICTIONARY.HLA_target` variable will contain the intersected HLA genes.
        """

        # Output(bmarker)
        self.bmarker = None

        # Flag
        self.f_save_intermediates = _f_save_intermediates


        ### Main Actions ###
        # print(self.__repr__())

        # Consenting HLA gene set.
        if self.CHPED.HLA_avail != self.HLA_DICTIONARY.HLA_avail:
            print(std_WARNING +
                  "Available HLA gene sets of given CHPED and HLA_DICTIONARY are different.\n"
                  "   (CHPED): {}\n"
                  "   (HLA_DICT): {}".format(self.CHPED.HLA_avail, self.HLA_DICTIONARY.HLA_avail))
            print(std_WARNING.lstrip("\n") +
                  "Next intersected HLA gene set will be used\n"
                  "   (Intersected): {}".format(self.HLA_DICTIONARY.HLA_target))
            """
            `HLA_DICTIONARY.HLA_target` is the result of intersection between CHPED and HLA_DICTIONARY.
            """

            self.CHPED = self.CHPED.subsetCHPED(join(self.out_dir, basename(_chped)+'.subset'), self.HLA_DICTIONARY.HLA_target)
            self.HLA_DICTIONARY = self.HLA_DICTIONARY.subsetHLA_DICTIONARY(self.out_dir, self.HLA_DICTIONARY.HLA_target)

            # print(std_MAIN + "New CHPED subsetted to the target genes.")
            # print(self.CHPED)
            # print(std_MAIN + "New HLA Dictionary files subsetted to the target genes.")
            # print(self.HLA_DICTIONARY)



        if not (self.CHPED.field_format == self.HLA_DICTIONARY.field_format):
            raise HATK_InputPreparation_Error(
                std_ERROR + "Each field format of CHPED and HLA dictionary files doesn't match.\n"
                            "-(CHPED): {}\n"
                            "-(HLA Dictionary): {}"\
                .format(FieldFormat2Label(self.CHPED.field_format), FieldFormat2Label(self.HLA_DICTIONARY.field_format))
            )

        self.bmarker = \
            bMarkerGenerator(self.CHPED, self.out_prefix, self.hg,
                             self.HLA_DICTIONARY.HLA_dict_AA_seq, self.HLA_DICTIONARY.HLA_dict_AA_map,
                             self.HLA_DICTIONARY.HLA_dict_SNPS_seq, self.HLA_DICTIONARY.HLA_dict_SNPS_map,
                             _variants=self.bfile.file_prefix, _f_save_intermediates=self.f_save_intermediates,
                             _plink=self.plink)

        # print(self.__repr__())



    def __repr__(self):

        str_PLINK = \
            "=====< INPUT(1) - PLINK Genotype >=====\n{}\n".format(self.bfile) if bool(self.bfile) else ""

        str_CHPED = \
            "\n=====< INPUT(2) - CHPED(HLA type) >=====\n{}\n".format(self.CHPED)

        str_HLA_DICTIONARY = \
            "\n=====< INPUT(3) - HLA Dictionary(IMGT2Seq) >=====\n{}\n".format(self.HLA_DICTIONARY)

        str_bMG = \
            "\n< bMarkerGenerator >\n"
        str_hg = \
            "- Human Genome version: hg{}\n".format(self.hg)
        str_out = \
            "- Output prefix: {}\n".format(self.out_prefix)
        str_bmarker = \
            "- Generated bMarker:\n{}\n".format(self.bmarker)


        str_summary = ''.join([
            str_PLINK, str_CHPED, str_HLA_DICTIONARY,
            str_bMG, str_hg, str_out, str_bmarker
        ]).rstrip('\n')

        return str_summary



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='bMarkerGenerator',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        bMarkerGenerator.py

        Generating markers based on HLA sequence information dictionary(generated by "IMGT2Seq").


        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--bfile", help="\nInput variants data file(.bed/.bim/.fam)\n\n")
    parser.add_argument("--chped", help="\nHLA Type Data(.chped)\n\n", required=True)
    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg",
                        default="19", required=True)
    parser.add_argument("--out", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", required=True)
    parser.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')

    ##### <for Test> #####

    # 2019. 01. 10
    # args = parser.parse_args(["-chped", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU_HLA.imgt370.4field.chped",
    #                           "-variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
    #                           "-hg", "18",
    #                           "-dict-AA", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_AA.hg18.imgt370",
    #                           "-dict-SNPS", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_SNPS.hg18.imgt370"
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    # # Implementing Main Function.
    # bMarkerGenerator(_CHPED=args.chped, _OUT=args.out, _hg=args.hg, _variants=args.variants,
    #                  _dictionary_AA=args.dict_AA, _dictionary_SNPS=args.dict_SNPS,
    #                  __save_intermediates=args.save_intermediates)

    HATK_bMarkerGenertor(args.chped, args.out, args.hg, args.dict_AA, args.dict_SNPS, _bfile=args.bfile,
                         _f_save_intermediates=args.save_intermediates)