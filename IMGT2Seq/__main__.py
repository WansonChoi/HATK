import os, sys, re
import argparse, textwrap

from IMGT2Seq.IMGT2Seq import IMGT2Seq
from IMGT2Seq.src.IMGT2SeqError import IMGT2SeqError, ArgumentCheck



########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % ("IMGT2Seq")
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % ("IMGT2Seq")
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % ("IMGT2Seq")

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



class HATK_IMGT2Seq(object):

    # @ArgumentCheck
    def __init__(self, _imgt, _hg, _out, _HLA, *args, **kwargs):

        """

        """


        ########## < Assigning arguments > ##########

        # Main positional arguments
        self.imgt = _imgt
        self.hg = _hg
        self.out = _out
        self.HLA = _HLA

        # Output Field format (*args)
        oneF = args[0]
        twoF = args[1]
        threeF = args[2]
        fourF = args[3]
        Ggroup = args[4]
        Pgroup = args[5]

        # Optional arguments (**kwargs)
        self.no_Ins = kwargs["_no_Ins"]
        self.multiprocess = kwargs["_multiprocess"]
        self.save_intermediates = kwargs["_save_intermediates"]
        self.imgt_dir = kwargs["_imgt_dir"]


        # Main Outputs
        self.__HAT__ = None
        self.__dict_AA__ = None
        self.__dict_SNPS__ = None
        self.__d_MapTable__ = None

        # Summary string
        self.summary_string = ""





        ########## < Assigning arguments > ##########

        ### Check output nomenclature
        if oneF:
            print(std_MAIN_PROCESS_NAME + "Selected Output Nomenclature : 1-field")
            Nfield_OUTPUT_FORMAT = 1
        elif twoF:
            print(std_MAIN_PROCESS_NAME + "Selected Output Nomenclature : 2-field")
            Nfield_OUTPUT_FORMAT = 2
        elif threeF:
            print(std_MAIN_PROCESS_NAME + "Selected Output Nomenclature : 3-field")
            Nfield_OUTPUT_FORMAT = 3
        elif fourF:
            print(std_MAIN_PROCESS_NAME + "Selected Output Nomenclature : 4-field")
            Nfield_OUTPUT_FORMAT = 4
        elif Ggroup:
            print(std_MAIN_PROCESS_NAME + "Selected Output Nomenclature : G-group")
            Nfield_OUTPUT_FORMAT = 5
        elif Pgroup:
            print(std_MAIN_PROCESS_NAME + "Selected Output Nomenclature : P-group")
            Nfield_OUTPUT_FORMAT = 6
        else:
            # Wrong
            raise IMGT2SeqError(std_ERROR_MAIN_PROCESS_NAME + "Wrong Output Nomenclature.")



        ### No more using existing result.
        self.__dict_AA__, self.__dict_SNPS__, self.__HAT__, self.__d_MapTable__ = \
            IMGT2Seq(self.imgt, self.hg, self.out, _imgt_dir=self.imgt_dir, _no_Indel=self.no_Ins,
                     _MultiP=self.multiprocess, _save_intermediates=self.save_intermediates,
                     _p_data="IMGT2Seq/data", __Nfield_OUTPUT_FORMAT=Nfield_OUTPUT_FORMAT)


        self.setSummaryString() # Inplace summary string setting.
        print(self.summary_string)

        # self.summary_string = \
        #     ''.join([self.summary_string,
        #              "< IMGT2Sequence(Newly generated.) >\n" \
        #              "- HLA Amino Acids : {}\n" \
        #              "- HLA SNPs : {}\n" \
        #              "- HLA Allele Table : {}\n" \
        #              "- Maptables for heatmap : \n" \
        #              "   A   : {A}\n" \
        #              "   B   : {B}\n" \
        #              "   C   : {C}\n" \
        #              "   DPA1: {DPA1}\n" \
        #              "   DPB1: {DPB1}\n" \
        #              "   DQA1: {DQA1}\n" \
        #              "   DQB1: {DQB1}\n" \
        #              "   DRB1: {DRB1}\n".format(self.__dict_AA__, self.__dict_SNPS__, self.__HAT__,
        #                                         **self.__d_MapTable__)
        #              ])

    # def __bool__(self):
    #     return (self.f__dict_AA__ and self.f__dict_SNPS__ and self.f__HAT__ and self.f__d_MapTable__)

    def __str__(self):
        return self.summary_string

    def getResult(self):
        return [self.__dict_AA__, self.__dict_SNPS__, self.__HAT__, self.__d_MapTable__]

    def setSummaryString(self):

        self.summary_string = "< IMGT2Sequence Output >\n"

        str_main_output = \
            "- HLA Amino Acids : {}\n" \
            "- HLA SNPs : {}\n" \
            "- HLA Allele Table : {}\n".format(self.__dict_AA__, self.__dict_SNPS__, self.__HAT__)

        self.summary_string = ''.join([self.summary_string, str_main_output])

        return 0



if __name__ == "__main__":

    """
    Main wrapper for IMGT2Seq.
    """

    parser = argparse.ArgumentParser(prog='IMGT2Seq',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        IMGT2Sequence.py



    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--out", "-o", help="\nOutput File Prefix\n\n", required=True)
    parser.add_argument("--hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)
    parser.add_argument("--imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                        required=True)
    parser.add_argument("--imgt-dir", help="\nIn case User just want to specify the directory of IMGT data folder.\n\n",
                        required=True)

    parser.add_argument("--HLA", help="\nHLA genes to include.\n\n", nargs='+', required=True)

    parser.add_argument("--mp", help="\nSetting off parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    parser.add_argument("--no-Ins", help="\nNo Insertion Markers in output.\n\n", action="store_true")
    parser.add_argument("--include-UTR", help="\nInclude UTR parts(5-prime, 3-prime) in output.\n\n", action="store_true")
    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files. (DEBUG)\n\n", action='store_true')

    # Output format selection
    format_selection = parser.add_mutually_exclusive_group(required=True)
    format_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n", action="store_true", dest="oneF")
    format_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n", action="store_true", dest="twoF")
    format_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n", action="store_true", dest="threeF")
    format_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields.\n\n", action="store_true", dest="fourF")
    format_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n", action="store_true")
    format_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n", action="store_true")



    ##### < for Test > #####

    # hg18 / imgt3320
    # _args = parser.parse_args(["--hg", "18",
    #                            "--out", "tests/20210706_IMGT2Seq/dummy_prefix",
    #                            "--imgt", "3320",
    #                            "--imgt-dir", "example/IMGTHLA3320",
    #                            "--mp", "8"])



    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    ##### < Main function Execution. > #####

    myIMGT2Seq = HATK_IMGT2Seq(args.imgt, args.hg, args.out, args.HLA,
                               args.oneF, args.twoF, args.threeF, args.fourF, args.Ggroup, args.Pgroup,
                               _no_Ins=args.no_Ins, _multiprocess=args.mp,
                               _save_intermediates=args.save_intermediates,
                               _imgt_dir=args.imgt_dir)
