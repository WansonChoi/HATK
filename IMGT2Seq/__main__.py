import os, sys, re
import argparse, textwrap
from IMGT2Seq.IMGT2Seq import IMGT2Seq


########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
raw_HLA_names = ["A", "B", "C", "DPA", "DPB", "DQA", "DQB", "DRB"]



class HATK_IMGT2Seq(object):

    def __init__(self, _imgt, _hg, _out, *args, **kwargs):

        """

        Either (1) generate a new Dictionary files or (2) find existing dictionaries.

        """

        if not _imgt:
            print(std_ERROR_MAIN_PROCESS_NAME + "IMGT version information wasn't given.\n"
                                                "Please check '-imgt' argument again.")
            sys.exit()

        if not _hg:
            print(std_ERROR_MAIN_PROCESS_NAME + "HG(Human Genome) version information wasn't given.\n"
                                                "Please check '-hg' argument again.")
            sys.exit()

        if not kwargs['_imgt_dir']:
            print(std_ERROR_MAIN_PROCESS_NAME + "Path to IMGT/HLA raw data directory wasn't given.\n"
                                                "Please check '--imgt-dir' argument again.")
            sys.exit()
        else:
            if not os.path.isdir(kwargs['_imgt_dir']):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given directory('{}') for IMGT-HLA raw data is not a directory.\n"
                                                    "Please check '--imgt-dir' argument again.".format(
                    kwargs['_imgt_dir']))
                sys.exit()

        ########## < Assigning arguments > ##########

        self.imgt = _imgt
        self.hg = _hg
        self.out = _out

        # Main Outputs
        self.__dict_AA__ = None
        self.__dict_SNPS__ = None
        self.__HAT__ = None
        self.__d_MapTable__ = None

        # Flags to check
        self.f__dict_AA__ = False
        self.f__dict_SNPS__ = False
        self.f__HAT__ = False
        self.f__d_MapTable__ = False

        # Optional arguments
        self.no_indel = kwargs["_no_indel"]
        self.multiprocess = kwargs["_multiprocess"]
        self.save_intermediates = kwargs["_save_intermediates"]
        self.imgt_dir = kwargs["_imgt_dir"]

        # Summary string
        self.summary_string = ""

        dir_path = os.path.dirname(self.out)
        version_label = "hg{}.imgt{}".format(self.hg, self.imgt)

        ########## < Checking existing results > ##########

        # (1) dict_AA
        dict_AA_txt = os.path.join(dir_path, "HLA_DICTIONARY_AA.{}.txt".format(version_label))
        dict_AA_map = os.path.join(dir_path, "HLA_DICTIONARY_AA.{}.map".format(version_label))

        if os.path.exists(dict_AA_txt) and os.path.exists(dict_AA_map):
            self.__dict_AA__ = os.path.join(dir_path, "HLA_DICTIONARY_AA.{}".format(version_label))
            self.f__dict_AA__ = True

        # (2) dict_SNPS
        dict_SNPS_txt = os.path.join(dir_path, "HLA_DICTIONARY_SNPS.{}.txt".format(version_label))
        dict_SNPS_map = os.path.join(dir_path, "HLA_DICTIONARY_SNPS.{}.map".format(version_label))

        if os.path.exists(dict_SNPS_txt) and os.path.exists(dict_SNPS_map):
            self.__dict_SNPS__ = os.path.join(dir_path, "HLA_DICTIONARY_SNPS.{}".format(version_label))
            self.f__dict_SNPS__ = True

        # (3) hat
        hat = os.path.join(dir_path, "HLA_ALLELE_TABLE.imgt{}.hat".format(self.imgt))

        if os.path.exists(hat):
            self.__HAT__ = hat
            self.f__HAT__ = True

        # (4) maptables
        self.__d_MapTable__ = {hla_name: None for hla_name in HLA_names}

        f_temp = True

        for hla_name in HLA_names:
            t_maptable = os.path.join(dir_path, "HLA_MAPTABLE_{}.{}.txt".format(hla_name, version_label))
            if os.path.exists(t_maptable):
                self.__d_MapTable__[hla_name] = t_maptable
            else:
                f_temp = False

        self.f__d_MapTable__ = f_temp

        # Output Field format
        oneF = args[0]
        twoF = args[1]
        threeF = args[2]
        fourF = args[3]
        Ggroup = args[4]
        Pgroup = args[5]

        Nfield_OUTPUT_FORMAT = 1 if oneF else 2 if twoF else 3 if threeF else 4

        if oneF:
            Nfield_OUTPUT_FORMAT = 1
        elif twoF:
            Nfield_OUTPUT_FORMAT = 2
        elif threeF:
            Nfield_OUTPUT_FORMAT = 3
        else:
            Nfield_OUTPUT_FORMAT = 4

            if Ggroup or Pgroup:
                print(std_WARNING_MAIN_PROCESS_NAME + "Given '--{}' argument will be overridden to '--4field'.".format(
                    'Ggroup' if Ggroup else 'Pgroup'))



        ### No more using existing result.
        self.__dict_AA__, self.__dict_SNPS__, self.__HAT__, self.__d_MapTable__ = \
            IMGT2Seq(self.imgt, self.hg, self.out, _imgt_dir=self.imgt_dir, _no_Indel=self.no_indel,
                     _MultiP=self.multiprocess, _save_intermediates=self.save_intermediates,
                     _p_data="IMGT2Seq/data", __Nfield_OUTPUT_FORMAT=Nfield_OUTPUT_FORMAT)

        self.f__dict_AA__ = True
        self.f__dict_SNPS__ = True
        self.f__HAT__ = True
        self.f__d_MapTable__ = True

        self.summary_string = \
            ''.join([self.summary_string,
                     "< IMGT2Sequence(Newly generated.) >\n" \
                     "- HLA Amino Acids : {}\n" \
                     "- HLA SNPs : {}\n" \
                     "- HLA Allele Table : {}\n" \
                     "- Maptables for heatmap : \n" \
                     "   A   : {A}\n" \
                     "   B   : {B}\n" \
                     "   C   : {C}\n" \
                     "   DPA1: {DPA1}\n" \
                     "   DPB1: {DPB1}\n" \
                     "   DQA1: {DQA1}\n" \
                     "   DQB1: {DQB1}\n" \
                     "   DRB1: {DRB1}\n".format(self.__dict_AA__, self.__dict_SNPS__, self.__HAT__,
                                                **self.__d_MapTable__)
                     ])

    def __bool__(self):
        return (self.f__dict_AA__ and self.f__dict_SNPS__ and self.f__HAT__ and self.f__d_MapTable__)

    def __str__(self):
        return self.summary_string

    def getResult(self):
        return [self.__dict_AA__, self.__dict_SNPS__, self.__HAT__, self.__d_MapTable__]





if __name__ == "__main__":

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
    parser.add_argument("--hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)
    parser.add_argument("--out", "-o", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser.add_argument("--imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                        required=True)

    parser.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    parser.add_argument("--multiprocess", "--mp", help="\nSetting off parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files. (DEBUG)\n\n", action='store_true')
    parser.add_argument("--imgt-dir", help="\nIn case User just want to specify the directory of IMGT data folder.\n\n", required=True)

    # Output format selection
    format_selection = parser.add_mutually_exclusive_group()
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

    _args = parser.parse_args()
    print(_args)


    ##### < Main function Execution. > #####

    myIMGT2Seq = HATK_IMGT2Seq(_args.imgt, _args.hg, _args.out,
                               _args.oneF, _args.twoF, _args.threeF, _args.fourF, _args.Ggroup, _args.Pgroup,
                               _no_indel=_args.no_indel, _multiprocess=_args.multiprocess,
                               _save_intermediates=_args.save_intermediates,
                               _imgt_dir=_args.imgt_dir)
