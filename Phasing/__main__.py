#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap

from Phasing.BEAGLE_Phase import Phasing_wrapper

std_MAIN_PROCESS_NAME = "\n[Phasing]: "
std_ERROR_MAIN_PROCESS_NAME = "\n[Phasing::ERROR]: "
std_WARNING_MAIN_PROCESS_NAME = "\n[Phasing::WARNING]: "

class HATK_Phasing(object):

    def __init__(self):
        pass



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='Phasing',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        [Phasing]

    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--bfile", help="\nThe prefix for PLINK genotype(SNP) data file(.bed/.bim/.fam)\n\n")

    parser.add_argument("--bgl", help="\nBeagle(v3.x.x) file.\n\n")
    parser.add_argument("--markers", help="\nThe '*.markers' file along with Beagle(v3.x.x) file.\n\n")

    parser.add_argument("--vcf", help="\nThe VCF file.\n\n")

    parser.add_argument("--out", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')

    ##### <for Test> #####

    # 2019. 01. 10
    # args = parser.parse_args(["--variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
    #                           "--out", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)
