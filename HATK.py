#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename
import argparse, textwrap

# from src.ArgRouter import ArgRouter
from src.implementAll import implementAll

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



if __name__ == "__main__":


    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        HATK.py (HLA Analysis Toolkit)

        (Created by Wanson Choi, 2018/05/27)


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )


    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"
    parser._optionals.description = '\n!-- COMMON ARGUMENTS TO BE USED IN ALL MODULES. --!\n\n'

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    # parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n") # Core argument of this program.
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True) # (***) required

    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="HG",
                        required=True) # (***) required

    # parser.add_argument("--HLA", help="\nHLA genes for association test.\n\n",
    #                     default=['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files. (DEBUG)\n\n", action='store_true') # It will be shared with 'b:MarkerGenerator'


    # PLINK related.
    parser.add_argument("--bfile", help="\nNormal SNP marker data. (Plink v1.9)\n\n", required=True) # Formerly, "--variants".

    parser.add_argument("--covar", help="\nCovariate file(e.g. '*.covar') (Plink v1.9).\n\n")
    parser.add_argument("--covar-name", help="\nThe column name(s) to use in .covar file for Association test. (Plink v1.9)\n\n")

    parser.add_argument("--pheno", help="\nPhenotype file(e.g. '*.phe' (Plink v1.9).\n\n", required=True)
    parser.add_argument("--pheno-name", help="\nA single column name for association test. (Plink v1.9)\n\n", required=True)
    parser.add_argument("--pheno-name-dtype", help="\nData type(either 'Continuous' or 'Binary') of the '--pheno-name' argument.\n\n",
                        choices=["Continuous", "Binary"])

    CondVars = parser.add_mutually_exclusive_group()
    # CondVars.add_argument("--condition", help="\nSpecify variant ID(s) to condition(i.e. To set it/them as covariate(s)). (Plink v1.9)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify a file of multiple variant IDs to condition(i.e. To set it as covariates). (Plink v1.9)\n\n")

    # Output format selection
    format_selection = parser.add_mutually_exclusive_group()
    format_selection.add_argument("--1field", help="\nHLA alleles will have maximum 1 field.\n\n", action="store_true", dest="F_one")
    format_selection.add_argument("--2field", help="\nHLA alleles will have maximum 2 fields.\n\n", action="store_true", dest="F_two") # default format. (2022.10.12.)
    format_selection.add_argument("--3field", help="\nHLA alleles will have maximum 3 fields.\n\n", action="store_true", dest="F_three")
    format_selection.add_argument("--4field", help="\nHLA alleles will have maximum 4 fields.\n\n", action="store_true", dest="F_four")
    format_selection.add_argument("--Ggroup", help="\nHLA alleles will have G code names.\n\n", action="store_true", dest="F_Ggroup")
    format_selection.add_argument("--Pgroup", help="\nHLA alleles will have P code names.\n\n", action="store_true", dest="F_Pgroup")



    ##### < IMGT2Sequence > #####
    parser.add_argument("--imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                        required=True)
    parser.add_argument("--multiprocess", help="\n# of CPU cores to use. (Mainly for IMGT2Seq.)\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    parser.add_argument("--imgt-dir", help="\nThe path for downloaded IMGT-HLA database folder.\n\n", required=True)



    ##### < NomenCleaner > #####
    parser.add_argument("--rhped", help="\nOutput result(s) from other HLA software to be processed by \'HLA2HPED\'.\n\n", nargs='*')
    HLA_type = parser.add_mutually_exclusive_group()
    HLA_type.add_argument("--hped", help="\nHLA Type Data not yet processed by the 'NomenCleaner'. (*.hped)\n\n")
    HLA_type.add_argument("--chped", help="\nHLA Type Data processed by the 'NomenCleaner' (*.chped)\n\n")

    # Additional utility flags
    # parser.add_argument("--NoCaption", help="\nMake converted HLA alleles NOT have HLA gene prefix(ex. \"A*\").\n\n", action='store_true')
    # parser.add_argument("--leave-NotFound", help="\nLeaving HLA alleles which can't be found in given *.hat file(Novel or Erroneous allele) intact.\n\n", action='store_true')



    ##### < Manhattan Plot > #####
    # parser.add_argument("--point-color", "-pc", help="\nPoint color(ex. \"#778899\").\n\n", default="#778899")
    # parser.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n", default="#FF0000")
    #
    # parser.add_argument("--point-size", "-ps", help="\nGeneral point size (default: 15).\n\n", default="15")
    # parser.add_argument("--yaxis-unit", "-yau", help="\nY-axis value(-log10(x)) unit (default : 10).\n\n", default="10")





    ##### < for Testing > #####

    # args = parser.parse_args('--omnibus --out /Users/wansun/Git_Projects/HATK/tests/_4_HLA_Analysis/Omnibus/merged/20190618_merged_covar --input /Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/data/HLA_Analysis/sample/Merged/merged --covar-name GWAS --pheno-name All_UC'.split(' '))


    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)

    hStudy_all = implementAll(args.out, args.hg, args.imgt, args.imgt_dir, args.hped, args.bfile,
                              args.pheno, args.pheno_name, args.pheno_name_dtype,
                              args.covar, args.covar_name, None, args.condition_list,
                              args.F_one, args.F_two, args.F_three, args.F_four, args.F_Ggroup, args.F_Pgroup,
                              args.multiprocess, args.save_intermediates)

    print(hStudy_all)
