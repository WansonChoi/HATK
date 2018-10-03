#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))





if __name__ == "__main__":


    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        HATK.py (HLA Analysis Toolkit)

        (Created by Wanson Choi, 2018/05/27)


        <SPECIFICATION> 


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )


    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"
    parser._optionals.description = '\n!-- COMMON ARGUMENTS TO BE USED IN ALL MODULES. --!\n\n'

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n") # Core argument of this program.
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n")

    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="HG")

    parser.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n")
    parser.add_argument("-hped", help="\nHLA Type Data processed by \'Nomencleaner\'(.hped)\n\n")



    ### IMGT2Sequence

    g_IMGT2Sequence = parser.add_argument_group(title='IMGT2Sequence', description='')

    g_IMGT2Sequence.add_argument("--imgt2sequence", help="\nGive this argument to implement \"IMGT2Sequence\" sub-module.\n\n", action='store_true')

    g_IMGT2Sequence.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n")

    g_IMGT2Sequence.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    g_IMGT2Sequence.add_argument("--no-multiprocess", help="\nSetting off parallel multiprocessing.\n\n", action='store_true')



    ### HLA2MARKER

    g_HLA2MARKER = parser.add_argument_group(title='HLA2MARKER', description='')

    g_HLA2MARKER.add_argument("--hla2marker", help="\nGive this argument to implement \"HLA2MARKER\" sub-module.\n\n",
                              action='store_true')

    g_HLA2MARKER.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n")
    g_HLA2MARKER.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")



    ### NomenCleaner

    g_NomenCleaner = parser.add_argument_group(title='NomenCleaner', description='')

    g_NomenCleaner.add_argument("--nomencleaner", help="\nGive this argument to implement \"NomenCleaner\" sub-module.\n\n",
                                action='store_true')

    # Additional input ped file type.
    g_NomenCleaner.add_argument("-ped-Ggroup", help="\nHLA Type Data(G-group allele \"*.ped\" file).\n\n", dest="ped_G")
    g_NomenCleaner.add_argument("-ped-Pgroup", help="\nHLA Type Data(P-group allele \"*.ped\" file).\n\n", dest="ped_P")

    g_NomenCleaner.add_argument("-iat", help="\nIntegrated Allele Table file(*.iat).\n\n")

    # Output format selection
    format_selection = g_NomenCleaner.add_mutually_exclusive_group()


    format_selection.add_argument("--1field", help="\nOutput ped file as '1-field' format.\n\n",
                                        action="store_true", dest="oneF")
    format_selection.add_argument("--2field", help="\nOutput ped file as '2-field' format.\n\n",
                                        action="store_true", dest="twoF")
    format_selection.add_argument("--3field", help="\nOutput ped file as '3-field' format.\n\n",
                                        action="store_true", dest="threeF")
    format_selection.add_argument("--4field",
                                        help="\nOutput ped file as '4-field(Current Standard Names)' format.\n\n",
                                        action="store_true", dest="fourF")
    format_selection.add_argument("--G-group", help="\nOutput ped file as 'G-group' format.\n\n",
                                        action="store_true")
    format_selection.add_argument("--P-group", help="\nOutput ped file as 'P-group' format.\n\n",
                                        action="store_true")

    # Additional utility flags
    g_NomenCleaner.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*\").\n\n", action='store_true')



    ### HLA-Analysis

    g_HLA_Analysis = parser.add_argument_group(title='HLA_Analysis', description='')

    g_HLA_Analysis.add_argument("--hla-analysis", help="\nGive this argument to implement \"HLA-Analysis\" sub-module.\n\n", action='store_true')

    g_HLA_Analysis.add_argument("--covar", help="\nSpecify .covar file (Plink v1.07).\n\n")
    g_HLA_Analysis.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use.(Plink v1.07).\n\n")

    g_HLA_Analysis.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.07).\n\n")
    g_HLA_Analysis.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use.(Plink v1.07).\n\n")

    CondVars = g_HLA_Analysis.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify Marker name(s) to condition as comma-separated(\",\") (Plink v1.07).\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify the Marker name(s) to condition as a file (Plink v1.07).\n\n")

    g_HLA_Analysis.add_argument("--reference-allele", help="\nSpecify the Reference Allele file (Plink v1.07).\n\n")



    ### Plotting

    g_Plotting = parser.add_argument_group(title='Plotting', description='')

    g_Plotting.add_argument("--plotting", help="\nGive this argument to implement \"Plotting\" sub-module.\n\n", action='store_true')

    g_Plotting.add_argument("--heatmap", help="\nGenerate Heatmap Plot.\n\n", action="store_true")
    g_Plotting.add_argument("--manhattan", help="\nGenerate Manhattan Plot.\n\n", action="store_true")



    ### Converter

    g_Converter = parser.add_argument_group(title='Converter', description='')

    g_Converter.add_argument("--converter", help="\nGive this argument to implement \"Converter\" sub-module.\n\n", action='store_true')


    g_Converter.add_argument("--AXIOM", help="\nAXIOM output file format.\n\n", action="store_true")
    g_Converter.add_argument("--HIBAG", help="\nHIBAG output file format.\n\n", action="store_true")
    g_Converter.add_argument("--xHLA", help="\nxHLA output file format.\n\n", action="store_true")









    ##### < for Testing > #####

    # args = parser.parse_args(["--imgt2sequence", "-imgt", "370", "-o", "TEST/TEST", "-hg", "18"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)




    ########## < Argument Chekcing > ##########


    ### Checking indispensable common arguments

    if not bool(args.out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.'.format("--out"))
        sys.exit()

    if not bool(args.input) and not args.imgt2sequence:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.'.format("--input"))
        sys.exit()

    if not (args.imgt2sequence or args.hla2marker or args.nomencleaner or args.hla_analysis or args.plotting or args.converter):
        # if none of flag for sub-module is given, then it is false state.
        print(std_ERROR_MAIN_PROCESS_NAME + "You should give at least on flag for sub-module to use.\n")
        sys.exit()





    if args.imgt2sequence:

        ##### IMGT2Sequence #####

        print(std_MAIN_PROCESS_NAME + "Implementing IMGT2Sequence.")

        """
        List of necessary arguments.

        1. -hg 
        2. -o (*)
        3. -imgt
        """

        if not bool(args.hg):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.'.format("-hg"))
            sys.exit()

        if not bool(args.imgt):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.'.format("-imgt"))
            sys.exit()


        from src.IMGT2Sequence.IMGT2Seqeunce import MakeDictionary

        MakeDictionary(_HG=args.hg, _OUTPUT=args.out, _IMGT=args.imgt,
                       _no_Indel=args.no_indel, _no_MultiP=args.no_multiprocess)



    if args.hla2marker:

        ##### HLA2MARKER #####

        print(std_MAIN_PROCESS_NAME + "Implementing HLA2MARKER.")

        """
        List of necessary arguments.

        1. -input
        2. -hped
        3. -hg
        4. -o
        5. --dict-AA
        6. --dict-SNPS
        """


    if args.nomencleaner:

        ##### NomenCleaner #####

        print(std_MAIN_PROCESS_NAME + "Implementing NomenCleaner.")

        """
        List of necessary arguments.

        1. either -ped, -ped-Ggroup or -ped-Pgroup
        2. -iat
        3. -o
        4. output format(ex. --1field, etc.)

        (optionals)
        5. --No-caption
        """


    elif args.hla_analysis:

        ##### HLA_Analysis #####

        print(std_MAIN_PROCESS_NAME + "Implementing HLA_Analysis.")

        """
        List of necessary arguments.

        1. -input(--bfile)
        2. -o
        3. either -lr or -ob

        (optionals)
        4. -covar (with -covar-name)
        5. -phe (with -phe-name)
        6. either --condition or --condition-list
        7. --reference-allele
        8. threshold
        9. phased.

        """


    if args.plotting:

        ##### Plotting #####

        print(std_MAIN_PROCESS_NAME + "Implementing Plotting.")

        """
        List of necessary arguments.
        """


    if args.converter:

        ##### Converter #####

        print(std_MAIN_PROCESS_NAME + "Implementing Converter.")

        """
        List of necessary arguments.
        """



