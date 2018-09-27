#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap


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

    # parser._optionals.description = '\n\"\"\" Common arguments to be used in all modules. \"\"\"\n\n'

    # parser._optionals.description =     \
    #     '---------------------------------------------------------------------------------\n' \
    #     '- Common arguments to be used in all modules.\n' \
    #     '---------------------------------------------------------------------------------\n'

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n") # Core argument of this program.
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n")

    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")

    parser.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n")
    parser.add_argument("-hped", help="\nHLA Type Data processed by \'Nomencleaner\'(.hped)\n\n")



    ### MakeDictionary

    g_MakeDictionary = parser.add_argument_group(title='MakeDictionary', description='')

    g_MakeDictionary.add_argument("--make-dictionary", help="\nGive this argument to implement \"MakeDictionary\" function.\n\n", action='store_true')
    g_MakeDictionary.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n")



    ### HLA2MARKER

    g_HLA2MARKER = parser.add_argument_group(title='HLA2MARKER', description='')

    g_HLA2MARKER.add_argument("--hla2marker", help="\nGive this argument to implement \"HLA2MARKER\" function.\n\n",
                              action='store_true')

    g_HLA2MARKER.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", default="Not_given")
    g_HLA2MARKER.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", default="Not_given")



    ### NomenCleaner

    g_NomenCleaner = parser.add_argument_group(title='NomenCleaner', description='')

    g_NomenCleaner.add_argument("--nomencleaner", help="\nGive this argument to implement \"NomenCleaner\" function.\n\n",
                                action='store_true')

    # Additional input ped file type.
    g_NomenCleaner.add_argument("-ped-Ggroup", help="\nHLA Type Data(G-group allele \"*.ped\" file).\n\n", dest="ped_G",
                          default="Not_given")
    g_NomenCleaner.add_argument("-ped-Pgroup", help="\nHLA Type Data(P-group allele \"*.ped\" file).\n\n", dest="ped_P",
                          default="Not_given")

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

    g_HLA_Analysis.add_argument("--hla-analysis", help="\nGive this argument to implement \"HLA-Analysis\" function.\n\n", action='store_true')

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

    g_Plotting.add_argument("--plotting", help="\nGive this argument to implement \"Plotting\" function.\n\n", action='store_true')

    g_Plotting.add_argument("--heatmap", help="\nGenerate Heatmap Plot.\n\n", action="store_true")
    g_Plotting.add_argument("--manhattan", help="\nGenerate Manhattan Plot.\n\n", action="store_true")



    ### Converter

    g_Converter = parser.add_argument_group(title='Converter', description='')

    g_Converter.add_argument("--converter", help="\nGive this argument to implement \"Converter\" function.\n\n", action='store_true')


    g_Converter.add_argument("--AXIOM", help="\nAXIOM output file format.\n\n", action="store_true")
    g_Converter.add_argument("--HIBAG", help="\nHIBAG output file format.\n\n", action="store_true")
    g_Converter.add_argument("--xHLA", help="\nxHLA output file format.\n\n", action="store_true")









    ##### < for Testing > #####

    # args = parser.parse_args(["--make-dictionary", "-imgt", "370", "-o", "TEST/TEST", "-hg", "18"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)



    ### Argument Chekcing


    # # [Error] : When more than 1 function flag is given.
    # if sum([int(args.)])
    #
    # _which_function = -1



