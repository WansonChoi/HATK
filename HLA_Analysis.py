# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import src.HLA_Analysis.HLA_Analysis_modules as hla_m


def HLA_Analysis(_bfile, _file, _out,
                 _covar, _covar_names,
                 _phe, _phe_name,
                 _condition, _condition_list,
                 _ref_allele,
                 _lr, _ob, _bt, _meta,
                 _heatmap, _manhattan):


    ########## < Core Varialbes > ##########

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))


    ########## < Argument Checking > ##########



    ########## < Conducting Association Test > ##########

    if _lr:

        print(std_MAIN_PROCESS_NAME + "Conducting Logistic Regression.\n")

        hla_m.__hla__Logistic_Regression(_bfile, _file, _out,
                                         _covar, _covar_names, _phe, _phe_name, _condition, _condition_list, _ref_allele)

    if _ob:

        print(std_MAIN_PROCESS_NAME + "Conducting Omnibus Test.\n")

    if _bt:

        print(std_MAIN_PROCESS_NAME + "Conducting Binary Test.\n")


    ########## < Conducting Meta Analysis > ##########

    if _meta:

        print(std_MAIN_PROCESS_NAME + "Conducting Meta Analysis.\n")

    ########## < Conducting Plotting > ##########

    if _heatmap:

        print(std_MAIN_PROCESS_NAME + "Generating HeatMap plot.\n")

    if _manhattan:

        print(std_MAIN_PROCESS_NAME + "Generating Manhattan plot.\n")


    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###############################################################################       
        
        HLA_Analysis.py
        
        To implment modules for HLA Analysis.
        (ex. OmnibusTest, Manhattan Plotting, etc.)
        
        
        made by Wanson Choi.
        
    ###############################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-i", help="\nInput Data file.\n\n", required=True)
    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)


    # for Publish
    args = parser.parse_args()
    # print(args)

    # for Testing


    HLA_Analysis()