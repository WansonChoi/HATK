# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import src.HLA_Analysis.HLA_Analysis_modules as hla_m



########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))


def ASSOCIATION_TEST(_bfile, _input, _out,
                     _covar, _covar_names,
                     _phe, _phe_name,
                     _condition, _condition_list,
                     _ref_allele, _phased, _threshold,
                     _lr, _ob):


    ### General Argument Checking for Association Tests.

    if bool(_input):

        # When "--input" argument is given.
        print(std_MAIN_PROCESS_NAME + "Necessary input files are given as common prefix(\"--input(-i)\").")

        _bfile = _input
        _covar = _input + ".covar"
        _phe = _input + ".phe"
        _phased = _input + ".aa"

    else:

        if (not bool(_bfile)):
            print(std_ERROR_MAIN_PROCESS_NAME + "You didn't give input file prefix.\n"
                                          "Please check \"--bfile\" option again.\n")
            sys.exit()



    if not bool(_out):
        # When "--out" argument is not given.
        print(std_ERROR_MAIN_PROCESS_NAME + "You didn't give output file prefix. "
                                            "Please check \"-o\" option again.")
        sys.exit()



    if _lr:

        ### Argument Checking for "Logistic Regression".

        if not bool(_ref_allele):
            # Making default reference allele.
            print(std_MAIN_PROCESS_NAME + "Using default reference allele.")
            if bool(_bfile):
                _ref_allele = hla_m.MakeDefaultReferenceAllele(_bfile)



        ### Conducting Logistic Regression.
        print(std_MAIN_PROCESS_NAME + "Conducting Logistic Regression(Plink v1.07).\n")

        hla_m.__hla__Logistic_Regression(_bfile, _input, _out,
                                         _covar, _covar_names, _phe, _phe_name, _condition, _condition_list, _ref_allele)



    if _ob:

        ### Argument Checking for "Omnibus Test".

        if not bool(_condition) and bool(_condition_list):

            # Omnibus doesn't use "--condition-list" but "--condition" only.
            print(std_MAIN_PROCESS_NAME + "Error. Omnibus Test should use \"--condition\" not \"--condition-list\".")
            sys.exit()

        if not bool(_threshold):
            # If user doesn't give "--rare-threshold" option, then set it 0 by default.
            _threshold = "0"



        ### Conducting Omnibus Test.
        print(std_MAIN_PROCESS_NAME + "Conducting Omnibus Test.\n")

        hla_m.__hla__Omnibus_Test(_bfile, _phased, _phe, _covar,
                                  _out, _phe_name, _threshold, _condition)


    return 0



def META_ANALYSIS(_out, *rassoc):

    print(std_MAIN_PROCESS_NAME + "function META_ANALYSIS().")

    # `rassoc` is given as tuple itself.


    ########## < Argument Checking > ##########

    # (1) Output prefix
    if not bool(_out):
        print(std_MAIN_PROCESS_NAME + "Error! You didn't give output file prefix.\n"
                                      "Please check \"-o\" option again.")

    # (2) Existence of input association result files.
    for f in rassoc:

        if not os.path.exists(f):
            print(std_MAIN_PROCESS_NAME + "Error. The file {0} doesn't exist.".format(f))
            sys.exit()

    # (2018. 8. 6.) As far as input association result files have extension defined by plink(ex. .assoc, .fisher, .assoc.logstic, etc.), Meta-analysis can be conducted.
    # If the result of Omnibus test also could be used in meta-analysis maybe...


    ########## < Conducting Meta-Analysis > ##########

    hla_m.__hla__Meta_Analysis(_out, *rassoc)


    return 0




def Pipeline1():


    return 0


def Pipeline2():

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###############################################################################       
        
        HLA_Analysis.py
        
        The modules prepared in "HLA_Analysis_modules.py" will be used in this 
        script. Checking and Preprocessing arguments will be conducted too.
        
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


    # HLA_Analysis()