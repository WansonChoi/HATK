# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import src.HLA_Analysis.HLA_Analysis_modules as hla_m



########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


# def ASSOCIATION_TEST(_bfile, _input, _out,
#                      _covar, _covar_names,
#                      _phe, _phe_name,
#                      _condition, _condition_list,
#                      _ref_allele, _phased, _threshold,
#                      _lr, _ob):
#
#
#     ### General Argument Checking for Association Tests.
#
#     if bool(_input):
#
#         # When "--input" argument is given.
#         print(std_MAIN_PROCESS_NAME + "Necessary input files are given as common prefix(\"--input(-i)\").")
#
#         _bfile = _input
#         _covar = _input + ".covar"
#         _phe = _input + ".phe"
#         _phased = _input + ".aa"
#
#     else:
#
#         if (not bool(_bfile)):
#             print(std_ERROR_MAIN_PROCESS_NAME + "You didn't give input file prefix.\n"
#                                           "Please check \"--bfile\" option again.\n")
#             sys.exit()
#
#
#
#     if not bool(_out):
#         # When "--out" argument is not given.
#         print(std_ERROR_MAIN_PROCESS_NAME + "You didn't give output file prefix. "
#                                             "Please check \"-o\" option again.")
#         sys.exit()
#
#
#
#     if _lr:
#
#         ### Argument Checking for "Logistic Regression".
#
#         if not bool(_ref_allele):
#             # Making default reference allele.
#             print(std_MAIN_PROCESS_NAME + "Using default reference allele.")
#             if bool(_bfile):
#                 _ref_allele = hla_m.MakeDefaultReferenceAllele(_bfile)
#
#
#
#         ### Conducting Logistic Regression.
#         print(std_MAIN_PROCESS_NAME + "Conducting Logistic Regression(Plink v1.07).\n")
#
#         hla_m.__hla__Logistic_Regression(_bfile, _input, _out,
#                                          _covar, _covar_names, _phe, _phe_name, _condition, _condition_list, _ref_allele)
#
#
#
#     if _ob:
#
#         ### Argument Checking for "Omnibus Test".
#
#         if not bool(_condition) and bool(_condition_list):
#
#             # Omnibus doesn't use "--condition-list" but "--condition" only.
#             print(std_MAIN_PROCESS_NAME + "Error. Omnibus Test should use \"--condition\" not \"--condition-list\".")
#             sys.exit()
#
#         if not bool(_threshold):
#             # If user doesn't give "--rare-threshold" option, then set it 0 by default.
#             _threshold = "0"
#
#
#
#         ### Conducting Omnibus Test.
#         print(std_MAIN_PROCESS_NAME + "Conducting Omnibus Test.\n")
#
#         hla_m.__hla__Omnibus_Test(_bfile, _phased, _phe, _covar,
#                                   _out, _phe_name, _threshold, _condition)
#
#
#     return 0



def HATK_ASSOC1_Logistic_Regression(_input, _out, _covar=None, _covar_names=None, _phe=None, _phe_name=None,
                                    _condition=None, _condition_list=None, _ref_allele=None, _b_marker_panels=None):


    ##### < Core Variable > #####

    bFILE = None
    COVAR = None
    PHE = None
    REFERENCE_ALLELE = None


    # bed, bim, fam (Main panels for association tests)
    if bool(_b_marker_panels) and \
            (os.path.exists(_b_marker_panels+".bed") and os.path.exists(_b_marker_panels+".bim") and os.path.exists(_b_marker_panels+".fam")):

        bFILE = _b_marker_panels
        REFERENCE_ALLELE = (_b_marker_panels + ".refallele") if os.path.exists(_b_marker_panels + ".refallele") else hla_m.MakeDefaultReferenceAllele(_b_marker_panels)
        # _out = _b_marker_panels

    else:

        if bool(_input) and \
           (os.path.exists(_input+".bed") and os.path.exists(_input+".bim") and os.path.exists(_input+".fam")):

            bFILE = _input
            REFERENCE_ALLELE = (_input + ".refallele") if os.path.exists(_input + ".refallele") else hla_m.MakeDefaultReferenceAllele(_input)

        else:

            print(std_ERROR_MAIN_PROCESS_NAME + "The output panel of \"b:MarkerGenerator\" can't be found. Please check them again.\n")
            sys.exit()


    # _input (*.covar, covar-name, *.phe, etc...)
    if bool(_input):

        COVAR = (_input + ".covar") if os.path.exists(_input + ".covar") else (_input + ".cov") if os.path.exists(_input + ".cov") else None
        PHE = (_input + ".phe") if os.path.exists(_input + ".phe") else (_input + ".pheno") if os.path.exists(_input + ".pheno") else None
        # REFERENCE_ALLELE = (_input + ".refallele") if os.path.exists(_input + ".refallele") else hla_m.MakeDefaultReferenceAllele(_input)

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--input(-i)"))
        sys.exit()


    # _covar
    if bool(_covar):

        # When the argument "--covar" is given. (with or without "--covar-name" argument)
        if os.path.exists(_covar):
            COVAR = (_covar)
        elif os.path.exists(_covar):
            COVAR = (_covar)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The covariate file \"{0}\" given by the argument \"--covar\" doesn't exist. Please check that again.\n")
            sys.exit()

    # _pheno
    if bool(_phe):

        # When the argument "--pheno" is given. (with or without "--pheno-name" argument)
        if os.path.exists(_phe):
            PHE = (_phe)
        elif os.path.exists(_phe):
            PHE = (_phe)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The phenotype file \"{0}\" given by the argument \"--pheno\" doesn't exist. Please check that again.\n")
            sys.exit()


    # _covar_names
    if not bool(COVAR):

        print(std_MAIN_PROCESS_NAME + "No Covariate file has given.\n")

        if bool(_covar_names):
            print(std_WARNING_MAIN_PROCESS_NAME + "Covaraite file han't given but the argument \"--covar-names\" has given. "
                                                  "The argument \"--covar-names\" will be ignored.\n")
            _covar_names = None
    else:
        print(std_MAIN_PROCESS_NAME + "Given covariate file is \"{0}\"".format(COVAR))


    # _pheno_name
    if not bool(PHE):

        print(std_MAIN_PROCESS_NAME + "No Phenotype file has given.\n")

        if bool(_phe_name):
            print(std_WARNING_MAIN_PROCESS_NAME + "Phenotype file han't given but the argument \"--pheno-name\" has given. "
                                                  "The argument \"--pheno-name\" will be ignored.\n")
            _phe_name = None
    else:
        print(std_MAIN_PROCESS_NAME + "Given phenotype file is \"{0}\"".format(PHE))


    # _ref_allele
    if bool(_ref_allele):

        # When the argument "--reference-allele" is given.
        if os.path.exists(_ref_allele):
            REFERENCE_ALLELE = _ref_allele
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The reference allele file \"{0}\" given by the argument \"--reference-allele\" doesn't exist. Please check that again.\n")
            sys.exit()

    if not bool(REFERENCE_ALLELE):
        print(std_MAIN_PROCESS_NAME + "No Reference Allele file has given.\n")
    else:
        print(std_MAIN_PROCESS_NAME + "Given Reference Allele file is \"{0}\"".format(REFERENCE_ALLELE))


    # _codition and _condition_list
    if bool(_condition) and bool(_condition_list):
        print(std_ERROR_MAIN_PROCESS_NAME + "Either \"{0}\" or \"{1}\" should be given.\n".format("--condition", "--conditino-list"))
        sys.exit()

    elif not bool(_condition) and bool(_condition_list):

        # When _condition_list is given.
        if not os.path.isfile(_condition_list):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given file of the argument \"{0}\" doesn't exist. Please check it again.\n".format(_condition_list))
            sys.exit()


    # _out
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out(-o)"))
        sys.exit()
    else:
        # Intermediate path.
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        INTERMEDIATE_PATH = os.path.dirname(_out)

        if not os.path.exists(INTERMEDIATE_PATH):
            os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))



    print(std_MAIN_PROCESS_NAME + "Conducting Logistic Regression for \"{0}\" (Plink v1.07).\n".format(bFILE))

    return hla_m.__hla__Logistic_Regression(_bfile=bFILE, _out=_out, _covar=COVAR, _covar_names=_covar_names,
                                            _phe=PHE, _phe_name=_phe_name, _condition=_condition, _condition_list=_condition_list,
                                            _ref_allele=REFERENCE_ALLELE)






def ASSOC2_Omnibus_Test(_input, _out, _phased, _phe, _phe_name, _covar, _covar_names, _threshold, _condition, _condition_list):



    ##### < Core Variable > #####

    PHE = None
    COVAR = None
    PHASED = None


    # _input
    if bool(_input):
        COVAR = (_input + ".covar") if os.path.exists(_input + ".covar") else (_input + ".cov") if os.path.exists(_input + ".cov") else None
        PHE = (_input + ".phe") if os.path.exists(_input + ".phe") else (_input + ".pheno") if os.path.exists(_input + ".pheno") else None

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--input(-i)"))
        sys.exit()

    # _covar
    if bool(_covar):

        # When the argument "--covar" is given. (with or without "--covar-name" argument)
        # Change `COVAR` from the one given by "_input" argument to another one given by "_covar" argument.
        if os.path.exists(_covar):
            COVAR = (_covar)
        elif os.path.exists(_covar):
            COVAR = (_covar)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The covariate file \"{0}\" given by the argument \"--covar\" doesn't exist. Please check that again.\n")
            sys.exit()

    # _pheno
    if bool(_phe):

        # When the argument "--pheno" is given. (with or without "--pheno-name" argument)
        if os.path.exists(_phe):
            PHE = (_phe)
        elif os.path.exists(_phe):
            PHE = (_phe)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The phenotype file \"{0}\" given by the argument \"--pheno\" doesn't exist. Please check that again.\n")
            sys.exit()

    # _covar_names
    if not bool(COVAR):

        # Actually, this part isn't prepared for OmnibusTest.R

        print(std_MAIN_PROCESS_NAME + "No Covariate file has given.\n")

        if bool(_covar_names):
            print(std_WARNING_MAIN_PROCESS_NAME + "Covaraite file han't given but the argument \"--covar-names\" has given. "
                                                  "The argument \"--covar-names\" will be ignored.\n")
            _covar_names = None
    else:
        print(std_MAIN_PROCESS_NAME + "Given covariate file is \"{0}\"".format(COVAR))

    # _pheno_name
    if not bool(PHE):

        print(std_MAIN_PROCESS_NAME + "No Phenotype file has given.\n")

        if bool(_phe_name):
            print(std_WARNING_MAIN_PROCESS_NAME + "Phenotype file han't given but the argument \"--pheno-name\" has given. "
                                                  "The argument \"--pheno-name\" will be ignored.\n")
            _phe_name = None
    else:
        print(std_MAIN_PROCESS_NAME + "Given phenotype file is \"{0}\"".format(PHE))

    # _codition and _condition_list
    if bool(_condition) and bool(_condition_list):
        print(std_ERROR_MAIN_PROCESS_NAME + "Either \"{0}\" or \"{1}\" should be given.\n".format("--condition", "--conditino-list"))
        sys.exit()

    elif not bool(_condition) and bool(_condition_list):

        # When _condition_list is given.
        if not os.path.isfile(_condition_list):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given file of the argument \"{0}\" doesn't exist. Please check it again.\n".format(_condition_list))
            sys.exit()

    # _out
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out(-o)"))
        sys.exit()
    else:
        # Intermediate path.
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        INTERMEDIATE_PATH = os.path.dirname(_out)

        if not os.path.exists(INTERMEDIATE_PATH):
            os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))


    # _phased
    if bool(_phased):
        if os.path.exists(_phased):
            PHASED = _phased
        else:
            PHASED = hla_m.__hla__GetPhasedAlleles(_input, _out)
            ### (2018. 10. 11.) 여기 잠깐 보류해놓음. 나중에 각 implementation wrapper 각각 도입하고 여기좀 다시 돌아올것.



    return hla_m.__hla__Omnibus_Test





def META_ANALYSIS(_out, *rassoc):

    print(std_MAIN_PROCESS_NAME + "function META_ANALYSIS().")


    ##### < Argument Checking > #####

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