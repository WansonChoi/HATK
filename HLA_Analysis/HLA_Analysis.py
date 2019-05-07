# -*- coding: utf-8 -*-

import os, sys, re
from shutil import which
import argparse, textwrap
# import inspect
import pandas as pd



########## < Core Global Varialbes > ##########


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

GLOBAL_p_plink = "./dependency/plink" if os.path.isfile("./dependency/plink") else which("plink")
GLOBAL_p_Rscript = which("Rscript")
GLOBAL_p_JAVA = which("java")



########## < Association Test > ##########


### (Logistic Regression)
def Logistic_Regression(_bfile, _out,
                        _phe=None, _phe_name=None,
                        _covar=None, _covar_names=None,
                        _condition=None, _condition_list=None,
                        _ref_allele=None,  # (2018. 8. 6.) Generate default reference_allele file and make it use that file by default.
                        _from_mb = 29, _to_mb = 34, _hide_covar = True, _ci = 0.95,
                        _chr = 6, _allow_no_sex = True):


    ### Argument Processing && Generating Command.

    command = [GLOBAL_p_plink, "--noweb", "--logistic", "--out {0}".format(_out)]

    if bool(_bfile):
        command.append("--bfile {0}".format(_bfile))

    if bool(_covar):
        command.append("--covar " + _covar)
    if bool(_covar_names):
        command.append("--covar-name " + _covar_names)
    if bool(_phe):
        command.append("--pheno " + _phe)
    if bool(_phe_name):
        command.append("--pheno-name " + _phe_name)

    if bool(_condition):
        command.append("--condition " + _condition)
    elif bool(_condition_list):
        command.append("--condition-list " + _condition_list)

    if bool(_ref_allele):
        # command.append("--reference-allele " + _ref_allele)
        command.append("--a1-allele " + _ref_allele)

    if not (bool(_from_mb) != bool(_to_mb)): # Exclusive-NOR
        command.append("--from-mb {0} --to-mb {1}".format(_from_mb, _to_mb))

    if bool(_hide_covar):
        command.append("--hide-covar")
    if bool(_allow_no_sex):
        command.append("--allow-no-sex")

    if bool(_ci):
        command.append("--ci {0}".format(_ci))
    if bool(_chr):
        command.append("--chr {0}".format(_chr))

    command = ' '.join(command)


    ### Conducting Logistic Regression by Plink(v1.9b).

    print(command)
    os.system(command)


    return (_out + ".assoc.logistic")


def MakeDefaultReferenceAllele(_bfile):

    # step1. Load ".bim" file
    bim = pd.read_table(_bfile+".bim", sep='\t', header=None, usecols=[1,4,5], names=["Label", "Al1", "Al2"])

    l_Al1 = bim.iloc[:, 1].tolist()
    l_Al2 = bim.iloc[:, 2].tolist()

    for i in range(0, bim.shape[0]):

        if (l_Al1[i] == "a" and l_Al2[i] == "p"):
            l_Al1[i] = "p"

    # Making the reference allele DataFrame.

    df_Ref_Allele = pd.concat([bim.iloc[:, 0], pd.Series(l_Al1)], axis=1)
    df_Ref_Allele.to_csv(_bfile+".refallele", sep='\t', header=False, index=False)

    return (_bfile+".refallele")



### (OmnibusTest)
def GetPhasedAlleles(_input, _bgl_phased, _out):


    FAM = _input + ".fam"

    command = ' '.join([GLOBAL_p_Rscript, "src/HLA_Analysis/AssociationTest/AllCC_Get_Phased_AA_Calls.R", _bgl_phased, FAM, _out])
    os.system(command)

    return _out + ".aa"


def Omnibus_Test(_input, _out, _phased, _phe, _phe_name, _covar, _covar_name="NA", _condition="NA"):

    """

    In terms of .covar file, it won't take the column names for selectively using specific columns.
    In other words, if you pass the .covar file to this omnibus test, all columns will be used as covariates. so, if
    you want to use specific columns in input .covar file, preprocess it first and give it to this program as an input.


    """

    ### File existence Check.



    ### Argument Processing && Generating Command.

    command = [GLOBAL_p_Rscript, "src/HLA_Analysis/AssociationTest/OmnibusTest_BHv5.R",
               _out, _input + ".fam", _phased, _phe, _phe_name, _covar, _covar_name]

    if bool(_condition):
        command.append(_condition)
    else:
        command.append("NA")

    command = ' '.join(command)
    print(command)


    if not os.system(command):
        return (_out + ".omnibus")
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Omnibus Test failed.\n")
        sys.exit()










########## < Meta Analysis > ##########


### < Meta Analysis >
def Meta_Analysis(_out, _rassoc):


    ### Generating Command

    command = [GLOBAL_p_plink, "--noweb", "--out {0}".format(_out), "--meta-analysis", *_rassoc]
    command = ' '.join(command)


    # Conducting Meta-Analysis by Plink(v1.07)

    print(command)

    if not os.system(command):
        return (_out + ".meta")
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Meta-Analysis failed.\n")
        sys.exit()








if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###############################################################################
           
        HLA_Analysis_modules.py
        
        A lot of sub-modules used in "HLA_Analysis.py" conduct a bash execution.
        (ex. Rscript, plink, ...)
        
        Those sub-modules are prepared in this script so that bash execution doesn't
        need to be considered anymore. 

        This script is just a collection of those modules.
                
    ###############################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')


    ### Test

    # Logistic_Regression('data/example/LogisticRegression/20190327_WTCCC.AA.CODED', 'tests/20190503_LogisticR',
    #                     _phe='/Users/wansun/Git_Projects/HLA_Analysis/data/example/LogisticRegression/wtccc_filtered_58C_NBS_RA_T1D.phe', _phe_name='RA')