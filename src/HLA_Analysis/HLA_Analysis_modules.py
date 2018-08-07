# -*- coding: utf-8 -*-

import os, sys, re
from shutil import which
import argparse, textwrap
# import inspect
import pandas as pd



########## < Core Global Varialbes > ##########


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))

GLOBAL_p_plink = "./dependency/plink" if os.path.isfile("./dependency/plink") else which("plink")
GLOBAL_p_Rscript = which("Rscript")
GLOBAL_p_JAVA = which("java")



########## < Association Test > ##########


### (Logistic Regression)
def __hla__Logistic_Regression(_bfile, _file, _out,
                               _covar, _covar_names,
                               _phe, _phe_name,
                               _condition, _condition_list,
                               _ref_allele, # (2018. 8. 6.) Generate default reference_allele file and make it use that file by default.
                               _from_mb = 29, _to_mb = 34, _hide_covar = True, _ci = 0.95,
                               _chr = 6, _allow_no_sex = True):


    ### Argument Processing && Generating Command.

    command = [GLOBAL_p_plink, "--noweb", "--logistic", "--out {0}".format(_out)]

    if bool(_bfile):
        command.append("--bfile {0}".format(_bfile))
    elif bool(_file):
        command.append("--file {0}".format(_file))

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
        command.append("--reference-allele " + _ref_allele)

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


    ### Conducting Logistic Regression by Plink(v1.07).

    print(command)
    os.system(command)


    return 0


def MakeDefaultReferenceAllele(_bfile):

    # step1. Load ".bim" file
    bim = pd.read_table(_bfile+".bim", sep='\t', header=None, usecols=[1,4,5], names=["Label", "Al1", "Al2"])

    l_Al1 = bim.iloc[:, 1].tolist()
    l_Al2 = bim.iloc[:, 2].tolist()

    for i in range(0, bim.shape[0]):

        if (l_Al1[i] == "A" and l_Al2[i] == "P"):
            l_Al1[i] = "P"

    # Making the reference allele DataFrame.

    df_Ref_Allele = pd.concat([bim.iloc[:, 0], pd.Series(l_Al1)], axis=1)
    df_Ref_Allele.to_csv(_bfile+".refallele", sep='\t', header=False, index=False)

    return (_bfile+".refallele")



### (OmnibusTest)
def __hla__Omnibus_Test():

    return 0







########## < Meta Analysis > ##########


### < Meta Analysis >
def __hla__Meta_Analysis(_out, *assocs):


    ### Generating Command

    command = [GLOBAL_p_plink, "--noweb", "--out {0}".format(_out), "--meta-analysis"]
    command.extend(assocs)

    command = ' '.join(command)


    # Conducting Meta-Analysis by Plink(v1.07)

    print(command)
    os.system(command)

    return 0



########## < Plotting > ##########

### < heatmap >
def __hla__Heatmap():

    return 0

### < SNP Manhattan Plot >
def __hla__Manhattan_Plot():

    return 0

### < AA Manhanttan Plot >
def __hla__AA_Manhattan_Plot():

    return 0



"""

Actually, Only functions to conduct bash execution of external modules are supposed to be collected here,
below "if __name__ == "__main__": ~" aren't necessary.


"""


def main():

    std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))

    return 0




if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###############################################################################
           
        HLA_Analysis_modules.py
        
        A lot of modules used in "HLA_Analysis.py" is coded by other programming 
        languages(ex. R).
        
        The functions that conduct a bash execution of those external modules are
        prepared in this script. These functions are to be used in "HLA_Analysis.py".  
        
        
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

    # for Testing

    # print(args)

    main()