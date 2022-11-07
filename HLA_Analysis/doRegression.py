# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join, isdir

from src.HATK_Error import HATK_InputPreparation_Error, HATK_PLINK_Execution_Error, RaiseError
from src.PLINK_Bash import Bash_RUN_PLINK

from src.PLINK_ASSOC import ASSOC


std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



def doRegression(_plink_exec, _out_prefix, _bfile, _phe, _phe_name, _phe_type=None,
                 _covar=None, _covar_name:str=None, _condition_list=None,
                 _chr=6, _from_mb=29, _to_mb=34, _ci=0.95,
                 _f_hide_covar=True, _f_allow_no_sex=True, _f_keep_allele_order=True) -> ASSOC:

    """
    A wrapper function of the 'doRegression' class method of the 'HLA_Study'.
    """

    command = [_plink_exec]

    if _phe_type == "Continuous": command.append("--linear")
    elif _phe_type == "Binary": command.append("--logistic")
    else:
        RaiseError(
            HATK_PLINK_Execution_Error,
            std_ERROR + "HATK can't decide which regression(Linear/Logistic) to implement. "
                        "Please check the data type of given phenotype('{}') whether it is binary or continuous." \
                        .format(_phe_name))

    if _f_hide_covar: command.append("hide-covar")
    if _f_allow_no_sex: command.append("--allow-no-sex")
    if _f_keep_allele_order: command.append("--keep-allele-order")

    command.append("--bfile " + _bfile)
    command.append("--pheno " + _phe)
    command.append("--pheno-name " + _phe_name)

    if _covar and _covar_name:
        command.append("--covar " + _covar)
        command.append("--covar-name " + _covar_name)

    if _condition_list:
        command.append("--condition-list " + _condition_list)

    command.append("--chr {} --from-mb {} --to-mb {}".format(_chr, _from_mb, _to_mb))
    command.append("--ci {}".format(_ci))

    command.append("--out " + _out_prefix)

    command = ' '.join(command) # command to return.
    # print(command)


    ## PLINK bash execution.
    assoc_result = Bash_RUN_PLINK(command, _out_prefix, _f_save_log=True) # Run PLINK
    regr_type = "linear" if _phe_type == "Continuous" else "logistic"

    return ASSOC(assoc_result+".assoc.{}".format(regr_type), regr_type)



def checkList():

    """
    1. "0/1" or "1/2" in the target phenotype vector.
    2. The same sample number/set in 'GT', 'PHENO', 'COVAR' files.
    3. fam file has phenotype value, too.
    4.
    """

    return 0



if __name__ == '__main__':


    ## Linux

    ## Mac
    _plink = "/Users/wansonchoi/miniconda3/bin/plink"
    _out = "/Users/wansonchoi/Git_Projects/HATK/tests/20221007_bMG/asdf.doRegr"
    _bfile = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged"
    _phe = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.phe"
    _phe_name = "Dis_CD"
    _phe_type = "Binary"
    _covar = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/UC-CD-HLA/data/Merged/merged.covar"
    _covar_name = "GWAS"

    Assoc = doRegression(_plink, _out, _bfile, _phe, _phe_name, _phe_type, _covar, _covar_name)
    print(Assoc)




    pass