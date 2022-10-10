#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join
import subprocess as sbp
from shutil import which
import argparse, textwrap

from src.util import Exists, checkFile
from src.HATK_Error import HATK_InputPreparation_Error, HATK_R_Execution_Error, RaiseError
from src.PLINK import PHENO, COVAR, CONDITION

std_MAIN = "\n[OmnibusTest]: "
std_ERROR = "\n[OmnibusTest::ERROR]: "
std_WARNING = "\n[OmnibusTest::WARNING]: "


class HATK_OmibusTest():

    # External software
    Rscript = checkFile(which("Rscript"))

    def __init__(self, _out_prefix, _bfile, _pheno, _pheno_name:list, _aa=None, _bgl_phased=None,
                 _covar=None, _covar_name=None, _condition=None, _condition_list=None,
                 _f_save_intermediates=False):

        ### Main Variables ###
        self.out_prefix = _out_prefix
        self.out_dir = dirname(self.out_prefix)
        self.fam = checkFile(str(_bfile) +'.fam', std_ERROR + "The given fam file('{}') can't be found.".format(_bfile + '.fam'))

        self.aa = _aa
        self.bgl_phased = _bgl_phased

        self.PHENO = PHENO(_pheno, _pheno_name)
        self.COVAR = COVAR(_covar, _covar_name) if _covar and _covar_name else None
        self.CONDITION = CONDITION(_condition, _condition_list) if bool(_condition) != bool(_condition_list) else None

        self.OMNIBUS = None # Result.


        ### Main Actions ###
        # print(self.__repr__())

        # Phenotype check.
        if self.PHENO.N_pheno_name_target != 1:
            if self.PHENO.N_pheno_name_target == 0:
                str_temp = "Given phenotype('{}') is NOT IN given phenotype file('{}').".format(_pheno_name, _pheno)
            else:
                str_temp = "The Omnibus test takes only one phenotype name per test. ('{}')".format(self.PHENO.pheno_name_target)

            raise HATK_InputPreparation_Error(std_ERROR + str_temp)

        # '*.aa'
        if (not self.aa) and self.bgl_phased: # '*.aa' not given.
            self.aa = getPhasedAACalls(self.fam, self.bgl_phased, self.out_prefix, self.Rscript, _f_save_intermediates)
        if not (self.aa or self.bgl_phased): # Both not given.
            RaiseError(HATK_InputPreparation_Error,
                       std_ERROR + "Phased information is required. Please check '--bgl-phased' and '--aa' arguments again.")

        # Omnibus Test
        self.OMNIBUS = Omnibus_Test(self.out_prefix, self.fam, self.aa, self.PHENO, self.COVAR, self.CONDITION,
                                    self.Rscript, _f_save_intermediates)
        print(self.__repr__())


    def __repr__(self):
        str_main = "\n< Summary of the Omnibus Test >\n"
        str_out_prefix = \
            "- Output prefix: {}\n".format(self.out_prefix)
        str_fam = \
            "- PLINK fam file: {}\n".format(self.fam)
        str_bgl_phased = \
            "- Phased beagle file: {}\n".format(self.bgl_phased)
        str_aa = \
            "- Preprocessed phased beagle file: {}\n".format(self.aa)
        str_OMNIBUS = \
            "- Omnibus test Result: {}\n".format(self.OMNIBUS)

        str_PHENO = \
            "[ PLINK Phenotype file ]\n{}\n".format(self.PHENO)
        str_COVAR = \
            "[ PLINK Covariate file ]\n{}\n".format(self.COVAR) if self.COVAR else ""
        str_CONDITION = \
            "[ PLINK Condition file ]\n{}\n".format(self.CONDITION) if self.CONDITION else ""


        str_summary = ''.join([str_main, str_out_prefix, str_fam, str_bgl_phased, str_aa, str_OMNIBUS,
                               str_PHENO, str_COVAR, str_CONDITION]).rstrip('\n')
        return str_summary



def getPhasedAACalls(_fam, _bgl_phased, _out_prefix, _Rscript, _f_save_intermediates=False):

    print(std_MAIN + "Preprocessing given phased beagle file('{}').".format(_bgl_phased))

    command = "{Rscript} {src} {bgl_phased} {fam} {out}" \
                .format(Rscript=_Rscript, src="OmnibusTest/src/AllCC_Get_Phased_AA_Calls.R", bgl_phased=_bgl_phased,
                        fam=_fam, out=_out_prefix)
    # print(command)

    try:
        with open(_out_prefix + '.aa.log', 'w') as f_aa_log:
            sbp.run(command.split(), check=True, stdout=f_aa_log, stderr=f_aa_log)

    except sbp.CalledProcessError:
        raise HATK_R_Execution_Error(
            std_ERROR + "Failed to preprocess the phased beagle file('{}').\n".format(_bgl_phased) +
            "Please check the log file('{}').".format(_out_prefix + '.aa.log')
        )
    else:
        if Exists(_out_prefix + '.aa.log') and not _f_save_intermediates: os.remove(_out_prefix + '.aa.log')
        return _out_prefix+'.aa'



def Omnibus_Test(_out_prefix, _fam, _aa, _PHENO: PHENO, _COVAR: COVAR, _COND: CONDITION, _Rscript=which("Rscript"),
                 _f_save_intermediates=False):

    print(std_MAIN + "Performing Omnibus Test.")

    # Phenotype info.
    _phe = _PHENO.phe
    _phe_name = _PHENO.pheno_name_target[0]
    # Covariate info.
    if _COVAR:
        _covar = _COVAR.phe
        _covar_name = _COVAR.covar_name_target
    else:
        _covar = "NA"
        _covar_name = "NA"
    # Condition info.
    if _COND:
        _condition = ','.join(_COND.l_condition)
        _out = _out_prefix + ".{}.omnibus".format(_condition)
    else:
        _condition = "NA"
        _out = _out_prefix + ".omnibus"


    command = "{Rscript} {src} {out} {fam} {aa} {phe} {phe_name} {covar} {covar_name} {condition}" \
        .format(Rscript=_Rscript, src="OmnibusTest/OmnibusTest_BHv5.R", out=_out, fam=_fam, aa=_aa, phe=_phe,
                phe_name=_phe_name, covar=_covar, covar_name=_covar_name, condition=_condition)
    print(command)

    try:
        with open(_out + '.log', 'w') as f_OM_log:
            sbp.run(command.split(), check=True, stdout=f_OM_log, stderr=f_OM_log)

    except sbp.CalledProcessError:
        raise HATK_R_Execution_Error(
            std_ERROR + "Failed to perform the Omnibus Test.\n" +
            "Please check the log file('{}').".format(_out+'.log')
        )
    else:
        if Exists(_out+'.log') and not _f_save_intermediates: os.remove(_out+'.log')

        print(std_MAIN.lstrip('\n') + "Omnibus test successfully done.")
        return _out



if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='OmnibusTest',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        Omnibus Test

        - The likelihood ratio test to a amino acid position with several residues.
        
        
    #################################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--out", help="\nOutput file prefix\n\n", required=True)
    parser.add_argument("--bfile", help="\nThe prefix for PLINK genotype(SNP) data file(.bed/.bim/.fam)\n\n", required=True)

    Phased = parser.add_mutually_exclusive_group()
    Phased.add_argument("--bgl-phased", help="The phased beagle file(ex. \"*.bgl.phased\") to make \"*.aa\" file.\n\n")
    Phased.add_argument("--aa", help="\nThe \"*.aa\" which is processed from the phased beagle file('*.bgl.phased).\n\n")

    parser.add_argument("--covar", help="\nSpecify .covar file.\n\n")
    parser.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use.\n\n",
                        nargs='+')

    parser.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.9).\n\n")
    parser.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use.\n\n",
                        nargs='+')

    CondVars = parser.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify a single variant ID to condition(i.e. To set it as covariate).\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify a tsv file of multiple variant IDs to condition(i.e. To set it as covariates). (Plink v1.9)\n\n")

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')

    ##### <for Test> #####

    # 2019. 01. 10
    # args = parser.parse_args([
    #     "--out", "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/remove/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.tt.OM",
    #     "--bfile", "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset",
    #     "--bgl-phased", "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/remove/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.tt.bgl.phased",
    #     "--pheno", "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.chr6.hg18.29-34mb.phe",
    #     "--pheno-name", "RA"
    # ])

    """
    python -m OmnibusTest \
        --out '/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/remove/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.tt.OM' \
        --bfile '/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset'  \
        --bgl-phased '/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/remove/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.tt.bgl.phased' \
        --pheno '/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.chr6.hg18.29-34mb.phe' \
        --pheno-name RA
    """

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    HATK_OmibusTest(args.out, args.bfile, args.pheno, args.pheno_name, args.aa, args.bgl_phased,
                    args.covar, args.covar_name, args.condition, args.condition_list,
                    args.save_intermediates)