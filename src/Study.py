#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join
import subprocess
from shutil import which
import pandas as pd
import numpy as np

import src.HATK_Error as HATK_Error
import src.PLINK as PLINK
from HLA_Manhattan.manhattan import Manhattan

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)



class Study(object):

    """
    Wrapper class to point input PLINK files and perform tests.

    """

    ### Class global variables
    plink_exec = which("plink") if exists(which("plink")) else  \
                "./dependency/plink" if exists("./dependency/plink") else -1
    if plink_exec == -1:
        HATK_Error.HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "PLINK(v1.9b) binary executable file can't be found."
            "Please check 'dependency/' folder or Install PLINK in your system manually."
        )

    ### Class methods
    # @classmethod
    # def setPLINK(cls, _p_plink):
    #     cls.plink_exec = _p_plink


    def __init__(self, _out_prefix, _file_GT=-1, _file_PT=-1, _file_CV=-1, _pheno_name=-1, _covar_name=-1,
                 _condition=-1, _condition_list=-1):

        ### Main Variables ###
        self.GT = -1
        self.PT = -1
        self.CV = -1
        self.pheno_name = -1
        self.covar_name = -1
        self.condition = -1
        self.condition_list = -1

        self.assoc = -1
        self.assoc_sort = -1

        self.f_isPheBinary = -1  # True: Case-Control / False : Quantitative phenotype.

        ### Main Actions ###
        self.setGenotype(_file_GT)

        self.setPhenotype(_file_PT)
        self.setPhenoName(_pheno_name)

        self.setCovariate(_file_CV)
        self.setCovarName(_covar_name)

        self.setCondition(_condition, _condition_list, _out_prefix+'.condlist')


    def setGenotype(self, _file_GT):
        if _file_GT == -1:
            raise HATK_Error.HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME + "No Genotype file was given."
            )
        else:
            self.GT = PLINK.Genotype(_file_GT)


    def setPhenotype(self, _file_PT):
        """
        - is the Phenotype file given? ('*.phe')
        - Does the fam file has phenotype information, too? ('*.fam')
        """
        if _file_PT == -1:
            raise HATK_Error.HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME + "No Phenotype file was given."
            )

            # if self.GT.f_hasPheInfo:
            #     # (1) Use the one in '*.fam'
            #
            #     # (Not realized yet) => Generate a new phe file.
            #     pass
            # else:
            #     # (2) Error.
            #     raise HATK_Error.HATK_InputPreparation_Error(
            #         std_ERROR_MAIN_PROCESS_NAME +
            #         "Phenotype information can't be found in the given fam file('{}')." \
            #         .format(self.GT.fam)
            #     )

        else:
            ## set `self.PT`
            self.PT = PLINK.Phenotype(_file_PT)


    def setCovariate(self, _file_CV):
        if _file_CV == -1:
            pass
        else:
            # set `self.CV`
            self.CV = PLINK.Covariate(_file_CV)


    def setPhenoName(self, _pheno_name):
        if self.PT == -1:
            raise HATK_Error.HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME +
                "The phenotype name('{}') was given but No phenotype file was given." \
                .format(_pheno_name)
            )

        else:
            ## set `self.pheno_name`
            if _pheno_name in self.PT.l_phenotypes:
                self.pheno_name = _pheno_name
            else:
                raise HATK_Error.HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "The given phenotype name('{}') is not in the header of phenotype file('{}')" \
                    .format(_pheno_name, self.PT.file_phe)
                )

            ## set `self.f_isPheBinary`
            sr_phe_target = pd.read_csv(self.PT.file_phe, sep='\s+', header=0, usecols=[self.pheno_name], squeeze=True)

            # All NA values.
            isNA = np.vectorize(lambda x: (x == -9) or (x == 0) or (x == -1))
            arr_isNA = isNA(sr_phe_target.values)
            if np.all(arr_isNA):
                raise HATK_Error.HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "The given phenotype('{}') values are all NAs(-9, 0, -1)."
                )

            # Mixed ints and floats.
            isInt = np.vectorize(lambda x: isinstance(x, int))
            isFloat = np.vectorize(lambda x: isinstance(x, float))

            arr_isInt = isInt(sr_phe_target.values)
            arr_isFloat = isFloat(sr_phe_target.values)

            prop_Int = float(sum(arr_isInt)) / len(arr_isInt)
            prop_Float = float(sum(arr_isFloat)) / len(arr_isFloat)

            if prop_Int > 0.8 and prop_Float < 0.5:
                self.f_isPheBinary = True
            elif prop_Int < 0.5 and prop_Float > 0.8:
                self.f_isPheBinary = False
            else:
                # print("prop_Int: {}(%)".format(prop_Int * 100))
                # print("prop_Float: {}(%)".format(prop_Float * 100))
                raise HATK_Error.HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "HATK can't decide whether the target phenotype('{}') is binary(Case-Control) or continuous(Quantitative).\n"
                    "Please check that phenotype column in the given phenotype file('{}')." \
                    .format(self.pheno_name, self.PT.file_phe)
                )


    def setCovarName(self, _covar_name):
        if self.CV == -1:
            if _covar_name != -1: # Covariate name given.
                print(std_WARNING_MAIN_PROCESS_NAME +
                      "The given covariate names('{}') will be ignored, because No covariate file was given." \
                        .format(_covar_name.split(',')))

        else:
            if _covar_name == -1:
                raise HATK_Error.HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "No any Covariate name was given."
                )
            else:

                f_alliswell = np.all(np.isin(_covar_name.split(','), self.CV.l_covariates))

                if f_alliswell:
                    # set `self.covar_name
                    self.covar_name = _covar_name.split(',')
                else:
                    raise HATK_Error.HATK_InputPreparation_Error(
                        std_ERROR_MAIN_PROCESS_NAME +
                        "Some of given covariate names('{}') are not in the header('{}') of the given covariate file('{}')." \
                        .format(_covar_name.split(','), self.CV.l_covariates, self.CV.file_covar)
                    )


    def setCondition(self, _condition, _condition_list, _out):
        if _condition == -1 and _condition_list == -1: # No condition given.
            pass # Do Nothing.

        elif (_condition != -1) and _condition_list == -1: # condition given.
            with open(_out, 'w') as f_cond:
                for cond in _condition.split(','):
                    f_cond.write(cond+'\n')

                self.condition = _condition
                self.condition_list = _out # as condition list file anyway.

        elif _condition == -1 and (_condition_list != -1): # condition list given.

            if exists(_condition_list):
                with open(_condition_list, 'r') as f_cond:
                    l_conds = [line.rstrip('\n') for line in f_cond]

                    self.condition = ','.join(l_conds)
                    self.condition_list = _condition_list
            else:
                raise HATK_Error.HATK_InputPreparation_Error(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "The given condition list file('{}') can't be found." \
                    .format(_condition_list)
                )

        else:
            raise HATK_Error.HATK_InputPreparation_Error(
                std_ERROR_MAIN_PROCESS_NAME +
                "The '--condition' and '--condition-list' arguments can't be used at the same time."
            )


    def doAssociationTest(self, _out_prefix, _from_mb=29, _to_mb=34, _hide_covar=True,
                          _ci=0.95, _chr=6, _allow_no_sex=True):

        command = [self.plink_exec, "--logistic" if self.f_isPheBinary else "--linear"]

        if _hide_covar: command.append("hide-covar")
        if _allow_no_sex: command.append("--allow-no-sex")

        command.append("--keep-allele-order")
        command.append("--bfile " + self.GT.file_gt)
        command.append("--pheno " + self.PT.file_phe)
        command.append("--pheno-name " + self.pheno_name)
        if self.CV != -1:
            command.append("--covar " + self.CV.file_covar)
            command.append("--covar-name " + ','.join(self.covar_name))
        if self.condition_list != -1:
            command.append("--condition-list {}".format(self.condition_list))
        command.append("--chr {} --from-mb {} --to-mb {}".format(_chr, _from_mb, _to_mb))
        command.append("--ci {}".format(_ci))
        command.append("--out {}".format(_out_prefix))

        self.assoc = Bash_RUN_PLINK(' '.join(command), _out_prefix) + '.assoc.logistic'

        try:
            # Sort the result
            pd.read_csv(self.assoc, sep='\s+', header=0, na_values='NA') \
                .sort_values('P') \
                .fillna('NA') \
                .to_csv(self.assoc+'.sort', sep='\t', header=True, index=False)
        except FileNotFoundError:
            raise HATK_Error.HATK_PLINK_Execution_Error(
                std_ERROR_MAIN_PROCESS_NAME +
                "The association test result('{}') can't be found." \
                .format(self.assoc)
            )
        else:
            self.assoc_sort = self.assoc+'.sort'

        return self.assoc, self.assoc_sort


    def plotManhattan(self, _out_prefix):

        # Temporarily interfaced with this way.
        m_plot = Manhattan([self.assoc], self.assoc, '18',
                           _p_src='../HLA_Manhattan/src', _p_data='../HLA_Manhattan/data')
        return m_plot


    # def doStepwiseConditional(self):
    #     # 'doAssociationTest()' and 'plotManhattan()' will be used, here.
    #     pass


    def __repr__(self):
        str_plink_exec = \
            "[PLINK Binary Executable file]: {}\n".format(self.plink_exec)
        str_GT = \
            "=====< GENOTYPE >=====\n{}\n".format(self.GT)
        str_PT = \
            "=====< PHENOTYPE >=====\n{}\n".format(self.PT)
        str_pheno_name = \
            "[Target Phenotype Name]: '{}'\n".format(self.pheno_name)
        str_isPheBinary = \
            "[Target Phenotype DataType]: {}\n" \
                .format("Binary(Case-Control)" if self.f_isPheBinary else "Continuous(Quantitative)")
        str_CV = \
            "=====< COVARIATE >=====\n{}\n".format(self.CV)
        str_covar_name = \
            "[Target Covariate Name(s)]: {}\n".format(self.covar_name)
        str_condition = \
            "=====< CONDITION >=====\n" \
            "[Condition(s)]: {}\n" \
            "[Condition List File]: {}\n".format(self.condition, self.condition_list)
        str_assoc = \
            "=====< ASSOCIATION TEST >=====\n" \
            "[Association Test({}) Result]: {}\n" \
            "[Association Test({}) Result(Sorted)]: {}\n" \
                .format("logistic" if self.f_isPheBinary else "linear", self.assoc,
                        "logistic" if self.f_isPheBinary else "linear", self.assoc_sort)

        str_summary = ''.join([str_plink_exec, str_GT,
                               str_PT, str_pheno_name, str_isPheBinary,
                               str_CV, str_covar_name,
                               str_condition,
                               str_assoc]).rstrip('\n')
        return str_summary


    def __bool__(self): return self.assoc != -1



def Bash_RUN_PLINK(_command, _out_prefix):
    try:
        # print(_command)
        # print(_command.split())
        subprocess.run(_command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    except subprocess.CalledProcessError:
        raise HATK_Error.HATK_PLINK_Execution_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "PLINK execution for the command('{}') failed. "
            "Please refer to its log('{}') file.".format(_command, _out_prefix + '.log')
        )
    else:
        if exists(_out_prefix+'.nosex'): os.remove(_out_prefix+'.nosex')
        return _out_prefix





if __name__ == '__main__':

    # st = Study("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_20220127/test.All_CD",
    #            _file_GT='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged',
    #            _file_PT='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe',
    #            _file_CV='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar',
    #            _pheno_name="All_CD",
    #            _covar_name='GWAS,Immunochip2')
    #
    # print("===== Association Test =====")
    # assoc = st.doAssociationTest("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_20220127/test.All_CD", )
    # print(assoc)

    # with condition
    st = Study("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_20220127/test.All_CD.cond",
               _file_GT='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged',
               _file_PT='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.phe',
               _file_CV='/media/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.covar',
               _pheno_name="All_CD",
               _covar_name='GWAS,Immunochip2',
               _condition='rs2894066,rs3763338')

    print("===== Association Test =====")
    assoc = st.doAssociationTest("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_20220127/test.All_CD.cond")
    print(assoc)

    print("===== Manhattan Plot =====")
    m_plot = st.plotManhattan('asdf')
    print(m_plot)

    print(st)

    pass