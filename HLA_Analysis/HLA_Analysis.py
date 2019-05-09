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



class HATK_LogisticRegression():

    def __init__(self, _bfile, _out, *args, **kwargs):

        """
        (Checklist)
        1. _pheno or _fam (Phenotype Info)
        2. _covar_name : whether it is single value or multiple value(as csv)
        3. _condition : whether it is single value or multiple value(as csv)
        3.5. _condition_list : it must be a file.(don't need to its inner contents.)
        4. _ref_allele : whether it is given(exists) or not.
        """

        if not kwargs['_phe']:

            kwargs['_phe_name'] = None  # Remove unnecessarily given phenotype column name just in case.

            PheInFam = pd.read_csv(_bfile+'.fam', sep='\s+', header=None, usecols=[5]) # 6th column of *.ped file is phenotype information.

            hasPhenotypeInfo = ~((PheInFam == -9).all())    # If not all of them are -9, then it has phenotype information.
            hasPhenotypeInfo = hasPhenotypeInfo.iat[0]

            if not hasPhenotypeInfo:
                print(std_ERROR_MAIN_PROCESS_NAME + "To perform association test, Phenotype information must be given.\n"
                                                    "Please check '--pheno' and '--pheno-name' arguments again, or\n"
                                                    "Please check the 6th column of '{}.fam' again.".format(_bfile))
                sys.exit()

        else:

            if not os.path.exists(kwargs['_phe']):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given phenotype file('{}') doesn't exist.\n"
                                                    "Please check '--pheno' argument again.".format(kwargs['_phe']))
                sys.exit()


            # Check whether given pheno-name is in phenotype file.
            Phe = open(kwargs['_phe'])

            header =  Phe.readline().rstrip('\n')
            header_items = re.split(r'\s+', header)

            if len(header_items) < 3:
                # Less than 3 columns
                print(std_ERROR_MAIN_PROCESS_NAME + "Given phenotype file must have more than 2 columns.\n"
                                                    "Please check '--pheno' argument again.")
                sys.exit()

            elif len(header_items) == 3:

                # Setting only one column as phenotype information.
                kwargs['_phe_name'] = header_items[-1]

            else:
                # More than 3 columns
                if not kwargs['_phe_name'] in header_items:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Given phenotype name('{}') can't be found in the columns of the given phenotype file.\n"
                                                        "Please check '--pheno' and '--pheno-name' arguments again.".format(kwargs['_phe_name']))
                    sys.exit()

            Phe.close()


        if not kwargs['_covar']:

            # Covariate file not given.

            kwargs['_covar_name'] = None  # Remove unnecessarily given covariate column name(s) just in case.


        else:

            # Covariate file given.

            if not os.path.exists(kwargs['_covar']):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given covariate file('{}') doesn't exist.\n"
                                                    "Please check '--covar' argument again.".format(kwargs['_covar']))
                sys.exit()


            # Checking covariate columns names

            # p_sep = re.compile(r',\s*')
            #
            # header_covar = re.split(pattern=r'\s+', string=open(kwargs['_covar'], 'r').readline().rstrip('\n'))
            # items_covar = p_sep.split(string=kwargs['_covar_name'])



        if kwargs['_condition_list']:

            if not os.path.isfile(kwargs['_condition_list']):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given condition list('{}') is not a file.\n"
                                                    "Please check '--condition-list' argument again.".format(kwargs['_condition_list']))
                sys.exit()

            """
            (2019.05.09.)
            In Plink1.9, It seems '--condition' argument takes only one variant ID. i.e. multiple variants ID can't be
            passed as csv. 
            
            For example, 
            => "--condition AA_A_-15_29910359_exon1,AA_A_-11_29910371_exon1,AA_A_9_29910558_exon2_F"
            
            This kind of usage is no longer possible.
            
            If user wants to pass multiple variants as covariate, there is no choice but to give tsv file of those IDs 
            with '--condition-list'.
              
            """



        if not kwargs['_ref_allele']:

            kwargs['_ref_allele'] = MakeDefaultReferenceAllele(_bfile, _out)


        ################################################################################################################


        self.results = Logistic_Regression(_bfile, _out, _phe=kwargs['_phe'], _phe_name=kwargs['_phe_name'],
                                           _covar=kwargs['_covar'], _covar_name=kwargs['_covar_name'],
                                           _condition=kwargs['_condition'], _condition_list=kwargs['_condition_list'],
                                           _ref_allele=kwargs['_ref_allele'])


    def getResults(self):
        return self.results





class HATK_MetaAnalysis():

    def __init__(self, _out, _assoc_result):


        if not isinstance(_assoc_result, list):
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--assoc-result' argument again.")
            sys.exit()

        self.results = Meta_Analysis(_out, _assoc_result)


    def getResults(self):
        return self.results





########## < Association Test > ##########


### (Logistic Regression)
def Logistic_Regression(_bfile, _out,
                        _phe=None, _phe_name=None,
                        _covar=None, _covar_name=None,
                        _condition=None, _condition_list=None,
                        _ref_allele=None, _from_mb=29, _to_mb=34, _hide_covar=True, _ci=0.95,
                        _chr=6, _allow_no_sex=True):


    ### Argument Processing && Generating Command.

    command = [GLOBAL_p_plink, "--noweb", "--logistic", "--out {0}".format(_out)]

    if bool(_bfile):
        command.append("--bfile {0}".format(_bfile))

    if bool(_covar):
        command.append("--covar " + _covar)
    if bool(_covar_name):
        command.append("--covar-name " + _covar_name)
    if bool(_phe):
        command.append("--pheno " + _phe)
    if bool(_phe_name):
        command.append("--pheno-name " + _phe_name)

    if bool(_condition):
        command.append("--condition " + _condition)
    elif bool(_condition_list):
        command.append("--condition-list " + _condition_list)

    if bool(_ref_allele):
        # command.append("--reference-allele " + _ref_allele) # Plink v1.07
        command.append("--a1-allele " + _ref_allele) # Plink v1.9

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


    __RETURN__ = _out + ".assoc.logistic"

    if os.path.exists(__RETURN__):
        return __RETURN__
    else:
        return -1


def MakeDefaultReferenceAllele(_bfile, _out):

    # step1. Load ".bim" file
    bim = pd.read_csv(_bfile+".bim", sep='\t', header=None, usecols=[1,4,5], names=["Label", "Al1", "Al2"])

    l_Al1 = bim.iloc[:, 1].tolist()
    l_Al2 = bim.iloc[:, 2].tolist()

    for i in range(0, bim.shape[0]):

        if (l_Al1[i] == "a" and l_Al2[i] == "p"):
            l_Al1[i] = "p"

    # Making the reference allele DataFrame.

    _out2 = os.path.join(os.path.dirname(_out), os.path.basename(_bfile)+'.refallele')

    df_Ref_Allele = pd.concat([bim.iloc[:, 0], pd.Series(l_Al1)], axis=1)
    df_Ref_Allele.to_csv(_out2, sep='\t', header=False, index=False)

    return _out2



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
def Meta_Analysis(_out, _assoc_result):


    ### Generating Command

    command = [GLOBAL_p_plink, "--noweb", "--out {0}".format(_out), "--meta-analysis", *_assoc_result]
    command = ' '.join(command)


    # Conducting Meta-Analysis by Plink(v1.9)

    print(command)

    if not os.system(command):

        if os.path.exists(_out + ".meta"):
            return (_out + ".meta")
        else:
            return -1

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Meta-Analysis failed.")
        return -1




