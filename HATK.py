#-*- coding: utf-8 -*-

import os, sys, re
# import logging
import argparse, textwrap


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]




if __name__ == "__main__":


    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        HATK.py (HLA Analysis Toolkit)

        (Created by Wanson Choi, 2018/05/27)


        <SPECIFICATION> 
        
        [1] IMGT2Seq
        [2] b:MarkerGenerator
        [3] NomenCleaner
        [4] Logistic Regression
        [5] Omnibus Test
        [6] Meta-Analysis
        [7] Heatmap
        [8] Manhattan
        [9] HLA2HPED


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )


    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"
    parser._optionals.description = '\n!-- COMMON ARGUMENTS TO BE USED IN ALL MODULES. --!\n\n'

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n") # Core argument of this program.
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n")

    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="HG")

    parser.add_argument("--variants", help="\nNormal SNP marker data.\n\n")


    ## Major 3 ways to take HLA type data.
    HLA_TYPE = parser.add_mutually_exclusive_group()

    # (1) Output results from other HLA software platform. (No specific file extension).
    HLA_TYPE.add_argument("--rhped", help="\nOutput result(s) from other HLA software to be processed by \'HLA2HPED\'.\n\n", nargs='*')
    # (2) "*.hped"
    HLA_TYPE.add_argument("--hped", help="\nHLA Type Data not yet processed by \'NomenCleaner\'.\n\n")
    # (3) "*.chped"
    HLA_TYPE.add_argument("--chped", help="\nHLA Type Data processed by \'NomenCleaner\' (*.chped)\n\n")

    # (Meta-Analysis or Manhattan Plot)
    parser.add_argument("--ar", help="\nAssociation test result file(s) (ex. *.assoc.logistic, *.omnibus).\n\n", nargs='+')


    ##### < IMGT2Sequence > #####

    g_IMGT2Sequence = parser.add_argument_group(title='IMGT2Seq', description='')

    g_IMGT2Sequence.add_argument("--imgt2seq", help="\nGive this argument to implement \"IMGT2Seq\" sub-module.\n\n", action='store_true')

    g_IMGT2Sequence.add_argument("--imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n")

    g_IMGT2Sequence.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    g_IMGT2Sequence.add_argument("--multiprocess", help="\nSetting off parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    g_IMGT2Sequence.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true') # It will be shared with 'b:MarkerGenerator'
    g_IMGT2Sequence.add_argument("--imgt-dir", help="\nGiving direct path to the directory of IMGT data folder.\n\n")



    ##### < b:MarkerGenerator > #####

    g_b_MarkerGenerator = parser.add_argument_group(title='b:MarkerGenerator', description='')

    g_b_MarkerGenerator.add_argument("--bmarkergenerator", help="\nGive this argument to implement \"b:MarkerGenerator\" sub-module.\n\n",
                                     action='store_true')

    g_b_MarkerGenerator.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n")
    g_b_MarkerGenerator.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")



    ##### < NomenCleaner > #####

    g_NomenCleaner = parser.add_argument_group(title='NomenCleaner', description='')

    g_NomenCleaner.add_argument("--nomencleaner", help="\nGive this argument to implement \"NomenCleaner\" sub-module.\n\n",
                                action='store_true')

    g_NomenCleaner.add_argument("--hat", help="\nHLA Allele Table file(*.hat).\n\n")

    # Output format selection
    format_selection = g_NomenCleaner.add_mutually_exclusive_group()
    format_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n", action="store_true", dest="oneF")
    format_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n", action="store_true", dest="twoF")
    format_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n", action="store_true", dest="threeF")
    format_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields.\n\n", action="store_true", dest="fourF")
    format_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n", action="store_true")
    format_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n", action="store_true")


    # Additional utility flags
    g_NomenCleaner.add_argument("--NoCaption", help="\nMake converted HLA alleles NOT have HLA gene prefix(ex. \"A*\").\n\n", action='store_true')
    g_NomenCleaner.add_argument("--leave-NotFound", help="\nLeaving HLA alleles which can't be found in given *.hat file(Novel or Erroneous allele) intact.\n\n", action='store_true')



    ##### < HLA Analysis related > #####

    ### Common arguments for Association Tests

    g_ASSOC = parser.add_argument_group(title='(Association Test) Common arguments for Association Tests', description='')

    g_ASSOC.add_argument("--covar", help="\nSpecify .covar file (Plink v1.9).\n\n")
    g_ASSOC.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use. (Plink v1.9)\n\n")

    g_ASSOC.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.9).\n\n")
    g_ASSOC.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use. (Plink v1.9)\n\n")

    CondVars = g_ASSOC.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify a single variant ID to condition(i.e. To set it as covariate). (Plink v1.9)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify a tsv file of multiple variant IDs to condition(i.e. To set it as covariates). (Plink v1.9)\n\n")

    g_ASSOC.add_argument("--reference-allele", help="\nSpecify the Reference Allele file('--a1-allele'). (Plink v1.9)\n\n")



    ### Association Test 1 (Logistic Regression by Plink 1.9)

    g_ASSOC1_Logistic = parser.add_argument_group(title='(Association Test 1) Logistic Regression', description='')

    g_ASSOC1_Logistic.add_argument("--logistic", "-lr", help="\nGive this argument to implement \"Logistic Regression\". (Plink v1.9)\n\n", action='store_true')



    ### Association Test 2 (Omnibus Test by Buhm Han.)

    g_ASSOC2_Omnibus = parser.add_argument_group(title='(Association Test 2) Omnibus Test', description='')
    g_ASSOC2_Omnibus.add_argument("--omnibus", "-om", help="\nGive this argument to implement \"Omnibus Test\".\n\n", action='store_true')

    g_ASSOC2_Omnibus.add_argument("--fam", help="\nSpecify \"*.fam\" file.\n\n")
    g_ASSOC2_Omnibus.add_argument("--phased", "-ph", help="\nSpecify the phased beagle file(ex. \"*.bgl.phased\") to make \"*.aa\" file.\n\n")
    g_ASSOC2_Omnibus.add_argument("--aa", help="\nSpecify the \"*.aa\" file to be used in Omnibus test.\n\n")
    # g_ASSOC2_Omnibus.add_argument("--rare-threshold", "-rth", help="\nSpecify the threshold value for rare variants.\n\n")



    ### Meta-Analysis

    g_MetaAnalysis = parser.add_argument_group(title='MetaAnalysis', description='')

    g_MetaAnalysis.add_argument("--metaanalysis", help="\nPeforming MetaAnalysis(Inverse-Variance Weight Method(Fixed effect) on the results of HATK.\n\n", action='store_true')

    parser.add_argument("--s1-logistic-result", "-s1lr", help="\nLogistic Regression Result file of Study 1(ex. *.assoc.logistic).\n\n")
    parser.add_argument("--s1-bim", "-s1b", help="\nThe bim file used in Study 1 Logistic Regression(ex. *.bim).\n\n")

    parser.add_argument("--s2-logistic-result", "-s2lr", help="\nLogistic Regression Result file of Study 2(ex. *.assoc.logistic).\n\n")
    parser.add_argument("--s2-bim", "-s2b", help="\nThe bim file used in Study 2 Logistic Regression(ex. *.bim).\n\n")



    ##### < Plotting related. > #####


    ### heatmap

    g_heatmap = parser.add_argument_group(title='(Plotting 1) Heatmap', description='')
    g_heatmap.add_argument("--heatmap", help="\nGenerate Heatmap Plot.\n\n", action="store_true")

    g_heatmap.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n", choices = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')
    g_heatmap.add_argument("--maptable", "-mt", help="\nMarker Dictionary file(Maptable) generated by 'IMGTt2Sequence'.\n\n")



    ### manhattan

    g_manhattan = parser.add_argument_group(title='(Plotting 2) Manhattan', description='')
    g_manhattan.add_argument("--manhattan", help="\nGenerate Manhattan Plot.\n\n", action="store_true")

    # g_manhattan.add_argument("--title", "-pt", help="\nThe title of output plot.\n\n", default="Manhattan")
    g_manhattan.add_argument("--point-color", "-pc", help="\nPoint color(ex. \"#778899\").\n\n", default="#778899")
    g_manhattan.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n", default="#FF0000")

    g_manhattan.add_argument("--point-size", "-ps", help="\nGeneral point size (default: 15).\n\n", default="15")
    g_manhattan.add_argument("--yaxis-unit", "-yau", help="\nY-axis value(-log10(x)) unit (default : 10).\n\n", default="10")


    ##### < HLA2HPED > #####

    g_hla2hped = parser.add_argument_group(title='HLA2HPED', description='')

    g_hla2hped.add_argument("--hla2hped", help="\nGive this argument to implement \"HLA2HPED\" sub-module.\n\n", action='store_true')

    g_hla2hped.add_argument("--platform", "-p", help="\nSoftware platform.(ex. 'HIBAG', 'AXIOM', etc.)\n\n",
                            choices=["AXIOM", "HIBAG", "xHLA"])









    ##### < for Testing > #####

    # args = parser.parse_args('--omnibus --out /Users/wansun/Git_Projects/HATK/tests/_4_HLA_Analysis/Omnibus/merged/20190618_merged_covar --input /Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/data/HLA_Analysis/sample/Merged/merged --covar-name GWAS --pheno-name All_UC'.split(' '))




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)


    from src.HLA_Study import HLA_Study
    myStudy = HLA_Study(args)