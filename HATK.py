#-*- coding: utf-8 -*-

import os, sys, re
import logging
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

    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="HG")

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
    # parser.add_argument("--results-assoc", "-ra", help="\nResult files conducted by Logistic Regression Tests(ex.\"*.assoc.logistic\").\n\n",
    #                     nargs='+')
    parser.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(s)(ex. *.assoc.logistic).\n\n", nargs='+')


    ##### < IMGT2Sequence > #####

    g_IMGT2Sequence = parser.add_argument_group(title='IMGT2Seq', description='')

    g_IMGT2Sequence.add_argument("--imgt2seq", help="\nGive this argument to implement \"IMGT2Seq\" sub-module.\n\n", action='store_true')

    g_IMGT2Sequence.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n")

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

    # Additional input ped file type.
    PED_TYPE = g_NomenCleaner.add_mutually_exclusive_group()
    # PED_TYPE.add_argument("-ped", help="\nHLA Type Data.\n\n") # moved to main group.
    PED_TYPE.add_argument("--hped_Ggroup", help="\nHLA Type Data(G-group allele \"*.hped\" file).\n\n", dest="hped_G")
    PED_TYPE.add_argument("--hped_Pgroup", help="\nHLA Type Data(P-group allele \"*.hped\" file).\n\n", dest="hped_P")

    g_NomenCleaner.add_argument("-iat", help="\nIntegrated Allele Table file(*.iat).\n\n")

    # Output format selection
    format_selection = g_NomenCleaner.add_mutually_exclusive_group()


    format_selection.add_argument("--1field", help="\nOutput *.chped file as '1-field' format.\n\n",
                                        action="store_true", dest="oneF")
    format_selection.add_argument("--2field", help="\nOutput *.chped file as '2-field' format.\n\n",
                                        action="store_true", dest="twoF")
    format_selection.add_argument("--3field", help="\nOutput *.chped file as '3-field' format.\n\n",
                                        action="store_true", dest="threeF")
    format_selection.add_argument("--4field",
                                        help="\nOutput chped file as '4-field(Current Standard Names)' format.\n\n",
                                        action="store_true", dest="fourF")
    format_selection.add_argument("--G-group", help="\nOutput *.chped file as 'G-group' format.\n\n",
                                        action="store_true")
    format_selection.add_argument("--P-group", help="\nOutput *.chped file as 'P-group' format.\n\n",
                                        action="store_true")

    format_selection.add_argument("--old-format", help="\nOutput hped file which is compatible with original MakeReference. (ex. A:01:01)\n\n",
                                        action="store_true")

    # Additional utility flags
    g_NomenCleaner.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*01:01\" -> \"01:01\").\n\n", action='store_true')
    g_NomenCleaner.add_argument("--leave-NotFound", help="\nLeaving HLA alleles which can't be found in given *.iat file(Novel or Erroneous allele) intact.\n\n", action='store_true')



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



    ### Association Test 1 (Logistic Regression by Plink 1.07)

    g_ASSOC1_Logistic = parser.add_argument_group(title='(Association Test 1) Logistic Regression', description='')

    g_ASSOC1_Logistic.add_argument("--logistic", "-lr", help="\nGive this argument to implement \"Logistic Regression\". (Plink v1.07)\n\n", action='store_true')



    ### Association Test 2 (Omnibus Test by Buhm Han.)

    g_ASSOC2_Omnibus = parser.add_argument_group(title='(Association Test 2) Omnibus Test', description='')
    g_ASSOC2_Omnibus.add_argument("--omnibus", "-om", help="\nGive this argument to implement \"Omnibus Test\".\n\n", action='store_true')

    g_ASSOC2_Omnibus.add_argument("--phased", "-ph", help="\nSpecify the \"*.bgl.phased\" file.\n\n")
    g_ASSOC2_Omnibus.add_argument("--aa", help="\nSpecify the \"*.aa\" file.\n\n")
    # g_ASSOC2_Omnibus.add_argument("--rare-threshold", "-rth", help="\nSpecify the threshold value for rare variants.\n\n")



    ### Meta-Analysis (by Plink 1.07)

    g_MetaAnalysis = parser.add_argument_group(title='MetaAnalysis', description='')

    g_MetaAnalysis.add_argument("--meta-analysis", help="\nGive this argument to implement \"Meta-Analysis\" sub-module. (Plink v1.07)\n\n", action='store_true')

    # g_MetaAnalysis.add_argument("--rassoc", "-ra", help="\nSpecify the Result file(s) of association test for Meta-Analysis\n"
    #                                                     "(ex. \"*.asssoc.logistic\", etc.)\n\n", nargs='+')



    ##### < Plotting related. > #####


    ### heatmap

    g_heatmap = parser.add_argument_group(title='(Plotting 1) Heatmap', description='')
    g_heatmap.add_argument("--heatmap", help="\nGenerate Heatmap Plot.\n\n", action="store_true")

    g_heatmap.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n", choices = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'])
    g_heatmap.add_argument("--maptable", "-mt", help="\nMarker Dictionary file(Maptable) generated by 'IMGTt2Sequence'.\n\n")
    # g_heatmap.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(ex. *.assoc.logistic).\n\n", required=True)
    g_heatmap.add_argument("--as4field", help="\nShow HLA allele names in 4-field format\n\n", action='store_true')



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

    # args = parser.parse_args(["--imgt2sequence", "-imgt", "370", "-o", "TEST/TEST", "-hg", "18"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)


    from src.HLA_Study import HLA_Study
    myStudy = HLA_Study(args)


    # ########## < Main > ##########
    #
    # if not flag_MODULEs:
    #
    #     ##### < Whole Implementation > #####
    #
    #     # ### [1] IMGT2Sequence
    #     #
    #     # previous_dictionary = IMGT_Dictionary(hg=args.hg, imgt=args.imgt, out=args.out)
    #     #
    #     # if previous_dictionary:
    #     #
    #     #     log_HATK.info(previous_dictionary)
    #     #
    #     #     __dict_AA__, __dict_SNPS__, __IAT__, __d_MapTable__ = previous_dictionary.getDictionaries()
    #     #
    #     # else:
    #     #
    #     #
    #     #     __dict_AA__, __dict_SNPS__, __IAT__, __d_MapTable__ = HATK_IMGT2Sequence(args.hg, args.out, args.imgt,
    #     #                                                                              _imgt_dir=args.imgt_dir,
    #     #                                                                              _no_MultiP=args.no_multiprocess)
    #     #
    #     #     t_string = \
    #     #         "< IMGT2Sequence >\n" \
    #     #         "- HLA Amino Acids : {}\n" \
    #     #         "- HLA SNPs : {}\n" \
    #     #         "- Integrated Allele Table : {}\n" \
    #     #         "- Maptables for heatmap : \n" \
    #     #         "   A   : {A}\n" \
    #     #         "   B   : {B}\n" \
    #     #         "   C   : {C}\n" \
    #     #         "   DPA1: {DPA1}\n" \
    #     #         "   DPB1: {DPB1}\n" \
    #     #         "   DQA1: {DQA1}\n" \
    #     #         "   DQB1: {DQB1}\n" \
    #     #         "   DRB1: {DRB1}\n".format(__dict_AA__, __dict_SNPS__, __IAT__, **__d_MapTable__)
    #     #
    #     #     log_HATK.info(t_string)
    #     #
    #     #
    #     #
    #     # ### [2] Summary of Study Information
    #     # myHLAstudy = HLA_Analysis.Study(rhped=args.rhped, hped=args.hped, chped=args.chped, platform=args.platform,
    #     #                                 variants=args.variants, pheno=args.pheno, pheno_name=args.pheno_name,
    #     #                                 covar=args.covar, covar_name=args.covar_name, condition=args.condition, condition_list=args.condition_list,
    #     #                                 refallele=args.reference_allele, input=args.input, hg=args.hg, out=args.out,
    #     #                                 dict_AA=__dict_AA__, dict_SNPS=__dict_SNPS__, iat=__IAT__, d_MapTable=__d_MapTable__)
    #
    #     """
    #     (2019. 1. 14.)
    #     Argument 마지막 부분
    #     "   dict_AA=__dict_AA__, dict_SNPS=__dict_SNPS__, iat=__IAT__, d_MapTable=__d_MapTable__)"
    #
    #     이 부분을 나중에 그냥 HLA_Dictionary instance하나 만들어서 이 instance를 받는 식으로 바꿀것.
    #
    #     """
    #
    #     # log_HATK.info(myHLAstudy)
    #
    #
    #     pass
    #
    #
    #
    #
    # else:
    #     ##### < Single Implementation > #####
    #
    #     if args.imgt2sequence:
    #
    #         hla_dictionaries = IMGT_Dictionary(args.imgt, args.hg, args.out,
    #                                            _no_indel=args.no_indel, _multiprocess=args.multiprocess,
    #                                            _save_intermediates=args.save_intermediates, _imgt_dir=args.imgt_dir)
    #         # print(hla_dictionaries)
    #
    #     elif args.hla2hped:
    #         pass
    #
    #     elif args.nomencleaner:
    #
    #         from src.NomenCleaner.NomenCleaner import NomenCleaner
    #
    #         __chped__ = NomenCleaner(_hped=args.hped, _hped_G=args.hped_G, _hped_P=args.hped_P, _iat=args.iat, _out=args.out,
    #                                  _1field=args.oneF, _2field=args.twoF, _3field=args.threeF, _4field=args.fourF,
    #                                  _Ggroup=args.G_group, _Pgroup=args.P_group, _old_format=args.old_format,
    #                                  __NoCaption=args.NoCaption,
    #                                  _imgt=args.imgt)
    #
    #
    #     elif args.bmarkergenerator :
    #         from src.b_MarkerGenerator.b_MarkerGenerator import b_MarkerGenerator
    #
    #         __b_Marker__ = b_MarkerGenerator(args.chped, args.out, args.hg, args.dict_AA, args.dict_SNPS, _variants=args.variants,
    #                                          __save_intermediates=args.save_intermediates)
    #
    #
    #     elif args.logistic:
    #         from src.HLA_Analysis.HLA_Analysis_modules import Logistic_Regression
    #
    #
    #
    #
    #     elif args.manhattan:
    #         from src.HLA_Analysis.Plotting.manhattanplot.manhattan import HATK_manhattan
    #
    #         """
    #         HATK_manhattan(_results_assoc, _plot_label, _out, _hg, _pointcol="#778899", _topcol="#FF0000", _min_pos="29.60E6", _max_pos="33.2E6","""
    #
    #         mamhattan_plot = HATK_manhattan(args.results_assoc, args.title, args.out, args.hg, _pointcol=args.point_color, _topcol=args.top_color)
    #
    #     else:
    #         # print(std_ERROR_MAIN_PROCESS_NAME + "이상해~")
    #
    #         from src.HLA_Analysis.Plotting.heatmap.heatmap import HATK_HEATMAP
    #
    #         heatmap_plot = HATK_HEATMAP(args.HLA, args.out, args.maptable, args.result_assoc_AA, args.result_assoc_HLA)
