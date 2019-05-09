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
    parser.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(ex. *.assoc.logistic).\n\n", nargs='+')


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

    g_ASSOC.add_argument("--covar", help="\nSpecify .covar file (Plink v1.07).\n\n")
    g_ASSOC.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use. (Plink v1.07)\n\n")

    g_ASSOC.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.07).\n\n")
    g_ASSOC.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use. (Plink v1.07)\n\n")

    CondVars = g_ASSOC.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify Marker name(s) to condition as comma-separated(\",\"). (Plink v1.07)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify the Marker name(s) to condition as a file. (Plink v1.07)\n\n")

    g_ASSOC.add_argument("--reference-allele", help="\nSpecify the Reference Allele file. (Plink v1.07)\n\n")



    ### Association Test 1 (Logistic Regression by Plink 1.07)

    g_ASSOC1_Logistic = parser.add_argument_group(title='(Association Test 1) Logistic Regression', description='')

    g_ASSOC1_Logistic.add_argument("--logistic", "-lr", help="\nGive this argument to implement \"Logistic Regression\". (Plink v1.07)\n\n", action='store_true')



    ### Association Test 2 (Omnibus Test by Buhm Han.)

    g_ASSOC2_Omnibus = parser.add_argument_group(title='(Association Test 2) Omnibus Test', description='')
    g_ASSOC2_Omnibus.add_argument("--omnibus", "-om", help="\nGive this argument to implement \"Omnibus Test\".\n\n", action='store_true')

    g_ASSOC2_Omnibus.add_argument("--phased", "-ph", help="\nSpecify the \".aa\" file.\n\n")
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


    # ### Preprocessing arguments
    #
    # # "--input" arguments
    #
    # if args.input:
    #     #     """
    #     #     * --input := the common prefix to point out next files
    #     #         (1) *.hped or *.chped
    #     #         (2) normal_SNPs
    #     #         (3) phenotype file(*.phe)
    #     #         (4) covariate file(*.covar)
    #     #         (5) condition (*.condvars)
    #     #         (6) reference allele (*.refallele)
    #     #     """
    #     #     if not args.hped:
    #     #         args.hped = args.input + ".hped" if os.path.exists(args.input + ".hped") else None
    #     #     if not args.chped:
    #     #         args.chped = args.input + ".chped" if os.path.exists(args.input + ".chped") else None
    #     #
    #     #     # `normal_SNPs` => just as prefix.
    #     #     if not args.variants:
    #     #         args.variants = args.input if os.path.exists(args.input + ".bed") and os.path.exists(args.input + ".bim") and os.path.exists(args.input + ".fam") else None
    #     #     if not args.pheno:
    #     #         args.pheno = args.input + ".phe" if os.path.exists(args.input + ".phe") else args.input + ".pheno" if os.path.exists(args.input + ".pheno") else None
    #     #     if not args.covar:
    #     #         args.covar = args.input + ".covar" if os.path.exists(args.input + ".covar") else args.input + ".cov" if os.path.exists(args.input + ".cov") else None
    #     #
    #     #     if not args.condition_list:
    #     #         args.condition_list = args.input + ".condvars" if os.path.isfile(args.input + ".condvars") else None
    #     #     if not args.refallele:
    #     #         args.refallele = args.input + ".refallele" if os.path.isfile(args.input + ".refallele") else args.input + ".ref" if os.path.exists(args.input + ".ref") else None
    #
    #
    # # Which module?
    # flag_MODULEs = (args.imgt2sequence or args.bmarkergenerator or args.nomencleaner or
    #                 args.logistic or args.omnibus or args.meta_analysis or
    #                 args.heatmap or args.manhattan or args.hla2hped)
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # ### Setting Output path
    # if not args.out:
    #     print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out"))
    #     sys.exit()
    # else:
    #     args.out = args.out if not args.out.endswith('/') else args.out.rstrip('/')
    #     if bool(os.path.dirname(args.out)): os.makedirs(os.path.dirname(args.out), exist_ok=True)
    #
    #
    # ### Setting Logger
    #
    # file_handler1 = logging.FileHandler(args.out+'.hatk.log', mode='w')
    # file_handler1.setFormatter(logging.Formatter((std_MAIN_PROCESS_NAME+'%(message)s')))
    #
    # log_HATK = logging.getLogger("HATK")
    # log_HATK.setLevel(logging.INFO)
    # log_HATK.addHandler(file_handler1)
    # log_HATK.info(args)
    #
    #
    #
    #
    #
    #
    #
    #
    # from src.IMGT2Sequence.IMGT2Seqeunce import IMGT_Dictionary
    # import src.HLA_Analysis.HLA_Analysis as HLA_Analysis
    #
    #
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





    #     ### Single Module implementation.
    #
    #
    #     if args.imgt2sequence:
    #
    #         ##### IMGT2Sequence #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing IMGT2Sequence.")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. -hg
    #         2. -o (*)
    #         3. -imgt
    #         """
    #
    #         __dict_AA__, __dict_SNPS__, __IAT__, __d_MapTable__ \
    #             = HATK_IMGT2Sequence(args.hg, args.out, args.imgt,
    #                                  _no_Indel=args.no_indel, _no_MultiP=args.no_multiprocess,
    #                                  _save_intermediates=args.save_intermediates, _imgt_dir=args.imgt_dir)
    #
    #
    #
    #
    #     if args.nomencleaner:
    #
    #         ##### NomenCleaner #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing NomenCleaner.")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. either -ped, -ped-Ggroup or -ped-Pgroup
    #         2. -iat
    #         3. -o (*)
    #         4. output format(ex. --1field, etc.)
    #
    #         (optionals)
    #         5. --No-caption
    #         """
    #
    #         t_hped = ""
    #         t_hped_descriptor = -1
    #         t_FIELD_FORMAT = -1
    #
    #         # hped type check
    #         if not(bool(args.hped) or bool(args.hped_G) or bool(args.hped_P)):
    #             print(std_ERROR_MAIN_PROCESS_NAME + "You must give at least one argument among \"-hped\", \"-hped-Ggroup\", \"-hped-Pgroup\".\n")
    #             sys.exit()
    #         else:
    #
    #             if bool(args.hped):
    #                 t_hped_descriptor = 1
    #                 t_hped = args.hped
    #             elif bool(args.hped_G):
    #                 t_hped_descriptor = 2
    #                 t_hped = args.hped_G
    #             elif bool(args.hped_P):
    #                 t_hped_descriptor = 3
    #                 t_hped = args.hped_P
    #             else:
    #                 print(std_ERROR_MAIN_PROCESS_NAME + "The arguments related to ped(\"-hped\", \"-hped-Ggroup\" or \"-hped-Pgroup\") has wrong values."
    #                                                     "Please check them again.\n")
    #                 sys.exit()
    #
    #         # field format request check
    #         if not(bool(args.oneF) or bool(args.twoF) or bool(args.threeF) or bool(args.fourF) or bool(args.G_group) or bool(args.P_group)):
    #             print(std_ERROR_MAIN_PROCESS_NAME + "The argument related to field format has not given. Please check it again.\n")
    #             sys.exit()
    #         else:
    #             t_FIELD_FORMAT = (1 if bool(args.oneF) else 2 if bool(args.twoF) else 3 if bool(args.threeF) else 4
    #                                 if bool(args.fourF) else 5 if bool(args.G_group) else 6 if bool(args.P_group) else -1)
    #
    #
    #
    #         HATK_NomenCleaner(_p_hped=t_hped, _ped_descriptor=t_hped_descriptor, _p_iat=args.iat, _out=args.out,
    #                           _field_format=t_FIELD_FORMAT, _f_NoCaption=bool(args.NoCaption))
    #
    #
    #
    #
    #     if args.bmarkergenerator:
    #
    #         ##### b:MarkerGenerator #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing b:MarkerGenerator.")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. -input (*)
    #         2. -hped or -ped
    #         3. -hg
    #         4. -o (*)
    #         5. --dict-AA
    #         6. --dict-SNPS
    #         """
    #
    #         HATK_b_MarkerGenerator(_CHPED=args.chped, _OUT=args.out, _hg=args.hg,
    #                                _dictionary_AA=args.dict_AA, _dictionary_SNPS=args.dict_SNPS,
    #                                _plain_SNP_DATA=args.input)
    #
    #
    #
    #
    #     if args.logistic:
    #
    #         ##### (Association Test 1) Logistic regression(by Plink 1.07) #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing Logistic Regression(Association Test, plink 1.07).")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. -input(--bfile) (*)
    #         2. -o (*)
    #
    #         (optionals)
    #         3. -covar (with -covar-name)
    #         4. -phe (with -phe-name)
    #         5. either --condition or --condition-list
    #         6. --reference-allele
    #
    #         """
    #
    #         t_result_logistic = HLA_Analysis.HATK_ASSOC1_Logistic_Regression(
    #             args.input, args.out, args.covar, args.covar_name, args.pheno, args.pheno_name,
    #             args.condition, args.condition_list, args.reference_allele)
    #
    #
    #
    #     if args.omnibus:
    #
    #         ##### (Association Test 2) Omnibus Test #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing Omnibus Test(Association Test).")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. --input
    #         2. -o (*)
    #         3. --phased
    #         4. --pheno
    #         5. --pheno-name
    #         6. --covar
    #         7. --covar-names
    #         8. --rare-threshold (2018. 10. 04.) Postponed.
    #         9. --condition
    #         10. --condition-list
    #
    #
    #         phased파일이 주어질 수 있는 방식이 크게 두가지인거지
    #
    #         1. output prefix공유해서 같이 input으로 같이 활용되든지
    #         2. 아니면 --phased로 따로 주든지.
    #
    #         그리고 다시 한번 정리하자면 phased옵션과 관련해서 가장 중요한건 *.fam파일과 .bgl.phased파일임.
    #         """
    #
    #         HLA_Analysis.HATK_ASSOC2_Omnibus_Test(_input=args.input, _out=args.out, _phased=args.phased,
    #                                               _phe=args.pheno, _phe_name=args.pheno_name,
    #                                               _covar=args.covar, _covar_names=args.covar_name,
    #                                               _condition=args.condition, _condition_list=args.condition_list)
    #
    #
    #     if args.meta_analysis:
    #
    #         ##### Meta-Analysis #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing Meta-Analysis(Plink1.07).\n")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. -o (*)
    #         2. -ra
    #
    #         """
    #
    #         __meta__ = HLA_Analysis.HATK_META_ANALYSIS(args.out, args.results_assoc)
    #
    #
    #
    #     if args.manhattan:
    #
    #         ##### manhattan #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing Manhattan(Plotting).")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. --result-assoc
    #         2. --out (*)
    #         3. -hg
    #
    #         (optionals)
    #         4. --point-color
    #         5. --top-color
    #         """
    #
    #         __manhattan__ = HATK_manhattan(args.results_assoc, args.title, args.out, args.hg,
    #                                        _pointcol=args.point_color, _topcol=args.top_color)
    #
    #
    #
    #
    #     if args.heatmap:
    #
    #         ##### heatmap #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing Heatmap(Plotting).")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. --hla
    #         2. -o
    #         3. --maptable
    #         4. -lrA (AA logistic regression result)
    #         5. -lrH (HLA logistic regression result)
    #         """
    #
    #         __heatmap__ = HATK_HEATMAP(args.HLA, args.out, args.maptable,
    #                                    args.result_assoc_AA, args.result_assoc_HLA)
    #
    #
    #
    #
    #     if args.hla2hped:
    #
    #         ##### Converter #####
    #
    #         print(std_MAIN_PROCESS_NAME + "Implementing HLA2HPED.")
    #
    #         """
    #         List of necessary arguments.
    #
    #         1. --platform
    #         2. --rhped
    #         3. --out(*)
    #
    #         No optional arguments.
    #         """
    #
    #         __HLA_type__ = HATK_HLA2HPED(args.rhped, args.out, args.platform)

