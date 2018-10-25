#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd


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

    ## Major 3 ways to take HLA type data.
    HLA_TYPE = parser.add_mutually_exclusive_group()

    # (1) Output results from other HLA software platform. (No specific file extension).
    HLA_TYPE.add_argument("-rhped", help="\nOutput result(s) from other HLA software.\n\n", nargs='*')
    # (2) "*.hped"
    HLA_TYPE.add_argument("-hped", help="\nHLA Type Data not processed by \'NomenCleaner\' yet.\n\n")
    # (3) "*.chped"
    HLA_TYPE.add_argument("-chped", help="\nHLA Type Data processed by \'NomenCleaner\' (*.chped)\n\n")

    # (Meta-Analysis)
    parser.add_argument("--results-assoc", "-ra", help="\nResult files conducted by Association Tests(ex.\"*.assoc.logistic\").\n\n",
                        nargs='+')



    ##### < IMGT2Sequence > #####

    g_IMGT2Sequence = parser.add_argument_group(title='IMGT2Sequence', description='')

    g_IMGT2Sequence.add_argument("--imgt2sequence", help="\nGive this argument to implement \"IMGT2Sequence\" sub-module.\n\n", action='store_true')

    g_IMGT2Sequence.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n")

    g_IMGT2Sequence.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    g_IMGT2Sequence.add_argument("--no-multiprocess", help="\nSetting off parallel multiprocessing.\n\n", action='store_true')
    g_IMGT2Sequence.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')
    g_IMGT2Sequence.add_argument("--imgt-dir", help="\nIn case User just want to specify the directory of IMGT data folder.\n\n")



    ##### < HLA2MARKER > #####

    g_HLA2MARKER = parser.add_argument_group(title='HLA2MARKER', description='')

    g_HLA2MARKER.add_argument("--bmarkergenerator", help="\nGive this argument to implement \"HLA2MARKER\" sub-module.\n\n",
                              action='store_true')

    g_HLA2MARKER.add_argument("-dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n")
    g_HLA2MARKER.add_argument("-dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")



    ##### < NomenCleaner > #####

    g_NomenCleaner = parser.add_argument_group(title='NomenCleaner', description='')

    g_NomenCleaner.add_argument("--nomencleaner", help="\nGive this argument to implement \"NomenCleaner\" sub-module.\n\n",
                                action='store_true')

    # Additional input ped file type.
    PED_TYPE = g_NomenCleaner.add_mutually_exclusive_group()
    # PED_TYPE.add_argument("-ped", help="\nHLA Type Data.\n\n") # moved to main group.
    PED_TYPE.add_argument("-ped-Ggroup", help="\nHLA Type Data(G-group allele \"*.ped\" file).\n\n", dest="ped_G")
    PED_TYPE.add_argument("-ped-Pgroup", help="\nHLA Type Data(P-group allele \"*.ped\" file).\n\n", dest="ped_P")

    g_NomenCleaner.add_argument("-iat", help="\nIntegrated Allele Table file(*.iat).\n\n")

    # Output format selection
    format_selection = g_NomenCleaner.add_mutually_exclusive_group()


    format_selection.add_argument("--1field", help="\nOutput ped file as '1-field' format.\n\n",
                                        action="store_true", dest="oneF")
    format_selection.add_argument("--2field", help="\nOutput ped file as '2-field' format.\n\n",
                                        action="store_true", dest="twoF")
    format_selection.add_argument("--3field", help="\nOutput ped file as '3-field' format.\n\n",
                                        action="store_true", dest="threeF")
    format_selection.add_argument("--4field",
                                        help="\nOutput ped file as '4-field(Current Standard Names)' format.\n\n",
                                        action="store_true", dest="fourF")
    format_selection.add_argument("--G-group", help="\nOutput ped file as 'G-group' format.\n\n",
                                        action="store_true")
    format_selection.add_argument("--P-group", help="\nOutput ped file as 'P-group' format.\n\n",
                                        action="store_true")

    # Additional utility flags
    g_NomenCleaner.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*01:01\" -> \"01:01\").\n\n", action='store_true')



    ##### < HLA Analysis related > #####


    ### Association Test 1 (Logistic Regression by Plink 1.07)

    g_ASSOC1_Logistic = parser.add_argument_group(title='(Association Test 1) Logistic Regression.', description='')

    g_ASSOC1_Logistic.add_argument("--logistic", "-lr", help="\nGive this argument to implement \"Logistic Regression\". (Plink v1.07)\n\n", action='store_true')

    g_ASSOC1_Logistic.add_argument("--covar", help="\nSpecify .covar file (Plink v1.07).\n\n")
    g_ASSOC1_Logistic.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use. (Plink v1.07)\n\n")

    g_ASSOC1_Logistic.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.07).\n\n")
    g_ASSOC1_Logistic.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use. (Plink v1.07)\n\n")

    CondVars = g_ASSOC1_Logistic.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify Marker name(s) to condition as comma-separated(\",\"). (Plink v1.07)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify the Marker name(s) to condition as a file. (Plink v1.07)\n\n")

    g_ASSOC1_Logistic.add_argument("--reference-allele", help="\nSpecify the Reference Allele file. (Plink v1.07)\n\n")



    ### Association Test 2 (Omnibus Test by Buhm Han.)

    g_ASSOC2_Omnibus = parser.add_argument_group(title='(Association Test 2) Omnibus Test.', description='')
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

    # g_heatmap

    ### manhattan

    g_manhattan = parser.add_argument_group(title='(Plotting 2) Manhattan', description='')
    g_manhattan.add_argument("--manhattan", help="\nGenerate Manhattan Plot.\n\n", action="store_true")

    g_manhattan.add_argument("--point-color", "-pc", help="\nPoint color(ex. \"#778899\").\n\n", default="#778899")
    g_manhattan.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n", default="#FF0000")



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




    ########## < All in One Implementation. > ##########

    """
    [Indispensable Arguments]
    
    (1) Output Prefix
    (2) HLA-typed data (raw_hped, *.hped, *.chped)
    (3) Dictionary from IMGT2Sequence
    
    """

    ### < Core variables for indispensable arguments. > ###

    from src.HLA2HPED.HLA2HPED import HATK_HLA2HPED
    from src.IMGT2Sequence.IMGT2Seqeunce import HATK_IMGT2Sequence
    from src.NomenCleaner.NomenCleaner import NomenCleaner
    from src.b_MarkerGenerator.b_MarkerGenerator import HATK_b_MarkerGenerator
    import src.HLA_Analysis.HLA_Analysis as HLA_Analysis


    __dict_AA__ = args.dict_AA
    __dict_SNPS__ = args.dict_SNPS
    __IAT__ = args.iat
    # __MapTable__ = args # for heatmap.

    __HLA_type__ = None

    __b_MARKER_PANELS__ = None



    ########## < [0] : Indispensable Arguments > ##########


    ##### (1) Output Prefix

    if not bool(args.out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out"))
        sys.exit()

    print(std_MAIN_PROCESS_NAME + "Output prefix :\n\n{0}\n".format(args.out))




    ##### (2) HLA-type data (raw_hped, *.hped, *.chped)

    if bool(args.platform):

        print(std_MAIN_PROCESS_NAME + "HLA type data given as output results of other HLA software.\n".format(args.hped))

        # -rhped
        if bool(args.rhped):

            print("<Platform> : {0}\n".format(args.platform))

            for i in range(0, len(HLA_names)):
                print("{0} : {1}".format(HLA_names[i], args.rhped[i]))


            ### Call "HLA2HPED"

            __HLA_type__ = HATK_HLA2HPED(args.rhped, args.out, args.platform)


            ### Call "NomenCleaner" (with "IMGT2Sequence")

            if not bool(args.iat):

                ### Call IMGT2Sequence

                __dict_AA__, __dict_SNPS__, __IAT__ = HATK_IMGT2Sequence(args.hg, args.out, args.imgt, _imgt_dir=args.imgt_dir)

                print(std_MAIN_PROCESS_NAME + "(IMGT2Sequence)")
                print("Newly prepared \"AA_Dictionary\": {0}".format(__dict_AA__))
                print("Newly prepared \"SNPS_Dictionary\" : {0}".format(__dict_SNPS__))
                print("Newly prepared \"IAT\" for NomenCleaner : {0}".format(__IAT__))


            __HLA_type__ = NomenCleaner(__HLA_type__, 1, __IAT__, args.out, 4)

            print(std_MAIN_PROCESS_NAME + "(NomenCleaner)")
            print("\n\"Cleaned HPED(*.chped) file\" for b:MarkerGenerator : {0}".format(__HLA_type__))



        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "You specified the argument \"--platform\" to use \"HLA2HPED\", but the argument \"-rhped(--raw-hped)\" wansn't given.\n")
            sys.exit()

    else:

        """
        Actually, the arugment "args.hped", "args.chped" and "args.rhped" must not be given simultaneously.
        Checking this condition is undertaken by argparse "mutually_exclusive_group". 
        See "HLA_TYPE" argument group.
        """

        ## .hped
        if bool(args.hped):

            if os.path.exists(args.hped) and args.hped.endswith(".hped"):

                print(std_MAIN_PROCESS_NAME + "HLA type data as \"*.hped\" : {0}".format(args.hped))

                __HLA_type__ = args.hped


                ### Call "NomenCleaner" (with "IMGT2Sequence")

                if not bool(args.iat):

                    ### Call IMGT2Sequence

                    __dict_AA__, __dict_SNPS__, __IAT__ = HATK_IMGT2Sequence(args.hg, args.out, args.imgt, _imgt_dir=args.imgt_dir)

                    print(std_MAIN_PROCESS_NAME + "(IMGT2Sequence)")
                    print("Newly prepared \"AA_Dictionary\": {0}".format(__dict_AA__))
                    print("Newly prepared \"SNPS_Dictionary\" : {0}".format(__dict_SNPS__))

                print("Prepared \"IAT\" for NomenCleaner : {0}".format(__IAT__))

                __HLA_type__ = NomenCleaner(__HLA_type__, 1, __IAT__, args.out, 4)

                print("\nPrepared \"Cleaned HPED(*.chped) file\" for b:MarkerGenerator : {0}".format(__HLA_type__))


            else:

                print(std_ERROR_MAIN_PROCESS_NAME + "The file give with the argument \"-hped\" can't be found. Please check it again.\n")
                sys.exit()


        ## .chped
        elif bool(args.chped):

            if os.path.exists(args.chped) and args.chped.endswith(".chped"):

                print(std_MAIN_PROCESS_NAME + "HLA type data given as \"*.chped\" : {0}".format(args.chped))

                __HLA_type__ = args.chped

                print("\nPrepared \"Cleaned HPED(*.chped) file\" for b:MarkerGenerator : {0}".format(__HLA_type__))

            else:

                print(std_ERROR_MAIN_PROCESS_NAME + "The file give with the argument \"-chped\" can't be found. Please check it again.\n")
                sys.exit()



        ## --input
        else:
            # Given as prefix with "--input"

            if os.path.exists(args.input + ".hped") and os.path.exists(args.input + ".chped"):
                print(std_ERROR_MAIN_PROCESS_NAME + "There are both \"*.hped\" and \"*.chped\". The program can't decide what to use. Please check them again.\n")
                sys.exit()

            if os.path.exists(args.input + ".hped"):
                print(std_MAIN_PROCESS_NAME + "HLA type data(hped) was given as the prefix by \"--input\" argument.\n")
                print("Given hped file : {0}\n".format(args.input + ".hped"))

                __HLA_type__ = args.input + ".hped"


                ### Call "NomenCleaner" (with "IMGT2Sequence")

                if not bool(args.iat):

                    ### Call IMGT2Sequence

                    __dict_AA__, __dict_SNPS__, __IAT__ = HATK_IMGT2Sequence(args.hg, args.out, args.imgt, _imgt_dir=args.imgt_dir)

                    print(std_MAIN_PROCESS_NAME + "(IMGT2Sequence)")
                    print("Newly prepared \"AA_Dictionary\": {0}".format(__dict_AA__))
                    print("Newly prepared \"SNPS_Dictionary\" : {0}".format(__dict_SNPS__))

                print("Prepared \"IAT\" for NomenCleaner : {0}".format(__IAT__))

                __HLA_type__ = NomenCleaner(__HLA_type__, 1, __IAT__, args.out, 4)

                print("\nPrepared \"Cleaned HPED(*.chped) file\" for b:MarkerGenerator : {0}".format(__HLA_type__))




            elif os.path.exists(args.input + ".chped"):

                print(std_MAIN_PROCESS_NAME + "HLA type data(chped) was given as the prefix by \"--input\" argument.\n")
                print("Given hped file : {0}\n".format(args.input + ".chped"))

                __HLA_type__ = args.input + ".chped"

                print("\nPrepared \"Cleaned HPED(*.chped) file\" for b:MarkerGenerator : {0}".format(__HLA_type__))



            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "HLA Type data can't be found anyhow. Please read the manual or documentation for HATK first.\n")
                sys.exit()





    ##### (3) Dictionaries

    # HLA Sequence information(Dictionary)
    if not(bool(__dict_AA__) and bool(__dict_SNPS__)):

        print(std_WARNING_MAIN_PROCESS_NAME + "Sequence dictionary files are not prepared completely.\n"
                                              "Implementing \"IMGT2Sequence\" forcibly.\n")

        ### Call IMGT2Sequence

        __dict_AA__, __dict_SNPS__, __IAT__ = HATK_IMGT2Sequence(args.hg, args.out, args.imgt, _imgt_dir=args.imgt_dir)

        print(std_MAIN_PROCESS_NAME + "(IMGT2Sequence)")
        print("Newly prepared \"AA_Dictionary\": {0}".format(__dict_AA__))
        print("Newly prepared \"SNPS_Dictionary\" : {0}".format(__dict_SNPS__))



    if bool(__dict_AA__) and not (os.path.exists(__dict_AA__ + ".txt") and os.path.exists(__dict_AA__ + ".map")):

        print(std_ERROR_MAIN_PROCESS_NAME + "The file give by the argument \"-dict-AA\" doesn't exist. Please check it again.\n")
        sys.exit()

    if bool(__dict_SNPS__) and not (os.path.exists(__dict_SNPS__ + ".txt") and os.path.exists(__dict_SNPS__ + ".map")):

        print(std_ERROR_MAIN_PROCESS_NAME + "The file give by the argument \"-dict-SNPS\" doesn't exist. Please check it again.\n")
        sys.exit()




    print(std_MAIN_PROCESS_NAME + "SUMMARY - HLA Sequence Information(IMGT2Sequence)\n")
    print("< HLA-type Data >")
    print("chped : {0}".format(__HLA_type__))
    print("< AA >")
    print("seq : {0}".format(__dict_AA__+".txt"))
    print("map : {0}".format(__dict_AA__+".map"))
    print("< SNPS >")
    print("seq : {0}".format(__dict_SNPS__+".txt"))
    print("map : {0}".format(__dict_SNPS__+".map"))
    print("< IAT >")
    print("iat : {0}".format(__IAT__))
    # print("< Maptable file >")



    ### 여기다가 이제 b_MarkerGenerator, Association Test, Plotting이렇게만 뭉치면 됨.

    ########## < [1] : b:MarkerGenerator > ##########


    __b_MARKER_PANELS__ = HATK_b_MarkerGenerator(__HLA_type__, args.out, args.hg, __dict_AA__, __dict_SNPS__)

    print(std_MAIN_PROCESS_NAME + "(b:MarkerGenerator)")
    print(__b_MARKER_PANELS__)


    ########## < [2] : Association Test > ##########

    # Logistic Regression

    __Logistic_Result__ = [HLA_Analysis.HATK_ASSOC1_Logistic_Regression(args.input, item, _b_marker_panels=item) for item in __b_MARKER_PANELS__]
    print("\nLogistic_Reuslt : {0}".format(__Logistic_Result__))


    ########## < [3] : Plotting > ##########

    # Heatmap



    # Manhattan




    # if not bool(args.input) and not (args.imgt2sequence or args.hla2marker or args.nomencleaner or args.manhattan):
    #     print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--input"))
    #     sys.exit()

    # (2018. 10. 04.) Turning off it for a while.
    # if not (args.imgt2sequence or args.hla2marker or args.nomencleaner or args.manhattan or args.heatmap or args.converter):
    #     # if none of flag for sub-module is given, then it is false state.
    #     print(std_ERROR_MAIN_PROCESS_NAME + "You should give at least on flag for sub-module to use.\n")
    #     sys.exit()





    # if args.imgt2sequence:
    #
    #     ##### IMGT2Sequence #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing IMGT2Sequence.")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. -hg
    #     2. -o (*)
    #     3. -imgt
    #     """
    #
    #
    #     from src.IMGT2Sequence.IMGT2Seqeunce import HATK_IMGT2Sequence
    #
    #     __dict_AA__, __dict_SNPS__, _t_IAT = HATK_IMGT2Sequence(_HG=args.hg, _OUTPUT=args.out, _IMGT=args.imgt,
    #                                                             _no_Indel=args.no_indel, _no_MultiP=args.no_multiprocess,
    #                                                             _save_intermediates=args.save_intermediates,
    #                                                             _imgt_dir=args.imgt_dir)
    #
    #
    #
    #
    # if args.nomencleaner:
    #
    #     ##### NomenCleaner #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing NomenCleaner.")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. either -ped, -ped-Ggroup or -ped-Pgroup
    #     2. -iat
    #     3. -o (*)
    #     4. output format(ex. --1field, etc.)
    #
    #     (optionals)
    #     5. --No-caption
    #     """
    #
    #     t_ped = ""
    #     t_ped_descriptor = -1
    #     t_FILE_FORMAT = -1
    #
    #
    #     if not(bool(args.ped) or bool(args.ped_Ggroup) or bool(args.ped_Pgroup)):
    #         print(std_ERROR_MAIN_PROCESS_NAME + "You must give at least one argument among \"-ped\", \"-ped-Ggroup\", \"-ped-Pgroup\".\n")
    #         sys.exit()
    #     else:
    #
    #         if bool(args.ped):
    #             t_ped_descriptor = 1
    #             t_ped = args.ped
    #         elif bool(args.ped_Ggroup):
    #             t_ped_descriptor = 2
    #             t_ped = args.ped_Ggroup
    #         elif bool(args.ped_Pgroup):
    #             t_ped_descriptor = 3
    #             t_ped = args.ped_Pgroup
    #         else:
    #             print(std_ERROR_MAIN_PROCESS_NAME + "The arguments related to ped(\"-ped\", \"-ped-Ggroup\" or \"-ped-Pgroup\") has wrong values."
    #                                                 "Please check them again.\n")
    #             sys.exit()
    #
    #
    #
    #     if not bool(args.iat):
    #         print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-iat"))
    #         sys.exit()
    #
    #     if not(bool(args.oneF) or bool(args.twoF) or bool(args.threeF) or bool(args.fourF) or bool(args.G_group) or bool(args.P_group)):
    #         print(std_ERROR_MAIN_PROCESS_NAME + "The argument related to field format has not given. Please check it again.\n")
    #         sys.exit()
    #     else:
    #         t_FILE_FORMAT = (1 if bool(args.oneF) else 2 if bool(args.twoF) else 3 if bool(args.threeF) else 4
    #                             if bool(args.fourF) else 5 if bool(args.G_group) else 6 if bool(args.P_group) else -1)
    #
    #
    #     if t_ped_descriptor == 2 and t_FILE_FORMAT == 5:
    #         print(std_ERROR_MAIN_PROCESS_NAME + "Pointless transformation. (Transformation G-group to G-group is meaningless.)")
    #         print("Skip this Transformation Request.\n")
    #         sys.exit()
    #
    #     if t_ped_descriptor == 3 and t_FILE_FORMAT == 6:
    #         print(std_ERROR_MAIN_PROCESS_NAME + "Pointless transformation. (Transformation P-group to P-group is meaningless.)")
    #         print("Skip this Transformation Request.\n")
    #         sys.exit()
    #
    #
    #     from src.NomenCleaner.NomenCleaner import NomenCleaner
    #
    #     NomenCleaner(_p_ped=t_ped, _ped_descriptor=t_ped_descriptor, _p_iat=args.iat, _out=args.out,
    #                  _field_format=t_FILE_FORMAT, _f_NoCaption=bool(args.NoCaption))
    #
    #
    #
    #
    # if args.bmarkergenerator:
    #
    #     ##### b:MarkerGenerator #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing b:MarkerGenerator.")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. -input (*)
    #     2. -hped or -ped
    #     3. -hg
    #     4. -o (*)
    #     5. --dict-AA
    #     6. --dict-SNPS
    #     """
    #
    #     from src.b_MarkerGenerator.b_MarkerGenerator import HATK_b_MarkerGenerator
    #
    #     HATK_b_MarkerGenerator(HPED=args.hped, _OUT=args.out, _hg=args.hg,
    #                            _dictionary_AA=__dict_AA__, _dictionary_SNPS=__dict_SNPS__,
    #                            _plain_SNP_DATA=args.input)
    #
    #
    #
    #
    # if args.logistic:
    #
    #     ##### (Association Test 1) logistic regression(by Plink 1.07) #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing Logistic Regression(Association Test, plink 1.07).")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. -input(--bfile) (*)
    #     2. -o (*)
    #
    #     (optionals)
    #     3. -covar (with -covar-name)
    #     4. -phe (with -phe-name)
    #     5. either --condition or --condition-list
    #     6. --reference-allele
    #
    #     """
    #
    #     if not bool(args.covar) and bool(args.covar_name):
    #         print(std_WARNING_MAIN_PROCESS_NAME + 'You gave "{0}" argument while "{1}" wasn\'t given. The argument "{0}" will be ignored.\n'.format("--covar-name", "--covar"))
    #         args.covar_name = None
    #
    #     if not bool(args.pheno) and bool(args.pheno_name):
    #         print(std_WARNING_MAIN_PROCESS_NAME + 'You gave "{0}" argument while "{1}" wasn\'t given. The argument "{0}" will be ignored.\n'.format("--pheno-name", "--pheno"))
    #         args.pheno_name = None
    #
    #
    #
    #     from src.HLA_Analysis.HLA_Analysis import ASSOC1_Logistic_Regression
    #
    #     t_result_logistic = ASSOC1_Logistic_Regression(args.input, args.out, args.covar, args.covar_name, args.pheno, args.pheno_name,
    #                                                    args.condition, args.condition_list, args.reference_allele)
    #
    #
    #
    # if args.omnibus:
    #
    #     ##### (Association Test 2) Omnibus Test #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing Omnibus Test(Association Test).")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. --input
    #     2. -o (*)
    #     3. --phased
    #     4. --pheno
    #     5. --pheno-name
    #     6. --covar
    #     7. --covar-names
    #     8. --rare-threshold (2018. 10. 04.) Postponed.
    #     9. --condition
    #     10. --condition-list
    #
    #     """
    #
    #     # --phased(.aa)
    #     if not bool(args.phased):
    #         print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--phased"))
    #         sys.exit()
    #
    #     # --pheno
    #     if not bool(args.pheno):
    #         print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--pheno"))
    #         sys.exit()
    #     else:
    #         # Phenotype file is given but the name of phenotype to use isn't given.
    #         if not bool(args.pheno_name):
    #             print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--pheno-name"))
    #             sys.exit()
    #
    #
    #     # (2018. 10. 04.) Postponed.
    #     # # --rare-threshold
    #     # if not bool(args.rare_threshold):
    #     #     print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--rare-threshold"))
    #     #     sys.exit()
    #
    #
    #     from src.HLA_Analysis.HLA_Analysis import ASSOC2_Omnibus_Test
    #
    #     ASSOC2_Omnibus_Test(_input=args.input, _out=args.out, _phased=args.phased, _phe=args.pheno, _phe_name=args.pheno_name,
    #                         _covar=args.covar, _covar_names=args.covar_name)
    #
    #
    # if args.meta_analysis:
    #
    #     ##### Meta-Analysis #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing Meta-Analysis.\n")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. -o (*)
    #     2. -ra
    #
    #     """
    #
    #     if not bool(args.results_assoc):
    #         print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--results-assoc(-ra)"))
    #         sys.exit()
    #
    #
    #
    # if args.manhattan:
    #
    #     ##### manhattan #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing Manhattan(Plotting).")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. --result-assoc
    #     2. --out (*)
    #     3. -hg
    #
    #     (optionals)
    #     4. --point-color
    #     5. --top-color
    #     """
    #
    #     if not bool(args.results_assoc):
    #         print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--results-assoc(-ra)"))
    #         sys.exit()
    #
    #     if not bool(args.hg):
    #         print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hg"))
    #         sys.exit()
    #
    #
    #     from src.HLA_Analysis.Plotting.manhattanplot.manhattan import manhattan
    #
    #     manhattan(_results_assoc=args.results_assoc, _out=args.out, _hg=args.hg,
    #               _pointcol=args.point_color, _topcol=args.top_color)
    #
    #
    #
    # if args.heatmap:
    #
    #     ##### heatmap #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing Heatmap(Plotting).")
    #
    #     """
    #     List of necessary arguments.
    #     """
    #
    #
    # # if args.converter:
    # #
    # #     ##### Converter #####
    # #
    # #     print(std_MAIN_PROCESS_NAME + "Implementing Converter.")
    # #
    # #     """
    # #     List of necessary arguments.
    # #     """
    #
    #
    #
