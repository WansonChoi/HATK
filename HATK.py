#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))





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

    parser.add_argument("-hped", help="\nHLA Type Data processed by \'Nomencleaner\'(.hped)\n\n")

    parser.add_argument("--results-assoc", "-ra", help="\nResult files conducted by Association Tests(ex.\"*.assoc.logistic\").\n\n",
                        nargs='+')



    ##### < IMGT2Sequence > #####

    g_IMGT2Sequence = parser.add_argument_group(title='IMGT2Sequence', description='')

    g_IMGT2Sequence.add_argument("--imgt2sequence", help="\nGive this argument to implement \"IMGT2Sequence\" sub-module.\n\n", action='store_true')

    g_IMGT2Sequence.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n")

    g_IMGT2Sequence.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    g_IMGT2Sequence.add_argument("--no-multiprocess", help="\nSetting off parallel multiprocessing.\n\n", action='store_true')



    ##### < HLA2MARKER > #####

    g_HLA2MARKER = parser.add_argument_group(title='HLA2MARKER', description='')

    g_HLA2MARKER.add_argument("--hla2marker", help="\nGive this argument to implement \"HLA2MARKER\" sub-module.\n\n",
                              action='store_true')

    g_HLA2MARKER.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n")
    g_HLA2MARKER.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")



    ##### < NomenCleaner > #####

    g_NomenCleaner = parser.add_argument_group(title='NomenCleaner', description='')

    g_NomenCleaner.add_argument("--nomencleaner", help="\nGive this argument to implement \"NomenCleaner\" sub-module.\n\n",
                                action='store_true')

    # Additional input ped file type.
    PED_TYPE = g_NomenCleaner.add_mutually_exclusive_group()
    PED_TYPE.add_argument("-ped", help="\nHLA Type Data.\n\n")
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
    g_NomenCleaner.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*\").\n\n", action='store_true')



    ##### < HLA Analysis related > #####


    ### Association Test 1 (Logistic Regression by Plink 1.07)

    g_ASSOC1_Logistic = parser.add_argument_group(title='(Association Test 1) Logistic Regression.', description='')

    g_ASSOC1_Logistic.add_argument("--logistic", "-lr", help="\nGive this argument to implement \"HLA-Analysis\" sub-module. (Plink v1.07)\n\n", action='store_true')

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

    g_ASSOC2_Omnibus.add_argument("--phased", "-ph", help="\nSpecify the \".aa\" file.(for Omnibus Test).\n\n")
    g_ASSOC2_Omnibus.add_argument("--rare-threshold", "-rth", help="\nSpecify the \".aa\" file.(for Omnibus Test).\n\n")


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



    ##### < Converter > #####

    g_Converter = parser.add_argument_group(title='Converter', description='')

    g_Converter.add_argument("--converter", help="\nGive this argument to implement \"Converter\" sub-module.\n\n", action='store_true')

    g_Converter.add_argument("--AXIOM", help="\nAXIOM output file format.\n\n", action="store_true")
    g_Converter.add_argument("--HIBAG", help="\nHIBAG output file format.\n\n", action="store_true")
    g_Converter.add_argument("--xHLA", help="\nxHLA output file format.\n\n", action="store_true")









    ##### < for Testing > #####

    # args = parser.parse_args(["--imgt2sequence", "-imgt", "370", "-o", "TEST/TEST", "-hg", "18"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)




    ########## < Argument Chekcing > ##########

    ### < Temporary variables for argument checking > ###

    _t_dict_AA = None
    _t_dict_SNPS = None


    ### Checking indispensable common arguments

    if not bool(args.out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out"))
        sys.exit()

    if not bool(args.input) and not (args.imgt2sequence or args.hla2marker or args.nomencleaner or args.manhattan):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--input"))
        sys.exit()

    if not (args.imgt2sequence or args.hla2marker or args.nomencleaner or args.manhattan or args.heatmap or args.converter):
        # if none of flag for sub-module is given, then it is false state.
        print(std_ERROR_MAIN_PROCESS_NAME + "You should give at least on flag for sub-module to use.\n")
        sys.exit()





    if args.imgt2sequence:

        ##### IMGT2Sequence #####

        print(std_MAIN_PROCESS_NAME + "Implementing IMGT2Sequence.")

        """
        List of necessary arguments.

        1. -hg 
        2. -o (*)
        3. -imgt
        """

        if not bool(args.hg):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hg"))
            sys.exit()

        if not bool(args.imgt):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-imgt"))
            sys.exit()


        from src.IMGT2Sequence.IMGT2Seqeunce import MakeDictionary

        _t_dict_AA, _t_dict_SNPS = MakeDictionary(_HG=args.hg, _OUTPUT=args.out, _IMGT=args.imgt, _no_Indel=args.no_indel, _no_MultiP=args.no_multiprocess)



    if args.hla2marker:

        ##### HLA2MARKER #####

        print(std_MAIN_PROCESS_NAME + "Implementing HLA2MARKER.")

        """
        List of necessary arguments.

        1. -input (*)
        2. -hped or -ped
        3. -hg
        4. -o (*)
        5. --dict-AA
        6. --dict-SNPS
        """

        if not bool(args.hped):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hped"))
            sys.exit()
        else:

            # Checking whether it went through "NomenCleaner.py".

            if not args.hped.endswith(".hped"):
                print(std_ERROR_MAIN_PROCESS_NAME + 'Given ped file should be processed by "NomenCleaner.py". '
                                                    '(The file extension of its output is ".hped".)\n')

                """
                Add code to implement "NomenCleaner.py" here.
                """
                sys.exit()

        if not bool(args.hg):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hg"))
            sys.exit()


        if bool(args.imgt2sequence) and (bool(_t_dict_AA) and bool(_t_dict_SNPS)):

            # When "--imgt2sequence" given with "--hla2marker".
            # "--dict-AA" and "--dict-SNPS" will be overridden.
            print(std_MAIN_PROCESS_NAME + "Using newly generated Dictionary for \"--dict-AA\" and \"--dict-SNPS\".\n")

        else:

            if bool(args.dict_AA) and bool(args.dict_SNPS):
                _t_dict_AA = args.dict_AA
                _t_dict_SNPS = args.dict_SNPS
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + 'One of "--dict-AA" or "--dict-SNPS" options(or Both) are not given. Please check them again.\n')
                sys.exit()


        from src.HLA2MARKER.HLA2MARKER import MakeReference

        MakeReference(_HLA_ped=args.hped, _OUT=args.out, _hg=args.hg,
                      _dictionary_AA=_t_dict_AA, _dictionary_SNPS=_t_dict_SNPS,
                      _plain_SNP_DATA=args.input)



    if args.nomencleaner:

        ##### NomenCleaner #####

        print(std_MAIN_PROCESS_NAME + "Implementing NomenCleaner.")

        """
        List of necessary arguments.

        1. either -ped, -ped-Ggroup or -ped-Pgroup
        2. -iat
        3. -o (*)
        4. output format(ex. --1field, etc.)

        (optionals)
        5. --No-caption
        """

        t_ped = ""
        t_ped_descriptor = -1
        t_FILE_FORMAT = -1


        if not(bool(args.ped) or bool(args.ped_Ggroup) or bool(args.ped_Pgroup)):
            print(std_ERROR_MAIN_PROCESS_NAME + "You must give at least one argument among \"-ped\", \"-ped-Ggroup\", \"-ped-Pgroup\".\n")
            sys.exit()
        else:

            if bool(args.ped):
                t_ped_descriptor = 1
                t_ped = args.ped
            elif bool(args.ped_Ggroup):
                t_ped_descriptor = 2
                t_ped = args.ped_Ggroup
            elif bool(args.ped_Pgroup):
                t_ped_descriptor = 3
                t_ped = args.ped_Pgroup
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "The arguments related to ped(\"-ped\", \"-ped-Ggroup\" or \"-ped-Pgroup\") has wrong values."
                                                    "Please check them again.\n")
                sys.exit()



        if not bool(args.iat):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-iat"))
            sys.exit()

        if not(bool(args.oneF) or bool(args.twoF) or bool(args.threeF) or bool(args.fourF) or bool(args.G_group) or bool(args.P_group)):
            print(std_ERROR_MAIN_PROCESS_NAME + "The argument related to field format has not given. Please check it again.\n")
            sys.exit()
        else:
            t_FILE_FORMAT = (1 if bool(args.oneF) else 2 if bool(args.twoF) else 3 if bool(args.threeF) else 4
                                if bool(args.fourF) else 5 if bool(args.G_group) else 6 if bool(args.P_group) else -1)


        if t_ped_descriptor == 2 and t_FILE_FORMAT == 5:
            print(std_ERROR_MAIN_PROCESS_NAME + "Pointless transformation. (Transformation G-group to G-group is meaningless.)")
            print("Skip this Transformation Request.\n")
            sys.exit()

        if t_ped_descriptor == 3 and t_FILE_FORMAT == 6:
            print(std_ERROR_MAIN_PROCESS_NAME + "Pointless transformation. (Transformation P-group to P-group is meaningless.)")
            print("Skip this Transformation Request.\n")
            sys.exit()


        from src.NomenCleaner.NomenCleaner import NomenCleaner

        NomenCleaner(_p_ped=t_ped, _ped_descriptor=t_ped_descriptor, _p_iat=args.iat, _out=args.out,
                     _field_format=t_FILE_FORMAT, _f_NoCaption=bool(args.NoCaption))



    # elif args.hla_analysis:
    #
    #     ##### HLA_Analysis #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing HLA_Analysis.")
    #
    #     """
    #     List of necessary arguments.
    #
    #     1. -input(--bfile)
    #     2. -o (*)
    #     3. either -lr or -ob
    #
    #     (optionals)
    #     4. -covar (with -covar-name)
    #     5. -phe (with -phe-name)
    #     6. either --condition or --condition-list
    #     7. --reference-allele
    #     8. threshold
    #     9. phased.
    #
    #     """


    if args.meta_analysis:

        ##### Meta-Analysis #####

        print(std_MAIN_PROCESS_NAME + "Implementing Meta-Analysis.\n")

        """
        List of necessary arguments.

        1. -o (*)
        2. -ra

        """

        if not bool(args.results_assoc):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--results-assoc(-ra)"))
            sys.exit()



    if args.manhattan:

        ##### manhattan #####

        print(std_MAIN_PROCESS_NAME + "Implementing Manhattan(Plotting).")

        """
        List of necessary arguments.
        
        1. --result-assoc
        2. --out (*)
        3. -hg
        
        (optionals)
        4. --point-color
        5. --top-color
        """

        if not bool(args.results_assoc):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--results-assoc(-ra)"))
            sys.exit()

        if not bool(args.hg):
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hg"))
            sys.exit()


        from src.HLA_Analysis.Plotting.manhattanplot.manhattan import manhattan

        manhattan(_results_assoc=args.results_assoc, _out=args.out, _hg=args.hg,
                  _pointcol=args.point_color, _topcol=args.top_color)



    if args.heatmap:

        ##### heatmap #####

        print(std_MAIN_PROCESS_NAME + "Implementing Heatmap(Plotting).")

        """
        List of necessary arguments.
        """


    # if args.converter:
    #
    #     ##### Converter #####
    #
    #     print(std_MAIN_PROCESS_NAME + "Implementing Converter.")
    #
    #     """
    #     List of necessary arguments.
    #     """



