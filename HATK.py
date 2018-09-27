#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

"""

    Main Interface for HLA Analysis Toolkit(HATK).
    Users will implement next function through this script.
    
    (1) HLA-Analysis
    (2) MakeRefernce
    (3) NomenCleaner
    (4) MakeDictionary
    (5) Converter
    

"""





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

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')



    ########## < Creating sub-parser > ##########
    subparsers = parser.add_subparsers(title="SUB-COMMANDS", metavar="{SUB}", dest="subparser_name")






    ##### sub-parser for 'MakeDictionary'
    parser_MAKEDICTIONARY = subparsers.add_parser('MAKEDICTIONARY',
                                                  help='MakeDictionary help\n\n',
                                                  formatter_class=argparse.RawTextHelpFormatter,
                                                  add_help=False,
                                                  description=textwrap.dedent('''\
    ###########################################################################################
        
        < MakeDictionary >
    
        explanation...
    
    ###########################################################################################
    '''))

    parser_MAKEDICTIONARY._optionals.title = "OPTIONS"
    parser_MAKEDICTIONARY.add_argument("-h", "--help", help="\nShow this help message for 'MAKEDICTIONARY' and exit\n\n", action='help')
    parser_MAKEDICTIONARY.add_argument("-hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                                       metavar="HG", required=True)
    # cf) It was "dest" parameter that setting namespace.

    parser_MAKEDICTIONARY.add_argument("-o", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser_MAKEDICTIONARY.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                                       required=True)

    parser_MAKEDICTIONARY.add_argument("--type", "-t", help="\nSNPS or AA or Both.\n\n", choices=["AA", "SNPS", "BOTH"],
                                       metavar="TYPE", default="BOTH")








    ##### sub-parser for 'MakeReference'
    parser_MAKEREFERENCE = subparsers.add_parser('HLA2MARKER',
                                                  help='HLA2MARKER help\n\n',
                                                  formatter_class=argparse.RawTextHelpFormatter,
                                                  add_help=False,
                                                  description=textwrap.dedent('''\
    #################################################################################################
        
        < HLA2MARKER.py >
        
        This script helps to prepare a reference panel for HLA region.

        (ex.)

        : python3 HLA2MARKER.py  
            -ped ./data/MakeReference/HAPMAP_CEU_HLA.4field.ped 
            -hg 18 
            -o ./Trial_HAPMAP_CEU
            -dict-AA ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370
            -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370 

        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
    '''))

    parser_MAKEREFERENCE._optionals.title = "OPTIONS"

    parser_MAKEREFERENCE.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser_MAKEREFERENCE.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n", required=True)
    parser_MAKEREFERENCE.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
    parser_MAKEREFERENCE.add_argument("-o", help="\nOutput file prefix\n\n")

    hla_dict = parser_MAKEREFERENCE.add_argument_group(title='HLA_DICTIONARY',
                                         description='- Arguments to specify HLA_DICTIONARY Information to New version of MakeReference(2.0)\n'
                                                     '- If you\'re going to use previous version of MakeReference, then Don\'t care about these options.')

    hla_dict.add_argument("-dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", default="Not_given")
    hla_dict.add_argument("-dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", default="Not_given")


    # parser.add_argument("-dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", default="Not_given")
    # parser.add_argument("-dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", default="Not_given")



    # ##### sub-parser for 'MakeReference'
    # parser_MAKEREFERENCE = subparsers.add_parser('MAKEREFERENCE',
    #                                               help='MakeReference help\n\n',
    #                                               formatter_class=argparse.RawTextHelpFormatter,
    #                                               add_help=False,
    #                                               description=textwrap.dedent('''\
    # #################################################################################################
    #
    #     < MakeReference.py >
    #
    #     This script helps prepare a reference dataset for HLA imputation
    #
    #     Usage(1)
    #     : python3 MakeReference.py --previous-version -i ./data/MakeReference_old/HAPMAP_CEU
    #         -ped ./data/MakeReference_old/HAPMAP_CEU_HLA.ped -hg 18 -o ./Trial_HAPMAP_CEU
    #
    #     Usage(2)
    #     : python3 MakeReference.py -i ./data/MakeReference/HAPMAP_CEU
    #         -ped ./data/MakeReference/HAPMAP_CEU_HLA.4field.ped -hg 18 -o ./Trial_HAPMAP_CEU
    #         -dict-AA ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt
    #         -dict-AA-map ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map
    #         -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt
    #         -dict-SNPS-map ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map
    #
    #     HLA PED file should contain HLA alleles in the following (alphabetical) order:
    #     HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1
    #
    # #################################################################################################
    # '''))
    #
    # parser_MAKEREFERENCE._optionals.title = "OPTIONS"
    #
    # parser_MAKEREFERENCE.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
    #
    # parser_MAKEREFERENCE.add_argument("-i", help="\nInput Data file(.bed/.bim/.fam)\n\n", required=True)
    # parser_MAKEREFERENCE.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n", required=True)
    # parser_MAKEREFERENCE.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
    # parser_MAKEREFERENCE.add_argument("-o", help="\nOutput file prefix\n\n")
    #
    # parser_MAKEREFERENCE.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as old version.\n\n", action='store_true')
    # parser_MAKEREFERENCE.add_argument("--only-markers", help="\nGenerate a panel for only markers, not entire reference panel.\n\n", action='store_true') # (2018. 8. 1.) Introduced.
    #
    # hla_dict = parser_MAKEREFERENCE.add_argument_group(title='HLA_DICTIONARY',
    #                                      description='- Arguments to specify HLA_DICTIONARY Information to New version of MakeReference(2.0)\n'
    #                                                  '- If you\'re going to use previous version of MakeReference, then Don\'t care about these options.')
    #
    # hla_dict.add_argument("-dict-AA", help="\nInput HLA Dictionary file for AA Information.\n\n", default="Not_given")
    # hla_dict.add_argument("-dict-AA-map", help="\nInput HLA Dictionary .map file for AA Information.\n\n", default="Not_given")
    # hla_dict.add_argument("-dict-SNPS", help="\nInput HLA Dictionary file for SNPS Information.\n\n", default="Not_given")
    # hla_dict.add_argument("-dict-SNPS-map", help="\nInput HLA Dictionary .map file for SNPS Information\n\n", default="Not_given")




    # (2018. 6. 5) sub-parser for 'SNP2HLA' also will be introduced here, soon. => Prof. Han didn't want it. Canceled.




    ##### sub-parser for 'HLA-Analysis'
    parser_HLA_ANALYSIS = subparsers.add_parser('HLAANALYSIS',
                                                help="HLA_analysis help\n\n",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                add_help=False,
                                                description=textwrap.dedent('''\
    ###########################################################################################
        
        <HLA-Analysis>

            under development.     

    ###########################################################################################
    '''))

    parser_HLA_ANALYSIS._optionals.title = "OPTIONS"

    parser_HLA_ANALYSIS.add_argument("-h", "--help", help="\nShow this help message for 'HLAANALYSIS' and exit\n\n", action='help')


    ### < Arguments group for Input files. > ###
    group_INPUTS = parser_HLA_ANALYSIS.add_argument_group(title="INPUTs",
                                                          description=
    '---------------------------------------------------------------------------------\n'
    '- First, Give necessary input information just like using Plink.\n'
    '- Next, Conduct a Associaiton Test or Meta-Analysis.\n'
    '---------------------------------------------------------------------------------\n')

    # Necessary inputs for common.
    PEDorBED = group_INPUTS.add_mutually_exclusive_group()
    PEDorBED.add_argument("--bfile", "-bf", help="\nPrefix of \".bed\", \".fam\" and \".bim\" files (Plink v1.07).\n\n")
    PEDorBED.add_argument("--input", "-i", help="\nPrefix of \".bed\", \".covar\", or \".phe\" etc, when common prefix is prepared.\n\n")
    # (2018. 8. 7.) "--input" option deals with common prefix of (1) ".bed", (2) ".covar", (3) ".phe", (4) ".aa(phased)".

    group_INPUTS.add_argument("--out", "-o", help="\nSpecify the output prefix (Plink v1.07).\n\n\n\n\n")

    # Additional inputs if available.
    group_INPUTS.add_argument("--covar", help="\nSpecify .covar file (Plink v1.07).\n\n")
    group_INPUTS.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use.(Plink v1.07).\n\n")

    group_INPUTS.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.07).\n\n")
    group_INPUTS.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use.(Plink v1.07).\n\n")

    CondVars = group_INPUTS.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify Marker name(s) to condition as comma-separated(\",\") (Plink v1.07).\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify the Marker name(s) to condition as a file (Plink v1.07).\n\n")

    group_INPUTS.add_argument("--reference-allele", help="\nSpecify the Reference Allele file (Plink v1.07).\n\n")


    # Only for Omnibus Test
    group_INPUTS.add_argument("--phased", "-ph", help="\nSpecify the \".aa\" file.(for Omnibus Test).\n\n")
    group_INPUTS.add_argument("--rare-threshold", "-rth", help="\nSpecify the \".aa\" file.(for Omnibus Test).\n\n")

    # Only for Meta-Analysis
    group_INPUTS.add_argument("--rassoc", "-ra", help="\nSpecify the Result file(s) of association test for Meta-Analysis\n"
                                                    "(ex. \"*.asssoc.logistic\", etc.)\n\n", nargs='+')



    ### < Arguments group for Logistic Regression. > ###
    group_ASSOCIATION_TEST = parser_HLA_ANALYSIS.add_argument_group(title="ASSOCIATION TEST",
                                                                    description=
    '---------------------------------------------------------------------------------\n'
    '(1) Logistic Regression, \n'
    '(2) Omnibus Test. \n'
    '---------------------------------------------------------------------------------\n')

    group_ASSOCIATION_TEST.add_argument("--logistic-regression", "-lr", help="\nConducting Logistic Regression(Plink v1.07).\n\n",
                                        action="store_true")
    group_ASSOCIATION_TEST.add_argument("--omnibus-test", "-om", help="\nConducting Omnibus Test.\n\n",
                                        action="store_true")



    # (2018. 8. 2.) Later, Introduce the way that argument can take ped and map files seperately.


    # Later, by using exclusive-nor, Don't allow using "Association Test", "Meta-analysis" and "plotting"(maybe this can be fine.).




    ##### < Arguments for Meta-Analysis. > #####
    group_METAANALYSIS = parser_HLA_ANALYSIS.add_argument_group(title="META ANALYSIS",
                                                                description=
    '---------------------------------------------------------------------------------\n'
    '- Give the multiple \".ped\" files.\n'
    '---------------------------------------------------------------------------------\n')

    group_METAANALYSIS.add_argument("--meta-analysis", "-meta", help="\nConducting Meta Analysis(Plink v1.07)\n\n",
                                    action="store_true")

    # group_METAANALYSIS.add_argument("-i", help="\nThe results of association tests.\n\n", nargs='+')








    ##### sub-parser for 'PLOTTING'
    parser_PLOTTING = subparsers.add_parser('PLOTTING',
                                           help="The modules to plot association test results.\n\n",
                                           formatter_class=argparse.RawTextHelpFormatter,
                                           add_help=False,
                                           description=textwrap.dedent('''\
    ###########################################################################################

        <PLOTTING>

        Under construction.

    ###########################################################################################
    '''))

    parser_PLOTTING._optionals.title = "OPTIONS"
    parser_PLOTTING.add_argument("-h", "--help", help="Show this help message for 'PLOTTING' and exit\n\n", action='help')

    parser_PLOTTING.add_argument("--heatmap", help="\nGenerate Heatmap Plot.\n\n", action="store_true")
    parser_PLOTTING.add_argument("--manhattan", help="\nGenerate Manhattan Plot.\n\n", action="store_true")






    ##### sub-parser for 'CONVERT'
    parser_CONVERT = subparsers.add_parser('CONVERT',
                                           help="Coverting results from other imputation platforms to the '*.ped' format.\n\n",
                                           formatter_class=argparse.RawTextHelpFormatter,
                                           add_help=False,
                                           description=textwrap.dedent('''\
    ###########################################################################################
        
        <CONVERT>
            
     

    ###########################################################################################
    '''))


    parser_CONVERT._optionals.title = "OPTIONS"
    parser_CONVERT.add_argument("-h", "--help", help="Show this help message for 'CONVERT' and exit\n\n", action='help')











    ### sub-parser for 'NOMENCLEANER'
    parser_NOMENCLEANER = subparsers.add_parser('NOMENCLEANER',
                                                help="Transforming HLA allele name field format.\n\n",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                add_help=False,
                                                description=textwrap.dedent('''\
    ###########################################################################################
        
        <NOMENCLENAER>

        explanation...

    ###########################################################################################
    '''))

    parser_NOMENCLEANER._optionals.title = "OPTIONS"
    parser_NOMENCLEANER.add_argument("-h", "--help", help="Show this help message for 'COATING' and exit\n\n", action='help')

    # Input (1) : *.ped file
    PED_TYPE = parser_NOMENCLEANER.add_mutually_exclusive_group(required=True)
    PED_TYPE.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", dest="ped",
                          default="Not_given")
    PED_TYPE.add_argument("-ped-Ggroup", help="\nHLA Type Data(G-group allele \"*.ped\" file).\n\n", dest="ped_G",
                          default="Not_given")
    PED_TYPE.add_argument("-ped-Pgroup", help="\nHLA Type Data(P-group allele \"*.ped\" file).\n\n", dest="ped_P",
                          default="Not_given")

    # Input (2) : *.iat file
    parser_NOMENCLEANER.add_argument("-iat", help="\nIntegrated Allele Table file(*.iat).\n\n", required=True)

    # Ouptut Prefix
    parser_NOMENCLEANER.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    # Output format
    output_digit_selection = parser_NOMENCLEANER.add_mutually_exclusive_group(required=True)
    output_digit_selection.add_argument("--1field", help="\nOutput ped file as '1-field' format.\n\n",
                                        action="store_true", dest="oneF")
    output_digit_selection.add_argument("--2field", help="\nOutput ped file as '2-field' format.\n\n",
                                        action="store_true", dest="twoF")
    output_digit_selection.add_argument("--3field", help="\nOutput ped file as '3-field' format.\n\n",
                                        action="store_true", dest="threeF")
    output_digit_selection.add_argument("--4field",
                                        help="\nOutput ped file as '4-field(Current Standard Names)' format.\n\n",
                                        action="store_true", dest="fourF")
    output_digit_selection.add_argument("--G-group", help="\nOutput ped file as 'G-group' format.\n\n",
                                        action="store_true")
    output_digit_selection.add_argument("--P-group", help="\nOutput ped file as 'P-group' format.\n\n",
                                        action="store_true")

    # Flag to remove HLA gene caption.
    parser_NOMENCLEANER.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*\").\n\n", action='store_true')








    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)



    ##### argument passing

    if args.subparser_name == "MAKEDICTIONARY":

        from MakeDictionary import MakeDictionary

        print("[%s]: Conducting 'MAKEDICTIONARY'." % (__file__))

        MakeDictionary(_HG=args.hg, _OUTPUT=args.o, _IMGT=args.imgt, _TYPE=args.type)




    elif args.subparser_name == "MAKEREFERENCE":

        print("[%s]: Conducting 'MAKEREFERENCE'." % (__file__))

        ##### Additional Argument processing #####

        ### Temporary Variables for Arguments.
        t_dict_AA = ""
        t_dict_AA_map = ""
        t_dict_SNPS = ""
        t_dict_SNPS_map = ""

        if (args.dict_AA != "Not_given" and args.dict_AA_map != "Not_given" and args.dict_SNPS != "Not_given" and args.dict_SNPS_map != "Not_given"):

            # When all HLA DICTIONARY information is given properly,

            t_dict_AA = args.dict_AA
            t_dict_AA_map = args.dict_AA_map
            t_dict_SNPS = args.dict_SNPS
            t_dict_SNPS_map = args.dict_SNPS_map


        elif (args.dict_AA == "Not_given" and args.dict_AA_map == "Not_given" and args.dict_SNPS == "Not_given" and args.dict_SNPS_map == "Not_given"):

            # No values are given to HLA DICTIONARY related options.

            if not args.previous_version:
                # Abort
                print("\n[Error]: None of HLA DICTIONARY files are given. Please check them all again.")
                print('{"-dict-AA", "-dict-AA-map", "-dict-SNPS", "-dict-SNPS-map"}\n')
                sys.exit()

            else:
                # In case of old version of MakeReference, then use default dataset.
                t_dict_AA = './data/MakeReference_old/HLA_DICTIONARY_AA.txt'
                t_dict_AA_map = './data/MakeReference_old/HLA_DICTIONARY_AA.map'
                t_dict_SNPS = './data/MakeReference_old/HLA_DICTIONARY_SNPS.txt'
                t_dict_SNPS_map = './data/MakeReference_old/HLA_DICTIONARY_SNPS.map'

        else:
            # Abort
            print("\n[Error]: Not all of HLA DICTIONARY files are given. Please check them all again.")
            print('{"-dict-AA", "-dict-AA-map", "-dict-SNPS", "-dict-SNPS-map"}\n')
            sys.exit()

        # (2018. 7. 16.)
        if args.previous_version:

            """
            If the option `--previous-version` is given, then the value of `-hg` option will be fixed to "18".
            I took this measure because previous version of the framework only covers hg 18.  
            
            (2018. 8. 1.)
            Argument "--only-markers" has been introduced. This argument also won't work with "--previous-version".
            
            """

            args.hg = "18"
            print("\n[Warning]: Only hg18 is available for previous version of MakeReference. Human Genome version will be set to hg 18 by force.\n")

            # (2018. 8. 1.)
            args.only_markers = False
            print("\n[Warning]: \"--only-markers\" argument is not available for previous version of MakeReference. It will be excluded automatically.\n")



        ### Implementing MakeReference.

        from MakeReference import MakeReference

        MakeReference(_INPUT_DATA=args.i, _HLA_ped=args.ped, _OUT=args.o, _dictionary_AA_map=t_dict_AA_map,
                      _dictionary_AA=t_dict_AA, _dictionary_SNPS_map=t_dict_SNPS_map, _dictionary_SNPS=t_dict_SNPS,
                      _previous_version=args.previous_version, _hg=args.hg, _only_markers=args.only_markers)




    elif args.subparser_name == "HLAANALYSIS":

        print(std_MAIN_PROCESS_NAME + "Conducting 'HLAANALYSIS'.")

        ##### Additional Argument processing #####

        f_ASSOCIATION_TEST = (args.logistic_regression or args.omnibus_test)
        f_META_ANALYSIS = args.meta_analysis

        if (int(f_ASSOCIATION_TEST) + int(f_META_ANALYSIS)) == 1:

            ### Loading "HLA_Analysis.py"
            import HLA_Analysis

            if f_ASSOCIATION_TEST:

                ### Implementing Association Test.
                HLA_Analysis.ASSOCIATION_TEST(_bfile=args.bfile, _input=args.input, _out=args.out, _covar=args.covar,
                                              _covar_names=args.covar_name, _phe=args.pheno, _phe_name=args.pheno_name,
                                              _condition=args.condition, _condition_list=args.condition_list,
                                              _ref_allele=args.reference_allele, _phased=args.phased, _threshold=args.rare_threshold,
                                              _lr=args.logistic_regression, _ob=args.omnibus_test)

            elif f_META_ANALYSIS:

                ### Implementing Meta-Analysis.
                HLA_Analysis.META_ANALYSIS(args.out, *args.rassoc)


        else:
            print(std_MAIN_PROCESS_NAME + "Error. Please conduct one analysis at a time.")
            sys.exit()



    elif args.subparser_name == "CONVERT":

        print("[%s]: Conducting 'CONVERT'." % (__file__))


    elif args.subparser_name == "NOMENCLEANER":

        print("[%s]: Conducting 'NOMENCLEANER'." % (__file__))

        ## Output Format Flags
        FILE_FORMAT = 1 if args.oneF else 2 if args.twoF else 3 if args.threeF else 4 if args.fourF else 5 if args.G_group else 6 if args.P_group else -1

        ## Which type of ped file given?
        _p_ped = -1
        _p_ped_descriptor = -1

        if args.ped != "Not_given":
            # Standard 4-field *.ped file given
            _p_ped = args.ped
            _p_ped_descriptor = 1
        elif args.ped_G != "Not_given":
            # G-group *.ped file given
            _p_ped = args.ped_G
            _p_ped_descriptor = 2
        elif args.ped_P != "Not_given":
            # P-group *.ped file given
            _p_ped = args.ped_P
            _p_ped_descriptor = 3
        else:
            # Assuming at least three of them given, there won't be the case which comes to here.
            print("\nArgument for input *.ped file has a problem. Please check it again.\n")
            sys.exit()

        ## If input ped file is given as G-group(or P-group), then user can't choose output format as same G-group(as same to P-group)
        if (_p_ped_descriptor == 2 and FILE_FORMAT == 5):
            print("\nYou just asked pointless transformation(Transformation G-group to G-group is pointless).")
            print("Skip this Transformation Request.\n")
            sys.exit()

            # (2018. 7. 10.) P to P or G to G 여도 NoDoubleColon인 경우는 할 수 있게 해줘야할듯.

        if (_p_ped_descriptor == 3 and FILE_FORMAT == 6):
            print("\nYou just asked pointless transformation(Transformation P-group to P-group is pointless).")
            print("Skip this Transformation Request.\n")
            sys.exit()


        ### Implementing NomenCleaner.

        from src.NomenCleaner import NomenCleaner

        NomenCleaner(_p_ped, _p_ped_descriptor, args.iat, args.o, FILE_FORMAT, _f_NoCaption=args.NoCaption)
