#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))

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






    ### sub-parser for 'MakeDictionary'
    parser_MAKEDICTIONARY = subparsers.add_parser('MAKEDICTIONARY',
                                                  help='MakeDictionary help\n\n',
                                                  formatter_class=argparse.RawTextHelpFormatter,
                                                  add_help=False,
                                                  description=textwrap.dedent('''\
    ###########################################################################################
        
        <MakeDictionary>
    
        explanation...
    
    ###########################################################################################
    '''))
    # parser_makedictionary = subparsers.add_parser('MKDICT', help='makedictionary help', description="This is a desription for makedictionary")

    parser_MAKEDICTIONARY._optionals.title = "OPTIONS"
    parser_MAKEDICTIONARY.add_argument("-h", "--help", help="\nShow this help message for 'MAKEDICTIONARY' and exit\n\n", action='help')

    parser_MAKEDICTIONARY.add_argument("-hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                                       metavar="HG", required=True)
    # 테스트해보니 "-hg" 이거 sub-parser내부에서는 중복으로 선언해도 상관없는듯.
    # cf) namespace 이름 설정해주는건 dest파라미터였음.

    parser_MAKEDICTIONARY.add_argument("-o", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser_MAKEDICTIONARY.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                                       required=True)

    parser_MAKEDICTIONARY.add_argument("--type", help="\nSNPS or AA or Both\n\n", choices=["AA", "SNPS", "BOTH"],
                                       metavar="TYPE", default="BOTH")








    ### sub-parser for 'MakeReference'
    parser_MAKEREFERENCE = subparsers.add_parser('MAKEREFERENCE',
                                                  help='MakeReference help\n\n',
                                                  formatter_class=argparse.RawTextHelpFormatter,
                                                  add_help=False,
                                                  description=textwrap.dedent('''\
    #################################################################################################
        
        MakeReference.py
        
        This script helps prepare a reference dataset for HLA imputation

        Usage(1)
        : python3 MakeReference.py --previous-version -i ./data/MakeReference_old/HAPMAP_CEU 
            -ped ./data/MakeReference_old/HAPMAP_CEU_HLA.ped -hg 18 -o ./Trial_HAPMAP_CEU
        
        Usage(2)
        : python3 MakeReference.py -i ./data/MakeReference/HAPMAP_CEU 
            -ped ./data/MakeReference/HAPMAP_CEU_HLA.4field.ped -hg 18 -o ./Trial_HAPMAP_CEU
            -dict-AA ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt 
            -dict-AA-map ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map 
            -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt 
            -dict-SNPS-map ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map
                
        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
    '''))

    parser_MAKEREFERENCE._optionals.title = "OPTIONS"

    parser_MAKEREFERENCE.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser_MAKEREFERENCE.add_argument("-i", help="\nInput Data file(.bed/.bim/.fam)\n\n", required=True)
    parser_MAKEREFERENCE.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n", required=True)
    parser_MAKEREFERENCE.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
    parser_MAKEREFERENCE.add_argument("-o", help="\nOutput file prefix\n\n")

    parser_MAKEREFERENCE.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as old version.\n\n", action='store_true')
    parser_MAKEREFERENCE.add_argument("--only-markers", help="\nGenerate a panel for only markers, not entire reference panel.\n\n", action='store_true') # (2018. 8. 1.) Introduced.

    hla_dict = parser_MAKEREFERENCE.add_argument_group(title='HLA_DICTIONARY',
                                         description='- Arguments to specify HLA_DICTIONARY Information to New version of MakeReference(2.0)\n'
                                                     '- If you\'re going to use previous version of MakeReference, then Don\'t care about these options.')

    hla_dict.add_argument("-dict-AA", help="\nInput HLA Dictionary file for AA Information.\n\n", default="Not_given")
    hla_dict.add_argument("-dict-AA-map", help="\nInput HLA Dictionary .map file for AA Information.\n\n", default="Not_given")
    hla_dict.add_argument("-dict-SNPS", help="\nInput HLA Dictionary file for SNPS Information.\n\n", default="Not_given")
    hla_dict.add_argument("-dict-SNPS-map", help="\nInput HLA Dictionary .map file for SNPS Information\n\n", default="Not_given")




    # (2018. 6. 5) sub-parser for 'SNP2HLA' also will be introduced here, soon. => Prof. Han didn't want it. Canceled.




    ### sub-parser for 'HLA-Analysis'
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


    group_INPUT = parser_HLA_ANALYSIS.add_argument_group(title="INPUT" ,description='')
    # description = '- The ways to give input files. (ex. {".ped",".map"} or {".bed",".bim",."fam"})'

    PEDorBED = group_INPUT.add_mutually_exclusive_group(required=True)
    PEDorBED.add_argument("--file", "-f", help="\nInput file Prefixes of \".ped\" and \".map\" files (Plink v1.07).\n\n",
                          nargs = "+")
    PEDorBED.add_argument("--bfile", "-bf", help="\nInput file Prefixes of \".bed\", \".fam\" and \".bim\" files (Plink v1.07).\n\n",
                          nargs = "+")

    group_INPUT.add_argument("--covar", help="\nSpecify .covar file (Plink v1.07).\n\n")
    group_INPUT.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use.(Plink v1.07).\n\n")

    group_INPUT.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.07).\n\n")
    group_INPUT.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use.(Plink v1.07).\n\n")

    CondVars = group_INPUT.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify Marker name(s) to condition by comma-separated(\",\") (Plink v1.07).\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify Marker name(s) to condition as a file (Plink v1.07).\n\n")

    group_INPUT.add_argument("--reference-allele", help="\nSpecify Reference Allele file (Plink v1.07).\n\n")


    group_INPUT.add_argument("--out", "-o", help="\nSpecify Output Prefix (Plink v1.07).\n\n", required=True)




    # (2018. 8. 2.) Later, Introduce the way that argument can take ped and map files seperately.

    group_AssociationTest = parser_HLA_ANALYSIS.add_argument_group(title='ASSOCIATION TEST', description='')

    group_AssociationTest.add_argument("--logistic-regression", "-lr", help="\nConducting Logistic Regression(Plink v1.07)\n\n",
                                       action="store_true")
    group_AssociationTest.add_argument("--omnibus-test", "-ob", help="\nConducting OmnibusTest.\n\n",
                                       action="store_true")
    group_AssociationTest.add_argument("--binary-test", "-bt", help="\nConducting BinaryTest(ANOVA).\n\n",
                                       action="store_true")

    # Later, by using exclusive-nor, Don't allow using "Association Test", "Meta-analysis" and "plotting"(maybe this can be fine.).

    group_MetaAnalysis = parser_HLA_ANALYSIS.add_argument_group(title="META ANALYSIS", description='')

    group_MetaAnalysis.add_argument("--meta-analysis", "-meta", help="\nConducting Meta Analysis(Plink v1.07)\n\n",
                                    action="store_true")

    group_Plotting = parser_HLA_ANALYSIS.add_argument_group(title="PLOTTING", description='')

    group_Plotting.add_argument("--heatmap", help="\nGenerate Heatmap Plot.\n\n",
                                    action="store_true")
    group_Plotting.add_argument("--manhattan", help="\nGenerate Manhattan Plot.\n\n",
                                    action="store_true")

    # (2018. 8. 2.) deprecated.
    # # Basic necessary information for STUDY1 instance of HLAStudy class.
    # parser_HLA_ANALYSIS.add_argument("-s1t", "--study1-tag", help="\nSpecify study1 tag(name) for analysis\n\n",
    #                                  nargs=1, required=True)
    # parser_HLA_ANALYSIS.add_argument("-s1p", "--study1-phe", help="\nSpecify study1 phenotype name(column name in '*.phe' file) for analysis\n\n",
    #                                  nargs=1, required=True)
    # parser_HLA_ANALYSIS.add_argument("-s1f", "--study1-file", help="\nSpecify study1 file prefix for analysis\n\n",
    #                                  nargs=1, required=True)
    # parser_HLA_ANALYSIS.add_argument("-s1c", "--study1-condition", help="\nSpecify condition for study1\n\n\n\n",
    #                                  nargs=1) # 기본적으로는 파일로 받게, 이 외에 2-3개 정도는 그냥 csv string형태로도 받을 수 있게
    #
    # parser_HLA_ANALYSIS.add_argument("-s2t", "--study2-tag", help="\nSpecify study2 tag(name) for analysis\n\n",
    #                                  nargs=1)
    # parser_HLA_ANALYSIS.add_argument("-s2p", "--study2-phe", help="\nSpecify study2 phenotype name(column name in '*.phe' file) for analysis\n\n",
    #                                  nargs=1)
    # parser_HLA_ANALYSIS.add_argument("-s2f", "--study2-file", help="\nSpecify study2 file prefix for analysis\n\n",
    #                                  nargs=1)
    # parser_HLA_ANALYSIS.add_argument("-s2c", "--study2-condition", help="\nSpecify condition for study2\n\n",
    #                                  nargs=1) # 기본적으로는 파일로 받게, 이 외에 2-3개 정도는 그냥 csv string형태로도 받을 수 있게








    ### sub-parser for 'CONVERT'
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











    ### sub-parser for 'COATING'
    parser_COATING = subparsers.add_parser('COATING',
                                           help="Transforming HLA allele naming system(2-field vs. 4-field).\n\n",
                                           formatter_class=argparse.RawTextHelpFormatter,
                                           add_help=False,
                                           description=textwrap.dedent('''\
    ###########################################################################################
        
        <COATING>

        explanation...

    ###########################################################################################
    '''))

    parser_COATING._optionals.title = "OPTIONS"
    parser_COATING.add_argument("-h", "--help", help="Show this help message for 'COATING' and exit\n\n", action='help')





    # for Publish
    args = parser.parse_args()
    print(args)

    # for Test
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPPER_TEST/imgt370/WRAPER_TEST", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPPER_TEST/imgt3300_hg19/WRAPER_TEST", "-imgt", "3300"])
    # args = parser.parse_args(["-hg", "18", "-o", "WRAPPER_TEST/imgt3300_hg18/WRAPER_TEST", "-imgt", "3300"])
    # args = parser.parse_args(["-hg", "38", "-o", "WRAPPER_TEST/imgt3300_hg38/WRAPER_TEST", "-imgt", "3300"])

    # args = parser.parse_args(["-hg", "19", "-o", "WRAPER_TEST", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPER_TEST", "-imgt", "3300"])


    # argument passing



    if args.subparser_name == "MAKEDICTIONARY":

        from MakeDictionary import MakeDictionary

        print("[%s]: Conducting 'Makedictionary'." % (__file__))

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

        ### Temporary Variables for Arguments.
        t_input = None


        ### Checking Input file type(.ped vs. .bed)

        if bool(args.bfile):

            # When input is given as .bed, .bim, .fam.

            if len(args.bfile) == 1:
                print(std_MAIN_PROCESS_NAME + "Single input set.")
            elif len(args.bfile) > 1:
                print(std_MAIN_PROCESS_NAME + "Multiple Inputs sets.")
                #(2018. 8. 2.) Merged or not. Choose to introduce later.
                # Maybe multiple input will be excluded.

            for item in args.bfile:

                if not (os.path.exists(item + ".bed") and os.path.exists(item + ".bim") and os.path.exists(item + ".fam")):
                    print(std_MAIN_PROCESS_NAME + "Error. Not all of input files(\".bed\", \".bim\", \".fam\") exist.")
                    sys.exit()

            args.bfile = args.bfile.pop()


        elif bool(args.file):

            # When input is given as .ped and .map.

            if len(args.file) == 1:
                print(std_MAIN_PROCESS_NAME + "Single input set.")
            elif len(args.file) > 1:
                print(std_MAIN_PROCESS_NAME + "Multiple Inputs sets.")
                #(2018. 8. 2.) Merged or not. Choose to introduce later.

            for item in args.file:

                if not (os.path.exists(item + ".ped") and os.path.exists(item + ".map")):
                    print(std_MAIN_PROCESS_NAME + "Error. Not all of input files(\".ped\", \".map\") exist.")
                    sys.exit()

            args.file = args.file.pop()


        ### Implementing "HLA_Analysis()"
        from HLA_Analysis import HLA_Analysis

        HLA_Analysis(_bfile=args.bfile, _file=args.file, _out=args.out,
                     _covar=args.covar, _covar_names=args.covar_name,
                     _phe=args.pheno, _phe_name=args.pheno_name,
                     _condition=args.condition, _condition_list=args.condition_list, _ref_allele=args.reference_allele,
                     _lr=args.logistic_regression, _ob=args.omnibus_test, _bt=args.binary_test, _meta=args.meta_analysis,
                     _heatmap=args.heatmap, _manhattan=args.manhattan)



    elif args.subparser_name == "CONVERT":

        print("[%s]: Conducting 'CONVERT'." % (__file__))


    elif args.subparser_name == "COATING":

        print("[%s]: Conducting 'COATING'." % (__file__))
