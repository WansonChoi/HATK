#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap




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
    parser_HLAANALYSIS = subparsers.add_parser('HLAANALYSIS',
                                               help="HLA_analysis help\n\n",
                                               formatter_class=argparse.RawTextHelpFormatter,
                                               add_help=False,
                                               description=textwrap.dedent('''\
    ###########################################################################################
        
        <HLA-Analysis>
     
            1. 만들게 될 HLAStudy 클래스의 instance가 한 개일때
                - 그냥 study instance의 summary statistics만 곱게 리턴(Logistic_Regression, OmnibusTest, BinaryTest)
                - BinaryMarker_manhattan
                - AA_manhattan
                - HEATMAP
                - Conditional Analysis
            
            2. 만들게 될 HLAStudy 클래스의 instance가 두 개일때(다뤄야할 study가 두 개 이하라고 가정)
                - Meta-Analysis(얘는 어쩌면 진짜 파일로 받을 필요도 있겠다..)

    ###########################################################################################
    '''))

    parser_HLAANALYSIS._optionals.title = "OPTIONS"

    parser_HLAANALYSIS.add_argument("-h", "--help", help="\nShow this help message for 'HLAANALYSIS' and exit\n\n", action='help')

    # Basic necessary information for STUDY1 instance of HLAStudy class.
    parser_HLAANALYSIS.add_argument("-s1t", "--study1-tag", help="\nSpecify study1 tag(name) for analysis\n\n",
                                    nargs=1, required=True)
    parser_HLAANALYSIS.add_argument("-s1p", "--study1-phe", help="\nSpecify study1 phenotype name(column name in '*.phe' file) for analysis\n\n",
                                    nargs=1, required=True)
    parser_HLAANALYSIS.add_argument("-s1f", "--study1-file", help="\nSpecify study1 file prefix for analysis\n\n",
                                    nargs=1, required=True)
    parser_HLAANALYSIS.add_argument("-s1c", "--study1-condition", help="\nSpecify condition for study1\n\n\n\n",
                                    nargs=1) # 기본적으로는 파일로 받게, 이 외에 2-3개 정도는 그냥 csv string형태로도 받을 수 있게

    parser_HLAANALYSIS.add_argument("-s2t", "--study2-tag", help="\nSpecify study2 tag(name) for analysis\n\n",
                                    nargs=1)
    parser_HLAANALYSIS.add_argument("-s2p", "--study2-phe", help="\nSpecify study2 phenotype name(column name in '*.phe' file) for analysis\n\n",
                                    nargs=1)
    parser_HLAANALYSIS.add_argument("-s2f", "--study2-file", help="\nSpecify study2 file prefix for analysis\n\n",
                                    nargs=1)
    parser_HLAANALYSIS.add_argument("-s2c", "--study2-condition", help="\nSpecify condition for study2\n\n",
                                    nargs=1) # 기본적으로는 파일로 받게, 이 외에 2-3개 정도는 그냥 csv string형태로도 받을 수 있게


    group1 = parser_HLAANALYSIS.add_argument_group(title='group1', description = "Options for single study analysis.")
    group1.add_argument("-lr", "--logistic-regression", help="\nConducting logistic regression(Plink v1.07)\n\n",
                        action="store_true")
    group1.add_argument("-ot", "--omnibus-test", help="\nConducting OmnibusTest.\n\n",
                        action="store_true")
    group1.add_argument("-bt", "--binary-test", help="\nConducting BinaryTest(ANOVA)\n\n",
                        action="store_true")

    # 플로팅(AA_manhattan, Binary_manhattan, HEATMAP) 얘네들 어디다 집어넣을지 고민해볼것

    group2 = parser_HLAANALYSIS.add_argument_group(title='group2', description = "Options for two studies analysis.")
    group2.add_argument("-m", "--meta-analysis", help="\nConducting Meta-Analysis(Inverse-Variance Weighting) for two studies.\n\n",
                        action='store_true')



    # 아예 logistic regression, binary test, omnibus test, aa_manhattan, binary_manhattan, heatmap 이걸 단위로 group으로 묶는게 나을듯(하도 딸려 부가적으로 필요한 argument가 많아서)
    # 생각해보니 logistic_regression은 디폴트로 해줘야겠구나.. 굳이 옵션 안만들어도 되것다 가 아니라 이경우에는 logistic regression '만' 수행하게 플래그들 케이스분류해야것다.







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

    # for Test
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPPER_TEST/imgt370/WRAPER_TEST", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPPER_TEST/imgt3300_hg19/WRAPER_TEST", "-imgt", "3300"])
    # args = parser.parse_args(["-hg", "18", "-o", "WRAPPER_TEST/imgt3300_hg18/WRAPER_TEST", "-imgt", "3300"])
    # args = parser.parse_args(["-hg", "38", "-o", "WRAPPER_TEST/imgt3300_hg38/WRAPER_TEST", "-imgt", "3300"])

    # args = parser.parse_args(["-hg", "19", "-o", "WRAPER_TEST", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPER_TEST", "-imgt", "3300"])

    print(args)

    # argument passing



    if args.subparser_name == "MAKEDICTIONARY":

        from MakeDictionary import MakeDictionary

        print("[%s]: Conducting 'Makedictionary' module." % (__file__))

        MakeDictionary(_HG=args.hg, _OUTPUT=args.o, _IMGT=args.imgt, _TYPE=args.type)




    elif args.subparser_name == "MAKEREFERENCE":

        print("[%s]: Conducting 'MAKEREFERENCE' module." % (__file__))

        ##### Additional Argument processing

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

        print("[%s]: Conducting 'HLAANALYSIS' module." % (__file__))


    elif args.subparser_name == "CONVERT":

        print("[%s]: Conducting 'CONVERT' module." % (__file__))


    elif args.subparser_name == "COATING":

        print("[%s]: Conducting 'COATING' module." % (__file__))
