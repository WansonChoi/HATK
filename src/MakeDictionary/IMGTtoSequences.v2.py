# -*- coding: utf-8 -*-

import os, sys, subprocess
import pandas as pd
import argparse, textwrap
import re

def IMGTtoSequences(_out, _hla, _type, _hg_Table, _imgt,
                    _nuc, _gen = "Not_given", _prot = "Not_given"):


    ########## < Core Variables > ##########


    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting IMGTtoSequences.")


    ### General
    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}




    ### (2018/1/17) HLA position information(exon, intron, etc.)
    HLA_INTEGRATED_POSITIONS_hg = pd.read_table(_hg_Table, sep='\t', header=None,
                                                usecols = [0,1,2,3],
                                                names=['HLA_name', 's_pos', 'e_pos', 'Class'],
                                                index_col=[0, 3])
    # print(HLA_INTEGRATED_POSITIONS_hg.head())




    ########## < Argument Checking. > ##########

    # Preparing intermediate paths.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')

    INTERMEDIATE_PATH = os.path.dirname(_out)

    if not os.path.exists(INTERMEDIATE_PATH):
        os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))







if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        IMGTtoSequences.v2.py

        : Processing(Parsing HLA sequence information distributed by IMGT-HLA.
        
        Renewed in 2018.09.10


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("-HLA", help="\nHLA gene name which you will process.\n\n", required=True, metavar='HLA',
                        choices=["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"])

    parser.add_argument("--type", "-t", help="\nSequence type to deal with(Amino Acids[AA] or SNPs[SNPS]\n\n", required=True, choices=["AA", "SNPS"], metavar='TYPE')
    parser.add_argument("--hg-table", "-hg", help="\nHLA gene position information table.\n\n", required=True)
    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="imgt_version", required=True)

    parser.add_argument("-nuc", help="\nInput *_nuc.txt file.\n\n", required=True)
    parser.add_argument("-gen", help="\nInput *_gen.txt file.\n\n", default="Not_given")
    parser.add_argument("-prot", help="\nInput *_prot.txt file.\n\n", default="Not_given")



    ##### < for Test > #####

    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "A",
    #                           "--type", "SNPS",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt"
    #                           ])



    ##### < for Publish > #####
    args = parser.parse_args()

    print(args)


    # main function execution
    IMGTtoSequences(_out = args.o, _hla=args.HLA, _type=args.type, _hg_Table=args.hg_table, _imgt=args.imgt,
                    _nuc=args.nuc, _gen=args.gen, _prot=args.prot)
