# -*- coding: utf-8 -*-

import os, sys, re
from pandas import concat
import argparse, textwrap
from glob import glob


def MakeDictionary(_HG, _OUTPUT, _IMGT, _TYPE="BOTH"):

    """
    """

    ########## < Core Variables > ##########

    ### Module name.
    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting MakeDictionary.")


    ### General
    raw_HLA_names = ["A", "C", "B", "DRB", "DQA", "DQB", "DPA", "DPB"]
    HLA_names = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]


    ### Dictionaries for Raw files.
    TARGET_prot_files = {}
    # TARGET_nuc_files = {}
    TARGET_gen_files = {}


    ### Variables for Paths.
    p_data = "./data/MakeDictionary"
    p_src = "./src/MakeDictionary"


    ### OUTPUT prefix

    _OUTPUT = _OUTPUT if not _OUTPUT.endswith('/') else _OUTPUT.rstrip('/')

    INTERMEDIATE_PATH = os.path.dirname(_OUTPUT)


    if not os.path.exists(INTERMEDIATE_PATH):
        os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))




    ########## < Dependency Checking > ##########

    """
    1. Necessary Python source files.
    2. Necessary Static Data files
    3. Checking each files of "*_nuc.txt", "*_prot.txt", "*_gen.txt".
    """


    ##### < Necessary Python source files > #####

    # IMGTtoSequences.py
    if os.path.exists(os.path.join(p_src, "IMGTtoSequences.py")):
        from src.MakeDictionary.IMGTtoSequences import IMGTtoSequences
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "\"IMGTtoSequences.py\" doesn't exist!")
        sys.exit()

    # MakeMap.py
    if os.path.exists(os.path.join(p_src, "MakeMap.py")):
        from src.MakeDictionary.MakeMap import MakeMap
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "\"MakeMap.py\" doesn't exist!")
        sys.exit()



    ##### < Necessary Static Data files > #####

    HLA_INTEGRATED_POSITIONS_filename = os.path.join(p_data, "HLA_INTEGRATED_POSITIONS_hg{0}.txt".format(_HG))
    IMGTHLA_directory = os.path.join(p_data, "IMGTHLA{0}".format(_IMGT))

    print(HLA_INTEGRATED_POSITIONS_filename)
    print(IMGTHLA_directory)

    if not os.path.exists(HLA_INTEGRATED_POSITIONS_filename):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}\" not found!".format(HLA_INTEGRATED_POSITIONS_filename))
        sys.exit()

    if not os.path.isdir(IMGTHLA_directory):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}\" not found!".format(IMGTHLA_directory))
        sys.exit()



    ##### < "*_nuc.txt", "*_prot.txt", "*_gen.txt" > #####

    files_in_alignments = glob(IMGTHLA_directory+'/alignments/*')

    for i in range(0, len(HLA_names)):

        ### <*_prot.txt, *_nuc.txt>
        found_prot = list(filter(re.compile("/"+ raw_HLA_names[i]+"_prot.txt").search, files_in_alignments))
        # found_nuc = list(filter(re.compile("/" + raw_HLA_names[i] + "_nuc.txt").search, files_in_alignments))


        if len(found_prot) == 1:
            # When found exactly.
            TARGET_prot_files[HLA_names[i]] = found_prot.pop()
            # TARGET_nuc_files[HLA_names[i]] = found_nuc.pop()
        else:
            # If not, then try same job to the names with number "1".
            found_prot = list(filter(re.compile("/" + HLA_names[i] + "_prot.txt").search, files_in_alignments))
            # found_nuc = list(filter(re.compile("/" + HLA_names[i] + "_nuc.txt").search, files_in_alignments))

            if len(found_prot) == 1:
                # When found exactly in 2nd trial.
                TARGET_prot_files[HLA_names[i]] = found_prot.pop()
                # TARGET_nuc_files[HLA_names[i]] = found_nuc.pop()
            else:
                # Not found finally.
                print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}_prot.txt\" file can't be found.".format(HLA_names[i]))
                sys.exit()


        ### <*_gen.txt>

        found_gen = list(filter(re.compile("/" + raw_HLA_names[i] + "_gen.txt").search, files_in_alignments))

        if len(found_gen) == 1:
            # When found exactly.
            TARGET_gen_files[HLA_names[i]] = found_gen.pop()
        else:
            # If not, then try same job to the names with number "1".
            found_gen = list(filter(re.compile("/" + HLA_names[i] + "_gen.txt").search, files_in_alignments))

            if len(found_gen) == 1:
                # When found exactly in 2nd trial.
                TARGET_gen_files[HLA_names[i]] = found_gen.pop()
            else:
                # Not found finally.
                print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}_gen.txt\" file can't be found.".format(HLA_names[i]))
                sys.exit()

    print("\nTarget prot and gen files : \n")
    print(TARGET_prot_files)
    print(TARGET_gen_files)





    ########## < Control Flags > ##########

    MAKING_AA_DICTIONARY = 1
    MAKING_AA_MAPFILE = 1
    MAKING_SNPS_DICTIONARY = 1
    MAKING_SNPS_MAPFILE = 1
    CLEAN_UP = 0




    if MAKING_AA_DICTIONARY:

        ########## < 1. Making AA Dictionary. > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Making AA dictionary file.")

        HLA_DICTIONARY_AA = concat([IMGTtoSequences(_inputfile=TARGET_prot_files[HLA_names[i]],
                                                    _OUTPUT='_'.join([_OUTPUT, HLA_names[i]]), _TYPE="AA",
                                                    _HLA=HLA_names[i], _HG=_HG, _return_as_dataframe=True,
                                                    _p_HLA_INTEGRATED_POSITIONS=HLA_INTEGRATED_POSITIONS_filename)
                                    for i in range(0, len(HLA_names))], axis=0)


        HLA_DICTIONARY_AA.to_csv(os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}.txt'.format(_HG, _IMGT)),
                                 sep='\t', header=False, index=True)





    if MAKING_AA_MAPFILE:

        ########## < 2. Making AA Map file. > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Making AA Map file.")

        HLA_DICTIONARY_AA_map = concat([MakeMap(_inputfile=_OUTPUT+"_{0}.AA.forMAP.txt".format(HLA_names[i]),
                                                _TYPE="AA", _HLA=HLA_names[i], _return_as_dataframe=True)
                                        for i in range(0, len(HLA_names))], axis=0)


        HLA_DICTIONARY_AA_map.to_csv(os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
                                     sep='\t', header=False, index=False)




    if MAKING_SNPS_DICTIONARY:

        ########## < 3. Making SNPS Dictionary. > ##########

        print(std_MAIN_PROCESS_NAME + "[3] Making SNPS dictionary file.")


        HLA_DICTIONARY_SNPS = concat([IMGTtoSequences(_inputfile=TARGET_gen_files[HLA_names[i]],
                                                      _OUTPUT='_'.join([_OUTPUT, HLA_names[i]]), _TYPE="SNPS",
                                                      _HLA=HLA_names[i], _HG=_HG, _return_as_dataframe=True,
                                                      _p_HLA_INTEGRATED_POSITIONS=HLA_INTEGRATED_POSITIONS_filename)
                                      for i in range(0, len(HLA_names))], axis=0)



        HLA_DICTIONARY_SNPS.to_csv(os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.txt'.format(_HG, _IMGT)),
                                   sep='\t', header=False, index=True)





    if MAKING_SNPS_MAPFILE:

        ########## < 4. Making SNPS Map file. > ##########

        print(std_MAIN_PROCESS_NAME + "[4] Making SNPS Map file.")

        HLA_DICTIONARY_SNPS_map = concat([MakeMap(_inputfile=_OUTPUT+"_{0}.SNPS.forMAP.txt".format(HLA_names[i]),
                                                  _TYPE="SNPS", _HLA=HLA_names[i], _return_as_dataframe=True)
                                          for i in range(0, len(HLA_names))], axis=0)


        HLA_DICTIONARY_SNPS_map.to_csv(os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
                                       sep='\t', header=False, index=False)





    if CLEAN_UP:

        ########## < 5. Removing unnecessary files. > ##########

        print(std_MAIN_PROCESS_NAME + "[5] Removing unnecessary files")



    print("\n")
    print(std_MAIN_PROCESS_NAME + "Making Dictionary is done.\n")

    return 0




if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################
    
        MakeDictionary.py

        
    
    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
    parser.add_argument("-hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"], metavar="hg", required=True)
    parser.add_argument("-o", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version", required=True)
    parser.add_argument("--type", "-t", help="\nSNPS or AA or Both\n\n", choices=["AA", "SNPS", "BOTH"], metavar="Type", default="BOTH")


    ##### < for Test > #####

    # args = parser.parse_args(["-hg", "18", "-o", "TEST20180812/imgt370_hg18/HLA_HLA_HLA", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPPER_TEST/imgt370/WRAPER_TEST", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPPER_TEST/imgt3300_hg19/WRAPER_TEST", "-imgt", "3300"])
    # args = parser.parse_args(["-hg", "18", "-o", "WRAPPER_TEST/imgt3300_hg18/WRAPER_TEST", "-imgt", "3300"])
    # args = parser.parse_args(["-hg", "38", "-o", "WRAPPER_TEST/imgt3300_hg38/WRAPER_TEST", "-imgt", "3300"])

    # args = parser.parse_args(["-hg", "19", "-o", "WRAPER_TEST", "-imgt", "370"])
    # args = parser.parse_args(["-hg", "19", "-o", "WRAPER_TEST", "-imgt", "3300"])


    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    ##### < Main function Execution. > #####
    MakeDictionary(_HG = args.hg, _OUTPUT=args.o, _IMGT=args.imgt)

