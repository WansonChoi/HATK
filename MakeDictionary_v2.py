# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd
import argparse, textwrap
from glob import glob


def MakeDictionary(_HG, _OUTPUT, _IMGT, _TYPE="BOTH", _no_Indel=False):

    """
    """

    ########## < Core Variables > ##########

    ### Module name.
    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting MakeDictionary.")


    # (2018. 8. 28.) This order must be kept.
    raw_HLA_names = ["A", "B", "C", "DPA", "DPB", "DQA", "DQB", "DRB"]
    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


    ### Dictionaries for Raw files.
    TARGET_prot_files = {}
    TARGET_nuc_files = {}
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
    if os.path.exists(os.path.join(p_src, "IMGTtoSequences_v2.py")):
        from src.MakeDictionary.IMGTtoSequences_v2 import IMGTtoSequences
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "\"IMGTtoSequences.v2.py\" doesn't exist!")
        sys.exit()

    # # MakeMap.py
    # if os.path.exists(os.path.join(p_src, "MakeMap.py")):
    #     from src.MakeDictionary.MakeMap import MakeMap
    # else:
    #     print(std_ERROR_MAIN_PROCESS_NAME + "\"MakeMap.py\" doesn't exist!")
    #     sys.exit()


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

    files_in_alignments = glob(IMGTHLA_directory + '/alignments/*')

    for i in range(0, len(HLA_names)):

        ### <*_prot.txt, *_nuc.txt>
        found_prot = list(filter(re.compile("/" + raw_HLA_names[i] + "_prot.txt").search, files_in_alignments))
        found_nuc = list(filter(re.compile("/" + raw_HLA_names[i] + "_nuc.txt").search, files_in_alignments))

        if len(found_prot) == 1:
            # When found exactly.
            TARGET_prot_files[HLA_names[i]] = found_prot.pop()
            TARGET_nuc_files[HLA_names[i]] = found_nuc.pop()

        else:
            # If not, then try same job to the names with number "1".
            found_prot = list(filter(re.compile("/" + HLA_names[i] + "_prot.txt").search, files_in_alignments))
            found_nuc = list(filter(re.compile("/" + HLA_names[i] + "_nuc.txt").search, files_in_alignments))

            if len(found_prot) == 1:
                # When found exactly in 2nd trial.
                TARGET_prot_files[HLA_names[i]] = found_prot.pop()
                TARGET_nuc_files[HLA_names[i]] = found_nuc.pop()
            else:
                # Not found finally.
                print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}_prot.txt\" or \"{0}_nuc.txt\" file can't be found.".format(HLA_names[i]))
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

    MAKING_DICTIONARY = 1
    MAKING_MAPFILE = 0
    MAKING_SNPS_DICTIONARY = 0
    # MAKING_SNPS_MAPFILE = 0
    CLEAN_UP = 0


    if MAKING_DICTIONARY:

        ########## < 1. Making AA Dictionary. > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Making dictionary file.")

        l_df_Seqs_AA = []
        l_df_forMAP_AA = []

        l_df_Seqs_SNPS = []
        l_df_forMAP_SNPS = []


        for i in range(0, len(HLA_names)):

            t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA = IMGTtoSequences(
                _OUTPUT, HLA_names[i], HLA_INTEGRATED_POSITIONS_filename, _IMGT,
                TARGET_nuc_files[HLA_names[i]], TARGET_gen_files[HLA_names[i]], TARGET_prot_files[HLA_names[i]])

            # Dictionary AA.
            l_df_Seqs_AA.append(t_df_Seqs_AA)
            l_df_forMAP_AA.append(MakeMap(HLA_names[i], "AA", t_df_forMAP_AA))

            # Dictionary SNPS.
            l_df_Seqs_SNPS.append(t_df_Seqs_SNPS)
            l_df_forMAP_SNPS.append(MakeMap(HLA_names[i], "SNPS", t_df_forMAP_SNPS))



        HLA_DICTIONARY_AA = pd.concat(l_df_Seqs_AA, axis=0)
        HLA_DICTIONARY_AA_map = pd.concat(l_df_forMAP_AA, axis=0)

        HLA_DICTIONARY_SNPS = pd.concat(l_df_Seqs_SNPS, axis=0)
        HLA_DICTIONARY_SNPS_map = pd.concat(l_df_forMAP_SNPS, axis=0)


        # Exporting AA dictionary.
        HLA_DICTIONARY_AA.to_csv(
            os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}.txt'.format(_HG, _IMGT)),
            sep='\t', header=False, index=True)

        HLA_DICTIONARY_AA_map.to_csv(
            os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
            sep='\t', header=False, index=False
        )

        # Exporting SNPS dictionary.
        HLA_DICTIONARY_SNPS.to_csv(
            os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.txt'.format(_HG, _IMGT)),
            sep='\t', header=False, index=True)

        HLA_DICTIONARY_SNPS_map.to_csv(
            os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
            sep='\t', header=False, index=False
        )




    # if MAKING_MAPFILE:
    #
    #     ########## < 2. Making AA Map file. > ##########
    #
    #     print(std_MAIN_PROCESS_NAME + "[2] Making AA Map file.")
    #
    #     HLA_DICTIONARY_AA_map = concat([MakeMap(_inputfile=_OUTPUT + "_{0}.AA.forMAP.txt".format(HLA_names[i]),
    #                                             _TYPE="AA", _HLA=HLA_names[i], _return_as_dataframe=True)
    #                                     for i in range(0, len(HLA_names))], axis=0)
    #
    #     HLA_DICTIONARY_AA_map.to_csv(
    #         os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
    #         sep='\t', header=False, index=False)


    # if MAKING_SNPS_DICTIONARY:
    #
    #     ########## < 3. Making SNPS Dictionary. > ##########
    #
    #     print(std_MAIN_PROCESS_NAME + "[3] Making SNPS dictionary file.")
    #
    #     l_df_Seqs = []
    #     l_df_forMAPs = []
    #
    #     for i in range(0, len(HLA_names)):
    #
    #         t_df_Seqs, t_df_forMAPs = IMGTtoSequences(_out = _OUTPUT, _hla=HLA_names[i],
    #                                                 _hg_Table=HLA_INTEGRATED_POSITIONS_filename,
    #                                                 _imgt=_IMGT, _type="SNPS",
    #                                                 _gen=TARGET_gen_files[HLA_names[i]])
    #
    #         l_df_Seqs.append(t_df_Seqs)
    #         # l_df_forMAPs.append(MakeMap(t_df_forMAPs, _type="SNPS", _hla=HLA_names[i]))
    #
    #
    #     HLA_DICTIONARY_SNPS = pd.concat(l_df_Seqs, axis=0)
    #     # HLA_DICTIONARY_SNPS_map = pd.concat(l_df_forMAPs, axis=0)
    #
    #
    #     HLA_DICTIONARY_SNPS.to_csv(
    #         os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.txt'.format(_HG, _IMGT)),
    #         sep='\t', header=False, index=True)
    #
    #     # HLA_DICTIONARY_SNPS_map.to_csv(
    #     #     os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
    #     #     sep='\t', header=False, index=False)


    # if MAKING_SNPS_MAPFILE:
    #
    #     ########## < 4. Making SNPS Map file. > ##########
    #
    #     print(std_MAIN_PROCESS_NAME + "[4] Making SNPS Map file.")
    #
    #     HLA_DICTIONARY_SNPS_map = pd.concat([MakeMap(_inputfile=_OUTPUT + "_{0}.SNPS.forMAP.txt".format(HLA_names[i]),
    #                                               _TYPE="SNPS", _HLA=HLA_names[i], _return_as_dataframe=True)
    #                                       for i in range(0, len(HLA_names))], axis=0)
    #
    #     HLA_DICTIONARY_SNPS_map.to_csv(
    #         os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}.map'.format(_HG, _IMGT)),
    #         sep='\t', header=False, index=False)


    if CLEAN_UP:

        ########## < 5. Removing unnecessary files. > ##########

        print(std_MAIN_PROCESS_NAME + "[5] Removing unnecessary files")

    print("\n")
    print(std_MAIN_PROCESS_NAME + "Making Dictionary is done.\n")

    return 0






#################### < Core Functions > ####################

def MakeMap(_hla, _type, _df_forMAP):

    if not (_type == "AA" or _type == "SNPS"):
        return -1

    if _type == "AA" and (_df_forMAP.shape[1] != 3):
        return -1

    if _type == "SNPS" and (_df_forMAP.shape[1] != 6):
        return -1

    """
    (1) Chr
    (2) Label
    (3) GD
    (4) Genomic Positions
    """

    _df_forMAP = _df_forMAP.astype(str)

    sr_Chr = pd.Series(["6" for i in range(0, _df_forMAP.shape[0])])
    l_Label = []
    sr_GD = pd.Series(["0" for i in range(0, _df_forMAP.shape[0])])
    sr_GenPos = _df_forMAP.iloc[:, 1]

    # Processing `l_Label`.

    p = re.compile('-?\d+x-?\d+')

    for i in range(0, _df_forMAP.shape[0]):

        # t_rel_pos = str(_df_forMAP.iat[i, 0])
        # t_gen_pos = str(_df_forMAP.iat[i, 1])
        # t_type = str(_df_forMAP.iat[i, 2])

        t_rel_pos = _df_forMAP.iat[i, 0]
        t_gen_pos = _df_forMAP.iat[i, 1]
        t_type = _df_forMAP.iat[i, 2]

        main_label = '_'.join([("INDEL" if bool(p.match(t_rel_pos)) else "AA" if _type == "AA" else "SNPS"),
                               _hla, t_rel_pos, t_gen_pos, t_type])

        if _type == "SNPS":

            additional_label = ["AA"]

            t_AA_rel_pos = _df_forMAP.iat[i, 3]
            t_AA_gen_pos = _df_forMAP.iat[i, 4]
            #             t_AA_type = _df_forMAP.iat[i, 5]

            if t_AA_rel_pos != "nan" and t_AA_rel_pos != "NaN":
                # additional_label.append(str(int(t_AA_rel_pos)))
                additional_label.append(t_AA_rel_pos)

            if t_AA_gen_pos != "nan" and t_AA_gen_pos != "NaN":
                # additional_label.append(str(t_AA_gen_pos))
                additional_label.append(t_AA_gen_pos)

            #             if not pd.isna(t_AA_type):
            #                 additional_label.append(str(t_AA_type))

            if len(additional_label) > 1:
                additional_label = '_'.join(additional_label)
                main_label = '_'.join([main_label, additional_label])

        l_Label.append(main_label)

    df_MAP = pd.concat([sr_Chr, pd.Series(l_Label), sr_GD, sr_GenPos], axis=1)

    # Checking as a file.
    df_MAP.to_csv("Test_MAKEDICTIONARY_Prototype.forMAP.{0}.{1}.txt".format(_hla, _type), sep='\t', header=False, index=False)

    return df_MAP




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
    parser.add_argument("-hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)
    parser.add_argument("-o", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                        required=True)
    parser.add_argument("--type", "-t", help="\nSNPS or AA or Both\n\n", choices=["AA", "SNPS", "BOTH"], metavar="Type",
                        default="BOTH")


    ##### < for Test > #####

    # args = parser.parse_args(["-hg", "18", "-o", "MAKEDICTIONARY_v2/makedictionary", "-imgt", "370"])
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
    MakeDictionary(_HG=args.hg, _OUTPUT=args.o, _IMGT=args.imgt, _TYPE=args.type)
