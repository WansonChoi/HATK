# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd
from glob import glob
import multiprocessing as mp

from IMGT2Seq.src.NfieldDictionary import NfieldDictionary
from IMGT2Seq.src.GenerateHAT import GenerateHAT
from IMGT2Seq.src.ProcessIMGT import ProcessIMGT
from IMGT2Seq.src.XgroupDictionary import XgroupDictionary


########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
raw_HLA_names = ["A", "B", "C", "DPA", "DPB", "DQA", "DQB", "DRB"]



def IMGT2Seq(_imgt, _hg, _out, _imgt_dir, _no_Indel=False, _MultiP=False, _save_intermediates=False, _no_prime=True,
             _p_data="IMGT2Seq/data", __Nfield_OUTPUT_FORMAT=4):

    ### Dictionaries for Raw files.
    TARGET_prot_files = {}
    TARGET_nuc_files = {}
    TARGET_gen_files = {}

    ### raw MapTables
    d_MapTables = {}

    ### OUTPUT prefix
    # Preparing intermediate paths.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)):
        INTERMEDIATE_PATH = os.path.dirname(_out)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)

    _OUTPUT_AA_RETURN = os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}'.format(_hg, _imgt))
    _OUTPUT_SNPS_RETURN = os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}'.format(_hg, _imgt))
    _OUTPUT_HAT = os.path.join(INTERMEDIATE_PATH, 'HLA_ALLELE_TABLE')



    ########## < Dependency Checking > ##########

    # IMGTHLA directory
    IMGTHLA_directory = _imgt_dir

    if not os.path.isdir(IMGTHLA_directory):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given IMGT-HLA directory \"{0}\" can't be found!".format(IMGTHLA_directory))
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
                print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}_prot.txt\" or \"{0}_nuc.txt\" file can't be found.".format(
                    HLA_names[i]))
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



    ##### < "Allelelist.txt", "hla_nom_g.txt", "hla_nom_p.txt" > #####

    t_allelelist_2009 = os.path.join(IMGTHLA_directory, "Nomenclature_2009.txt")
    t_allelelist = os.path.join(IMGTHLA_directory, "Allelelist.txt")
    t_p_Group = os.path.join(IMGTHLA_directory, "wmda/hla_nom_g.txt")
    t_p_Proup = os.path.join(IMGTHLA_directory, "wmda/hla_nom_p.txt")

    if not os.path.exists(t_allelelist_2009):
        print(std_ERROR_MAIN_PROCESS_NAME + "'Nomenclature_2009.txt' file can't be found.\n")
        sys.exit()

    if not os.path.exists(t_allelelist):
        print(std_ERROR_MAIN_PROCESS_NAME + "'Allelelist.txt' file can't be found.\n")
        sys.exit()

    if not os.path.exists(t_p_Group):
        print(std_ERROR_MAIN_PROCESS_NAME + "'wmda/hla_nom_g.txt' file can't be found.\n")
        sys.exit()

    if not os.path.exists(t_p_Proup):
        print(std_ERROR_MAIN_PROCESS_NAME + "'wmda/hla_nom_p.txt' file can't be found.\n")
        sys.exit()



    _1_GENERATE_HAT = 1
    _2_MAKING_DICTIONARY = 1
    CLEAN_UP = 0



    if _1_GENERATE_HAT:
        ########## < 2. Making *.hat file. > ##########

        _OUTPUT_HAT = GenerateHAT(t_allelelist_2009, t_allelelist, t_p_Group, t_p_Proup, _imgt, _OUTPUT_HAT)



    if _2_MAKING_DICTIONARY:

        ########## < 1. Making HLA Dictionary. > ##########

        # print(std_MAIN_PROCESS_NAME + "[1] Making HLA dictionary file.")

        l_df_Seqs_AA = []
        l_df_forMAP_AA = []

        l_df_Seqs_SNPS = []
        l_df_forMAP_SNPS = []

        if _MultiP > 1:
            print(std_MAIN_PROCESS_NAME + "Multiprocessing.")

            pool = mp.Pool(processes=_MultiP)

            dict_Pool = {HLA_names[i]: pool.apply_async(ProcessIMGT, (_out, HLA_names[i], _hg, _imgt,
                                                                      TARGET_nuc_files[HLA_names[i]],
                                                                      TARGET_gen_files[HLA_names[i]],
                                                                      TARGET_prot_files[HLA_names[i]],
                                                                      _p_data, _no_Indel, _save_intermediates))
                         for i in range(0, len(HLA_names))}

            pool.close()
            pool.join()

        for i in range(0, len(HLA_names)):

            if _MultiP > 1:

                t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE = \
                    dict_Pool[HLA_names[i]].get()

            else:
                t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE \
                    = ProcessIMGT(_out, HLA_names[i], _hg, _imgt,
                                  TARGET_nuc_files[HLA_names[i]], TARGET_gen_files[HLA_names[i]],
                                  TARGET_prot_files[HLA_names[i]],
                                  _p_data, _no_Indel=_no_Indel, _save_intermediates=_save_intermediates)

            # Dictionary AA.
            l_df_Seqs_AA.append(t_df_Seqs_AA)
            l_df_forMAP_AA.append(MakeMap(HLA_names[i], "AA", t_df_forMAP_AA))

            # Dictionary SNPS.
            l_df_Seqs_SNPS.append(t_df_Seqs_SNPS)
            l_df_forMAP_SNPS.append(MakeMap(HLA_names[i], "SNPS", t_df_forMAP_SNPS))

            # raw MapTable
            d_MapTables[HLA_names[i]] = t_MAPTABLE

        # Amino acid sequence dictionaries
        HLA_DICTIONARY_AA = pd.concat(l_df_Seqs_AA, axis=0)
        HLA_DICTIONARY_AA_map = pd.concat(l_df_forMAP_AA, axis=0)

        # DNA sequence dictionaries
        HLA_DICTIONARY_SNPS = pd.concat(l_df_Seqs_SNPS, axis=0)
        HLA_DICTIONARY_SNPS_map = pd.concat(l_df_forMAP_SNPS, axis=0)

        ### Finalizing all output.

        # Exporting AA dictionary.
        HLA_DICTIONARY_AA.to_csv(_OUTPUT_AA_RETURN + ".txt", sep='\t', header=False, index=True)
        HLA_DICTIONARY_AA_map.to_csv(_OUTPUT_AA_RETURN + ".map", sep='\t', header=False, index=False)

        # Exporting SNPS dictionary.
        HLA_DICTIONARY_SNPS.to_csv(_OUTPUT_SNPS_RETURN + ".txt", sep='\t', header=False, index=True)
        HLA_DICTIONARY_SNPS_map.to_csv(_OUTPUT_SNPS_RETURN + ".map", sep='\t', header=False, index=False)


        if 0 < __Nfield_OUTPUT_FORMAT < 4:
            # 1,2,3-field
            NfieldDictionary(_OUTPUT_AA_RETURN + ".txt", _OUTPUT_SNPS_RETURN + ".txt", d_MapTables,
                             __Nfield_OUTPUT_FORMAT)

        elif __Nfield_OUTPUT_FORMAT == 5:
            # Ggroup
            XgroupDictionary(_OUTPUT_HAT, _OUTPUT_AA_RETURN + ".txt", _OUTPUT_SNPS_RETURN + ".txt", d_MapTables,
                             _out=None, _OUTPUT_FORMAT=__Nfield_OUTPUT_FORMAT) # Overide

        elif __Nfield_OUTPUT_FORMAT == 6:
            # Pgroup
            XgroupDictionary(_OUTPUT_HAT, _OUTPUT_AA_RETURN + ".txt", _OUTPUT_SNPS_RETURN + ".txt", d_MapTables,
                             _out=None, _OUTPUT_FORMAT=__Nfield_OUTPUT_FORMAT) # Overide



    if CLEAN_UP:
        ########## < 5. Removing unnecessary files. > ##########
        pass

    return [_OUTPUT_AA_RETURN, _OUTPUT_SNPS_RETURN, _OUTPUT_HAT, d_MapTables]





#################### < Core Functions > ####################



def MakeMap(_hla, _type, _df_forMAP, _no_prime=True):

    if not (_type == "AA" or _type == "SNPS"):
        return -1

    # if _type == "AA" and (_df_forMAP.shape[1] != 3):
    #     return -1
    #
    # if _type == "SNPS" and (_df_forMAP.shape[1] != 6):
    #     return -1

    """
    (1) Chr
    (2) Label
    (3) GD
    (4) Genomic Positions
    """

    if not _no_prime:

        _df_forMAP = _df_forMAP.astype(str)

        sr_Chr = pd.Series(["6" for i in range(0, _df_forMAP.shape[0])])
        l_Label = []
        sr_GD = pd.Series(["0" for i in range(0, _df_forMAP.shape[0])])
        sr_GenPos = _df_forMAP.iloc[:, 1]

        # Processing `l_Label`.

        p = re.compile('-?\d+x-?\d+')

        for i in range(0, _df_forMAP.shape[0]):

            t_rel_pos = _df_forMAP.iat[i, 0]
            t_gen_pos = _df_forMAP.iat[i, 1]
            t_type = _df_forMAP.iat[i, 2]

            main_label = '_'.join([("INS" if bool(p.match(t_rel_pos)) else "AA" if _type == "AA" else "SNPS"),
                                   _hla, t_rel_pos, t_gen_pos, t_type])

            if _type == "SNPS":

                additional_label = ["AA"]

                t_AA_rel_pos = _df_forMAP.iat[i, 3]
                t_AA_gen_pos = _df_forMAP.iat[i, 4]
                #             t_AA_type = _df_forMAP.iat[i, 5]

                if t_AA_rel_pos != "nan" and t_AA_rel_pos != "NaN":
                    additional_label.append(t_AA_rel_pos)

                if t_AA_gen_pos != "nan" and t_AA_gen_pos != "NaN":
                    additional_label.append(t_AA_gen_pos)

                if len(additional_label) > 1:
                    additional_label = '_'.join(additional_label)
                    main_label = '_'.join([main_label, additional_label])

            l_Label.append(main_label)

        df_MAP = pd.concat([sr_Chr, pd.Series(l_Label), sr_GD, sr_GenPos], axis=1)




    else:

        sr_Chr = pd.Series(("6" for i in range(0, _df_forMAP.shape[0])))
        l_Label = []
        sr_GD = pd.Series(("0" for i in range(0, _df_forMAP.shape[0])))
        sr_GenPos = _df_forMAP.iloc[:, 1]

        p = re.compile('-?\d+x-?\d+')  # The pattern for "INDEL"

        for row in _df_forMAP.astype(str).itertuples():

            """
            row[0] := index

            row[1] := rel_pos
            row[2] := gen_pos
            row[3] := type

            """

            # t = "INDEL" if bool(p.match(row[1])) else "AA" if _type == "AA" else "SNPS" # (2019. 03. 03.) SNP2HLA README.v2.txt수정하다가, Indel label 때문에 수정하기 전.

            if bool(p.match(row[1])):

                if _type == "AA":
                    t = "INS_AA"
                elif _type == "SNPS":
                    t = "INS_SNPS"

                main_label = '_'.join([t, _hla, *row[1:-1]])


            else:

                if _type == "AA":
                    t = "AA"
                elif _type == "SNPS":
                    t = "SNPS"

                main_label = '_'.join([t, _hla, *row[1:]])

            l_Label.append(main_label)

        l_Label = pd.Series(l_Label)

        # Output map file.
        df_MAP = pd.concat([sr_Chr, l_Label, sr_GD, sr_GenPos], axis=1)

    return df_MAP