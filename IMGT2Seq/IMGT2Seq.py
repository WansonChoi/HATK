# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join, exists
import pandas as pd
from glob import glob
import multiprocessing as mp

from IMGT2Seq.src.NfieldDictionary import NfieldDictionary
from IMGT2Seq.src.GenerateHAT import GenerateHAT
# from IMGT2Seq.src.ProcessIMGT import ProcessIMGT
from IMGT2Seq.src.ProcessIMGTv2 import ProcessIMGTv2
from IMGT2Seq.src.XgroupDictionary import XgroupDictionary
from IMGT2Seq.src.IMGT2SeqError import IMGT2SeqError


########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (basename(__file__))

# HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
# raw_HLA_names = ["A", "B", "C", "DPA", "DPB", "DQA", "DQB", "DRB"]



def IMGT2Seq(_out, _imgt_dir, _imgt, _hg, _HLA, _Nfield_OUTPUT_FORMAT=4,
             _p_data="IMGT2Seq/data", _multiprocess=1,
             _no_Ins=False, _include_UTR=False, _save_intermediates=False):

    """

    (1) Generate the HAT (HLA Allele Table)
    (2) check available HLA genes
    (3) Generate dictionary and maptable files restricted to the available HLA genes.


    """

    ### raw MapTables
    d_MapTables = {}

    ### OUTPUT prefix
    _out_dir = dirname(_out)
    os.makedirs(_out_dir, exist_ok=True)

    _OUTPUT_AA_RETURN = join(_out_dir, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}'.format(_hg, _imgt))
    _OUTPUT_SNPS_RETURN = join(_out_dir, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}'.format(_hg, _imgt))
    _OUTPUT_HAT = join(_out_dir, 'HLA_ALLELE_TABLE')





    ########## < 1. Making *.hat file. > ##########

    # Input file check. ("Nomenclature_2009.txt", "Allelelist.txt", "hla_nom_g.txt", "hla_nom_p.txt")
    t_allelelist_2009 = join(_imgt_dir, "Nomenclature_2009.txt")
    t_allelelist = join(_imgt_dir, "Allelelist.txt")
    t_p_Group = join(_imgt_dir, "wmda/hla_nom_g.txt")
    t_p_Proup = join(_imgt_dir, "wmda/hla_nom_p.txt")

    if not exists(t_allelelist_2009):
        raise IMGT2SeqError(
            std_ERROR_MAIN_PROCESS_NAME +
            "'Nomenclature_2009.txt' file can't be found.('{}')\n".format(t_allelelist_2009)
        )

    if not exists(t_allelelist):
        raise IMGT2SeqError(
            std_ERROR_MAIN_PROCESS_NAME +
            "'Allelelist.txt' file can't be found.('{}')\n".format(t_allelelist)
        )

    if not exists(t_p_Group):
        raise IMGT2SeqError(
            std_ERROR_MAIN_PROCESS_NAME +
            "'wmda/hla_nom_g.txt' file can't be found.('{}')\n".format(t_p_Group)
        )

    if not exists(t_p_Proup):
        raise IMGT2SeqError(
            std_ERROR_MAIN_PROCESS_NAME +
            "'wmda/hla_nom_p.txt' file can't be found.('{}')\n".format(t_p_Proup)
        )


    # HAT file generation.
    __HAT__, l_HLA_IMGT_available = GenerateHAT(t_allelelist_2009, t_allelelist, t_p_Group, t_p_Proup, _imgt, _OUTPUT_HAT)
    # print(l_HLA_IMGT_available)





    ########## < 2. Available HLA genes. > ##########

    # (1) get 'IMGT_avail' and 'BP_avail' HLA genes.
    df_BP_start_codon = pd.read_csv(join(_p_data, 'HLA_EXON1_START_CODON_POSITIONS_hg{}.txt'.format(_hg)),
                                    sep='\s+', header=None, dtype=str, names=['HLA_BP_available', 'strand', 'BP_start_codon'])
    # print("df_BP_start_codon:\n{}\n".format(df_BP_start_codon))


    # Union of 'IMGT_avail', 'BP_avail' and 'HLA requested'
    t_set_HLA_IMGT_available = set(l_HLA_IMGT_available)
    t_set_HLA_BP_available = set(df_BP_start_codon['HLA_BP_available'])
    t_set_requested = set(_HLA)

    l_union = list(
        t_set_HLA_IMGT_available | t_set_HLA_BP_available | t_set_requested
    )
    l_union.sort()

    sr_union = pd.Series(l_union, name='HLA')
    # print("sr_union:\n{}\n".format(sr_union))

    # as a Table
    df_HLA_avail_table = pd.DataFrame(
        {'HLA': sr_union,
         'HLA_IMGT_avail' : sr_union.isin(l_HLA_IMGT_available),
         'HLA_BP_avail' : sr_union.isin(df_BP_start_codon['HLA_BP_available'])
         }
    ).set_index('HLA')
    # print("df_HLA_avail_table:\n{}\n".format(df_HLA_avail_table))

    df_HLA_avail_table['HLA_both_avail'] = df_HLA_avail_table.apply(lambda x : x.all(), axis=1)
    # print("df_HLA_avail_table2:\n{}\n".format(df_HLA_avail_table))


    # (2) All requested HLA genes in the both avail?

    sr_HLA_requested = pd.Series(_HLA, name='HLA_requested')
    l_both_avail = df_HLA_avail_table[df_HLA_avail_table['HLA_both_avail']].index.tolist()
    # print("l_both_avail: ", l_both_avail)

    f_all_is_well = sr_HLA_requested.isin(l_both_avail)

    ## Main notification.
    if not f_all_is_well.all():
        print(std_WARNING_MAIN_PROCESS_NAME +
              "Next requested HLA genes can't be included in the output: {} (requested : {})" \
              .format(sr_HLA_requested[~f_all_is_well].tolist(), _HLA))
        print("\n{}\n".format(df_HLA_avail_table.loc[sr_HLA_requested[~f_all_is_well], :]))


    ## Subsetted Variables for the next step.
    HLA_names = sr_HLA_requested[f_all_is_well].tolist()
    df_BP_start_codon = df_BP_start_codon[df_BP_start_codon['HLA_BP_available'].isin(HLA_names)] \
                            .set_index('HLA_BP_available')
    # print("df_BP_start_codon(Subsetted):\n{}\n".format(df_BP_start_codon))





    ########## < 3. Making HLA Dictionary. > ##########

    ## Required input file check. (*.gen, *.nuc, *.prot files)
    for hla in HLA_names:

        f_gen = join(_imgt_dir, "alignments/{}_gen.txt".format(hla))
        f_nuc = join(_imgt_dir, "alignments/{}_nuc.txt".format(hla))
        f_prot = join(_imgt_dir, "alignments/{}_prot.txt".format(hla))

        if not exists(f_gen):
            raise IMGT2SeqError(
                std_ERROR_MAIN_PROCESS_NAME +
                "'{}' file can't be found in '{}' directory. Please check it again." \
                .format(f_gen, join(_imgt_dir, "alignments"))
            )

        if not exists(f_nuc):
            raise IMGT2SeqError(
                std_ERROR_MAIN_PROCESS_NAME +
                "'{}' file can't be found in '{}' directory. Please check it again." \
                .format(f_nuc, join(_imgt_dir, "alignments"))
            )

        if not exists(f_prot):
            raise IMGT2SeqError(
                std_ERROR_MAIN_PROCESS_NAME +
                "'{}' file can't be found in '{}' directory. Please check it again." \
                .format(f_prot, join(_imgt_dir, "alignments"))
            )





    l_df_Seqs_AA = []
    l_df_MAP_AA = []

    l_df_Seqs_SNPS = []
    l_df_MAP_SNPS = []


    if _multiprocess > 1:
        pass
    else:

        for hla in HLA_names:

            print("Processing HLA-{} (Serial).".format(hla))

            t_gen = join(_imgt_dir, join(_imgt_dir, "alignments/{}_gen.txt".format(hla)))
            t_nuc = join(_imgt_dir, join(_imgt_dir, "alignments/{}_nuc.txt".format(hla)))
            t_prot = join(_imgt_dir, join(_imgt_dir, "alignments/{}_prot.txt".format(hla)))

            t_BP_start_codon = df_BP_start_codon.loc[hla, 'BP_start_codon']
            t_isREVERSE = (df_BP_start_codon.loc[hla, 'strand'] == '-')


            t_output = ProcessIMGTv2(
                _out, hla, _imgt, t_BP_start_codon, t_isREVERSE,
                t_gen, t_nuc, t_prot,
                _no_Ins=_no_Ins, _include_UTR=_include_UTR,
                _save_intermediates=_save_intermediates
            )

            t_df_gen_dictionary, t_df_prot_dictionary, t_df_SNPS_MAP, t_df_AA_MAP, t_df_prot_markers = t_output

            print("t_df_gen_dictionary:\n{}\n".format(t_df_gen_dictionary))
            print("t_df_SNPS_MAP:\n{}\n".format(t_df_SNPS_MAP))


    # Gathering... (Need to refine ProcessIMGTv2 more...)







    sys.exit()
    ########## < 3. Making HLA Dictionary. > ##########

    # print(std_MAIN_PROCESS_NAME + "[1] Making HLA dictionary file.")

    l_df_Seqs_AA = []
    l_df_forMAP_AA = []

    l_df_Seqs_SNPS = []
    l_df_forMAP_SNPS = []

    if _multiprocess > 1:
        print(std_MAIN_PROCESS_NAME + "Multiprocessing.")

        pool = mp.Pool(processes=_multiprocess)

        dict_Pool = {HLA_names[i]: pool.apply_async(ProcessIMGT, (_out, HLA_names[i], _hg, _imgt,
                                                                  TARGET_nuc_files[HLA_names[i]],
                                                                  TARGET_gen_files[HLA_names[i]],
                                                                  TARGET_prot_files[HLA_names[i]],
                                                                  _p_data, _no_Ins, _save_intermediates))
                     for i in range(0, len(HLA_names))}

        pool.close()
        pool.join()

    for i in range(0, len(HLA_names)):

        if _multiprocess > 1:

            t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE = \
                dict_Pool[HLA_names[i]].get()

        else:
            t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE \
                = ProcessIMGT(_out, HLA_names[i], _hg, _imgt, TARGET_nuc_files[HLA_names[i]],
                              TARGET_gen_files[HLA_names[i]], TARGET_prot_files[HLA_names[i]], _p_data,
                              _no_Indel=_no_Ins, _save_intermediates=_save_intermediates)

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


    if 0 < _Nfield_OUTPUT_FORMAT < 4:
        # 1,2,3-field
        NfieldDictionary(_OUTPUT_AA_RETURN + ".txt", _OUTPUT_SNPS_RETURN + ".txt", d_MapTables,
                         _Nfield_OUTPUT_FORMAT)

    elif _Nfield_OUTPUT_FORMAT == 5:
        # Ggroup
        XgroupDictionary(__HAT__, _OUTPUT_AA_RETURN + ".txt", _OUTPUT_SNPS_RETURN + ".txt", d_MapTables,
                         _out=None, _OUTPUT_FORMAT=_Nfield_OUTPUT_FORMAT) # Overide

    elif _Nfield_OUTPUT_FORMAT == 6:
        # Pgroup
        XgroupDictionary(__HAT__, _OUTPUT_AA_RETURN + ".txt", _OUTPUT_SNPS_RETURN + ".txt", d_MapTables,
                         _out=None, _OUTPUT_FORMAT=_Nfield_OUTPUT_FORMAT) # Overide




    return [_OUTPUT_AA_RETURN, _OUTPUT_SNPS_RETURN, __HAT__, d_MapTables]





# #################### < Core Functions > ####################
#
#
#
# def MakeMap(_hla, _type, _df_forMAP, _no_prime=True):
#
#     if not (_type == "AA" or _type == "SNPS"):
#         return -1
#
#     # if _type == "AA" and (_df_forMAP.shape[1] != 3):
#     #     return -1
#     #
#     # if _type == "SNPS" and (_df_forMAP.shape[1] != 6):
#     #     return -1
#
#     """
#     (1) Chr
#     (2) Label
#     (3) GD
#     (4) Genomic Positions
#     """
#
#     if not _no_prime:
#
#         _df_forMAP = _df_forMAP.astype(str)
#
#         sr_Chr = pd.Series(["6" for i in range(0, _df_forMAP.shape[0])])
#         l_Label = []
#         sr_GD = pd.Series(["0" for i in range(0, _df_forMAP.shape[0])])
#         sr_GenPos = _df_forMAP.iloc[:, 1]
#
#         # Processing `l_Label`.
#
#         p = re.compile('-?\d+x-?\d+')
#
#         for i in range(0, _df_forMAP.shape[0]):
#
#             t_rel_pos = _df_forMAP.iat[i, 0]
#             t_gen_pos = _df_forMAP.iat[i, 1]
#             t_type = _df_forMAP.iat[i, 2]
#
#             main_label = '_'.join([("INS" if bool(p.match(t_rel_pos)) else "AA" if _type == "AA" else "SNPS"),
#                                    _hla, t_rel_pos, t_gen_pos, t_type])
#
#             if _type == "SNPS":
#
#                 additional_label = ["AA"]
#
#                 t_AA_rel_pos = _df_forMAP.iat[i, 3]
#                 t_AA_gen_pos = _df_forMAP.iat[i, 4]
#                 #             t_AA_type = _df_forMAP.iat[i, 5]
#
#                 if t_AA_rel_pos != "nan" and t_AA_rel_pos != "NaN":
#                     additional_label.append(t_AA_rel_pos)
#
#                 if t_AA_gen_pos != "nan" and t_AA_gen_pos != "NaN":
#                     additional_label.append(t_AA_gen_pos)
#
#                 if len(additional_label) > 1:
#                     additional_label = '_'.join(additional_label)
#                     main_label = '_'.join([main_label, additional_label])
#
#             l_Label.append(main_label)
#
#         df_MAP = pd.concat([sr_Chr, pd.Series(l_Label), sr_GD, sr_GenPos], axis=1)
#
#
#
#
#     else:
#
#         sr_Chr = pd.Series(("6" for i in range(0, _df_forMAP.shape[0])))
#         l_Label = []
#         sr_GD = pd.Series(("0" for i in range(0, _df_forMAP.shape[0])))
#         sr_GenPos = _df_forMAP.iloc[:, 1]
#
#         p = re.compile('-?\d+x-?\d+')  # The pattern for "INDEL"
#
#         for row in _df_forMAP.astype(str).itertuples():
#
#             """
#             row[0] := index
#
#             row[1] := rel_pos
#             row[2] := gen_pos
#             row[3] := type
#
#             """
#
#             # t = "INDEL" if bool(p.match(row[1])) else "AA" if _type == "AA" else "SNPS" # (2019. 03. 03.) SNP2HLA README.v2.txt수정하다가, Indel label 때문에 수정하기 전.
#
#             if bool(p.match(row[1])):
#
#                 if _type == "AA":
#                     t = "INS_AA"
#                 elif _type == "SNPS":
#                     t = "INS_SNPS"
#
#                 main_label = '_'.join([t, _hla, *row[1:-1]])
#
#
#             else:
#
#                 if _type == "AA":
#                     t = "AA"
#                 elif _type == "SNPS":
#                     t = "SNPS"
#
#                 main_label = '_'.join([t, _hla, *row[1:]])
#
#             l_Label.append(main_label)
#
#         l_Label = pd.Series(l_Label)
#
#         # Output map file.
#         df_MAP = pd.concat([sr_Chr, l_Label, sr_GD, sr_GenPos], axis=1)
#
#     return df_MAP