# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join
import pandas as pd
import multiprocessing as mp

from IMGT2Seq.src.NfieldDictionary import NfieldDictionary
from IMGT2Seq.src.GenerateHAT import GenerateHAT
from IMGT2Seq.src.ProcessIMGT import ProcessIMGT
from IMGT2Seq.src.AvailableHLAs import getTargetProtFiles
from src.util import Exists
from src.HATK_Error import HATK_InputPreparation_Error

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


def IMGT2Seq(_imgt, _hg, _out_dir, _imgt_dir, _no_Indel=False, _MultiP=False, _f_save_intermediates=False,
             _no_prime=True, _p_data='./data', __Nfield_OUTPUT_FORMAT=4,
             _HLA_target=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

    ### Main Variables ###
    dict_prot_files = getTargetProtFiles(_HLA_target, _imgt_dir, "prot")
    dict_nuc_files = getTargetProtFiles(_HLA_target, _imgt_dir, "nuc")
    dict_gen_files = getTargetProtFiles(_HLA_target, _imgt_dir, "gen")

    d_MapTables = {hla: None for hla in _HLA_target}

    p_allelelist_2009 = join(_imgt_dir, "Nomenclature_2009.txt")
    p_allelelist = join(_imgt_dir, "Allelelist.txt")
    p_wmda_Ggroup = join(_imgt_dir, "wmda/hla_nom_g.txt")
    p_wmda_Pgroup = join(_imgt_dir, "wmda/hla_nom_p.txt")

    _out_prefix_AA = join(_out_dir, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}'.format(_hg, _imgt))
    _out_prefix_SNPS = join(_out_dir, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}'.format(_hg, _imgt))
    _out_prefix_HAT = join(_out_dir, 'HLA_ALLELE_TABLE')


    ### Main Actions ###
    if len(_HLA_target) == 0:
        raise HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "No target HLA to generate a sequence dictionary.\n"
            "Please check the '--HLA' argument again."
        )

    checkRequiredFiles(_imgt_dir, dict_prot_files, dict_nuc_files, dict_gen_files, p_allelelist_2009,
                       p_allelelist, p_wmda_Ggroup, p_wmda_Pgroup) # Dependency check


    _1_MAKING_DICTIONARY = 1
    _2_GENERATE_HAT = 1
    CLEAN_UP = 0


    if _1_MAKING_DICTIONARY:

        ########## < 1. Making HLA Dictionary. > ##########

        # print(std_MAIN_PROCESS_NAME + "[1] Making HLA dictionary file.")
        dict_temp = {hla: None for hla in _HLA_target}

        if _MultiP > 1: # Parallel

            print(std_MAIN_PROCESS_NAME + "Multiprocessing.")

            pool = mp.Pool(processes=_MultiP)

            dict_Pool = {_HLA_target[i]: pool.apply_async(ProcessIMGT, (_out_dir, _HLA_target[i], _hg, _imgt,
                                                                        dict_nuc_files[_HLA_target[i]], dict_gen_files[_HLA_target[i]], dict_prot_files[_HLA_target[i]],
                                                                        _p_data, _no_Indel, _f_save_intermediates))
                         for i in range(0, len(_HLA_target))}

            pool.close()
            pool.join()


            for hla in _HLA_target:
                try:
                    dict_temp[hla] = dict_Pool[hla].get()
                except TypeError:
                    print(std_WARNING_MAIN_PROCESS_NAME +
                          "Generating sequence information dictionary for HLA-{} failed.".format(hla))
                    dict_temp[hla] = None
                else:
                    pass

        else: # Sequential

            for hla in _HLA_target:
                dict_temp[hla] = \
                    ProcessIMGT(_out_dir, hla, _hg, _imgt, dict_nuc_files[hla], dict_gen_files[hla],
                                dict_prot_files[hla], _p_data, _no_Indel=_no_Indel,
                                _save_intermediates=_f_save_intermediates)
        # done.


        l_df_Seqs_AA = []
        l_df_forMAP_AA = []

        l_df_Seqs_SNPS = []
        l_df_forMAP_SNPS = []

        for i in range(0, len(_HLA_target)):

            if not dict_temp[_HLA_target[i]]:
                print(std_MAIN_PROCESS_NAME +
                      "Generating sequence information dictionary for HLA-{} failed.".format(_HLA_target[i]))
                continue


            t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE = dict_temp[_HLA_target[i]]


            # Dictionary AA.
            l_df_Seqs_AA.append(t_df_Seqs_AA)
            l_df_forMAP_AA.append(MakeMap(_HLA_target[i], "AA", t_df_forMAP_AA))

            # Dictionary SNPS.
            l_df_Seqs_SNPS.append(t_df_Seqs_SNPS)
            l_df_forMAP_SNPS.append(MakeMap(_HLA_target[i], "SNPS", t_df_forMAP_SNPS))

            # raw MapTable
            d_MapTables[_HLA_target[i]] = t_MAPTABLE


        # Amino acid sequence dictionaries
        HLA_DICTIONARY_AA = pd.concat(l_df_Seqs_AA, axis=0)
        HLA_DICTIONARY_AA_map = pd.concat(l_df_forMAP_AA, axis=0)

        # DNA sequence dictionaries
        HLA_DICTIONARY_SNPS = pd.concat(l_df_Seqs_SNPS, axis=0)
        HLA_DICTIONARY_SNPS_map = pd.concat(l_df_forMAP_SNPS, axis=0)


        ### Finalizing all output.

        # Exporting AA dictionary.
        HLA_DICTIONARY_AA.to_csv(_out_prefix_AA + ".txt", sep='\t', header=False, index=True)
        HLA_DICTIONARY_AA_map.to_csv(_out_prefix_AA + ".map", sep='\t', header=False, index=False)

        # Exporting SNPS dictionary.
        HLA_DICTIONARY_SNPS.to_csv(_out_prefix_SNPS+".txt", sep='\t', header=False, index=True)
        HLA_DICTIONARY_SNPS_map.to_csv(_out_prefix_SNPS+".map", sep='\t', header=False, index=False)



        if 0 < __Nfield_OUTPUT_FORMAT < 4:

            NfieldDictionary(_out_prefix_AA + ".txt", _out_prefix_SNPS+".txt", d_MapTables, __Nfield_OUTPUT_FORMAT,
                             _HLA_target=_HLA_target)



    if _2_GENERATE_HAT:

        ########## < 2. Making *.hat file. > ##########

        _out_prefix_HAT = \
            GenerateHAT(p_allelelist_2009, p_allelelist, p_wmda_Ggroup, p_wmda_Pgroup, _imgt, _out_prefix_HAT)



    if CLEAN_UP:

        ########## < 5. Removing unnecessary files. > ##########
        pass


    return _out_prefix_AA, _out_prefix_SNPS, _out_prefix_HAT, d_MapTables




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


        p = re.compile('-?\d+x-?\d+') # The pattern for "INDEL"

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


def checkRequiredFiles(_imgt_dir, _dict_prot, _dict_nuc, _dict_gen,
                       _p_allelelist_2009, _p_allelelist, _p_wmda_Ggroup, _p_wmda_Pgroup):

    # # prot
    # for hla, file in _dict_prot.items():
    #     if not Exists(file):
    #         raise HATK_InputPreparation_Error(
    #             std_ERROR_MAIN_PROCESS_NAME +
    #             "The prot file for HLA-{} gene can't be found in '{}'.".format(hla, _imgt_dir+'/alignments')
    #         )
    # # nuc
    # for hla, file in _dict_nuc.items():
    #     if not Exists(file):
    #         raise HATK_InputPreparation_Error(
    #             std_ERROR_MAIN_PROCESS_NAME +
    #             "The nuc file for HLA-{} gene can't be found in '{}'.".format(hla, _imgt_dir+'/alignments')
    #         )
    # # gen
    # for hla, file in _dict_gen.items():
    #     if not Exists(file):
    #         raise HATK_InputPreparation_Error(
    #             std_ERROR_MAIN_PROCESS_NAME +
    #             "The gen file for HLA-{} gene can't be found in '{}'.".format(hla, _imgt_dir+'/alignments')
    #         )

    # 'Nomenclature_2009.txt'
    if not Exists(_p_allelelist_2009):
        raise HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "The 'Nomenclature_2009.txt' file can't be found in '{}'.".format(_imgt_dir)
        )

    # 'Allelelist.txt'
    if not Exists(_p_allelelist):
        raise HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "The 'Allelelist.txt' file can't be found in '{}'.".format(_imgt_dir)
        )

    # 'wmda/hla_nom_g.txt'
    if not Exists(_p_wmda_Ggroup):
        raise HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "The 'wmda/hla_nom_g.txt' file can't be found in '{}'.".format(_imgt_dir)
        )

    # 'wmda/hla_nom_p.txt'
    if not Exists(_p_wmda_Pgroup):
        raise HATK_InputPreparation_Error(
            std_ERROR_MAIN_PROCESS_NAME +
            "The 'wmda/hla_nom_p.txt' file can't be found in '{}'.".format(_imgt_dir)
        )