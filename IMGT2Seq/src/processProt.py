# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

from src.util import printDict
import IMGT2Seq.src.processModules as modules

def processProt(_prot, _hla, _rel_pos_start:int, _direction:str, _ref_allele:str):

    print("\n\n< Load raw seq. >")
    dict_chunks = modules.chopAlignmentFile(_prot)

    print("\n\n< Join seqs within chunks >")
    l_joined_within = modules.join_within_chunks_Prot(dict_chunks)
    for i in range(len(l_joined_within)):
        print("{}th chunk:".format(i))
        printDict(l_joined_within[i], 5)

    print("\n\n< Join seqs between chunks >")
    dict_joined_between = modules.join_between_chunks_Prot(l_joined_within)
    printDict(dict_joined_between, 10)
    # print(dict_joined_between['A*03:437Q'])
    # print("{} vs. {}".format(len(dict_joined_between['A*01:01:01:01']), len(dict_joined_between['A*03:437Q'])))

    print("\n\n< Ligate the right end. >")
    dict_ligated = modules.ligateRightEnd(dict_joined_between)
    printDict(dict_ligated, 10)

    print("\n\n< Substitute bases. >")
    dict_substituted = modules.substituteBase(dict_ligated)
    printDict(dict_substituted, 10)

    print("\n\n< Process Indel. >")
    dict_IndelProcessed = modules.processIndel(dict_substituted) # already joined.
    printDict(dict_IndelProcessed, 10)


    _ref_allele = getRefSeq_prot(_hla, _ref_allele, dict_IndelProcessed)
    _ref_seq = dict_IndelProcessed[_ref_allele]

    df_Map0_AA = allocateRelPos_prot(_ref_seq, _rel_pos_start)
    print("df_Map0_AA:\n{}\n".format(df_Map0_AA))

    # export for testing
    df_Map0_AA.to_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_AA.txt", sep='\t', header=True, index=False)

    return 0


def getRefSeq_prot(_hla, _ref_allele, _dict_Indel_Processed):

    _ref_allele = "{}*{}".format(_hla, _ref_allele)

    if _ref_allele in _dict_Indel_Processed.keys():
        return _ref_allele
    else:
        return list(_dict_Indel_Processed.keys())[0]


def allocateRelPos_prot(_ref_seq, _rel_pos_start:int):

    """
    """

    iter_ref_seq = iter(_ref_seq)
    acc_rel_pos = _rel_pos_start

    ## negative part
    l_negative_part = []
    for residue in iter_ref_seq:

        print("{} {}".format(acc_rel_pos, residue))

        if residue == 'z' or residue == 'Z':
            l_negative_part.append('{}x{}'.format(acc_rel_pos-1, acc_rel_pos))
        else:
            l_negative_part.append(acc_rel_pos)
            acc_rel_pos += 1

        if acc_rel_pos == 0: break

    acc_rel_pos += 1

    ## positive part
    l_positive_part = []
    for residue in iter_ref_seq:
        print("{} {}".format(acc_rel_pos, residue))

        if residue == 'z' or residue == 'Z':
            l_positive_part.append('{}x{}'.format(acc_rel_pos-1, acc_rel_pos))
        else:
            l_positive_part.append(acc_rel_pos)
            acc_rel_pos += 1


    l_rel_pos = l_negative_part + l_positive_part

    df_RETURN = pd.concat([
        pd.Series(l_rel_pos, name='AA_rel_pos'),
        pd.Series(list(_ref_seq), name='AA_residue')
    ], axis=1)

    ### check
    # for rel_pos, residue in zip(l_rel_pos, list(_ref_seq)):
    #     print("{} {}".format(rel_pos, residue))

    return df_RETURN


if __name__ == '__main__':

    # A, prot
    prot = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3500/alignments/DQB1_prot.txt"
    # prot = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3500/alignments/A_prot.txt"

    processProt(prot, "DQB1", -32, '-', '05:01:01:01')

    """
    (Prot)
    - Load raw sequences.
    - Join seqs within chunks.
    - Join seqs between chunks.

    - Ligate right end.
    - Substitute '-'(dash), '*', and 'x'.
    - Process indel.

    - Extract the ref. seq for MAP file generation.


    """

    pass