# -*- coding: utf-8 -*-
import os, sys, re
from os.path import join
import numpy as np
import pandas as pd

from src.util import printDict
import IMGT2Seq.src.processModules as modules

from Bio import pairwise2
from Bio.Seq import Seq

def processNuc(_nuc, _hla, _ref_allele:str, _p_data="IMGT2Seq/data"):

    """
    - Nuc file is just used for generating AA map.
    - Only ref seq is needed.

    (2022.12.11.)
    planned to jointly use the nuc's ref seq. in generating map file, but I postponed it again.
    Using nuc's ref seq to generate map file with AA and SNPS might cause more trouble.

    For now, likewise previous version, the relative position start information will be acquired. (ex. HLA-A: -24)

    - chopAlignmentFile
    -

    For each chunk
    - get AA rel_pos start
    - Only ref seq
        - as codons (Don't join)



    """

    ### Load gen raw seqs
    dict_chunks_gen = modules.chopAlignmentFile(_nuc)
    # printDict(dict_chunks_gen, 2)

    l_cDNA, l_AA_codon, l_dict_codons = func1(dict_chunks_gen)
    print("l_cDNA:\n{}".format(l_cDNA))
    print("l_AA_codon:\n{}".format(l_AA_codon))
    for i in range(len(l_dict_codons)):
        print()
        printDict(l_dict_codons[i], 10)

    print("\n\n< Join betw. codon list. (list.extend) >")
    dict_codons = join_betw_codonList(l_dict_codons)
    printDict(dict_codons, 10)

    print("\n\n< Substitute codon bases. >")
    dict_codons = substituteBase(dict_codons) # Base substituted
    printDict(dict_codons, 10)

    print("\n\n< Get ref codons. >")
    _ref_allele, _ref_codons = getRefCodons(_hla, _ref_allele, dict_codons)
    _ref_seq = ''.join(_ref_codons)
    print("ref allele: {}".format(_ref_allele))
    print("ref codons: {}".format(_ref_codons))
    print("ref seq: {}".format(_ref_seq))


    ## Manual test
    # l_vbar = [item for item in _ref_codons if '|' in item]
    # print(l_vbar)
    # print(len(l_vbar))

    print("\n\n< Trim ref codons for AA (Map0 for AA).>")
    df_Map0_AA = makeMap0_AA(_ref_seq, int(l_AA_codon[0]), _hla, _p_data)

    print("\n\n< Trim ref codons for SNPS (Map0 for SNPS).>")
    df_Map0_SNPS = makeMap0_SNPS(_ref_seq, _hla, _p_data)

    # Check as file.
    df_Map0_AA.to_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_AA.nuc.txt", sep='\t', header=True, index=False)
    df_Map0_SNPS.to_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_SNPS.nuc.txt", sep='\t', header=True, index=False)




    # temp_chunk = dict_chunks_gen[1] # 1st chunk
    #
    # # line1 = temp_chunk[0].split()
    # # print(line1)
    # line2 = temp_chunk[1].split()
    # # print(line2)
    #
    # AA_rel_pos_start = line2[2] # ex. '-32'

    return 0



def func1(_dict_chunks):

    l_cDNA = []
    l_AA_codon = []
    l_dict_codons = []

    for i in range(1, len(_dict_chunks) - 1): # Useless 1st and last chunks.

        temp_chunk = _dict_chunks[i]

        temp_cDNA = temp_chunk[0].split()[1]
        temp_AA_codon = temp_chunk[1].split()[2]
        temp_dict_codons = {}

        for line in temp_chunk[3:]:

            t_line = line.split()

            temp_dict_codons[t_line[0]] = t_line[1:]

        l_cDNA.append(temp_cDNA)
        l_AA_codon.append(temp_AA_codon)
        l_dict_codons.append(temp_dict_codons)


    return l_cDNA, l_AA_codon, l_dict_codons


# def ligateRightEnd_nuc(_dict_): return 0

def join_betw_codonList(_l_dict_codons):

    dict_RETURN = {allele: None for allele in _l_dict_codons[0].keys()}

    for allele in _l_dict_codons[0].keys():

        l_temp = []

        for dict_codons in _l_dict_codons:
            if allele in dict_codons:
                l_temp.extend(dict_codons[allele])

        # print(l_temp)

        dict_RETURN[allele] = l_temp

    return dict_RETURN


def substituteBase(_dict_codons):

    def SubstitutionRule(_base_virtual, _base_ref):
        if _base_virtual == '-':
            return _base_ref
        elif _base_virtual == '*':
            return 'x'
        elif _base_ref == '|':  # Only when 'gen' and 'nuc'.
            return _base_ref
        else:
            return _base_virtual


    dict_RETURN = {}

    iter_dict_codons = iter(_dict_codons.items())

    ref_allele, ref_l_codons = next(iter_dict_codons)


    ## 1st seq. (ref seq)
    dict_RETURN[ref_allele] = ref_l_codons # as it is

    ## 2nd seq. ~
    for allele, l_codons in iter_dict_codons:
        # l_codon_substituted = [SubstitutionRule(base_v, base_ref) for base_v, base_ref in zip(l_codons, ref_l_codons)]

        l_temp = [] # substituted `l_codon`

        for str_virtual, str_ref in zip(l_codons, ref_l_codons):

            str_temp = ''.join([SubstitutionRule(base_v, base_ref) for base_v, base_ref in zip(str_virtual, str_ref)])
            # ex) `str_virtual` = '---' / 'str_ref' = 'ATG'

            l_temp.append(str_temp)

        dict_RETURN[allele] = l_temp

    return dict_RETURN


def getRefCodons(_hla, _ref_allele, _dict_codons:dict):

    _ref_allele = "{}*{}".format(_hla, _ref_allele)

    if _ref_allele not in _dict_codons:
        _ref_allele = list(_dict_codons.keys())[0]

    return _ref_allele, _dict_codons[_ref_allele]


def makeMap0_AA(_ref_seq, _AA_rel_pos_start:int, _hla, _p_data):

    """
    A module equivalent to 'processIndel' in 'processModules.py'.

    (Note)
    - `_ref_seq` contains '|'.
    """

    # l_ref_seq_split = _ref_seq.split('|')
    # print("ref seq split: ", l_ref_seq_split)

    ## INS removal
    p_INS = re.compile(r'\.+')
    _ref_seq_InsRemoved = p_INS.sub('', _ref_seq) # used for `l_annots`.

    ## '|' removal
    l_ref_seq_vbar_split = _ref_seq_InsRemoved.split('|') # for `_ref_codon_AA`.
    _ref_seq_InsRemoved_2 = ''.join(l_ref_seq_vbar_split)
    print("ref seq split (INS removed): ", l_ref_seq_vbar_split)
    print("ref seq joined (INS removed): ", _ref_seq_InsRemoved_2)
    print('Length: ', len(_ref_seq_InsRemoved_2))


    ## ref_codons & annots_AA & AA_rel_pos
    _ref_codons_AA = [_ref_seq_InsRemoved_2[i:i+3] for i in range(0, len(_ref_seq_InsRemoved_2), 3)]
    print("_ref_codons_AA", _ref_codons_AA)
    print(len(_ref_codons_AA))

    l_annots = [(base, "exon{}".format(i+1)) for i in range(len(l_ref_seq_vbar_split)) for base in l_ref_seq_vbar_split[i]]
    l_annots = [l_annots[i][1] for i in range(0, len(l_annots), 3)]
    print("l_annots: ", l_annots)
    print(len(l_annots))

    AA_rel_pos = list(range(_AA_rel_pos_start, 0)) + list(range(1, len(_ref_codons_AA) +1 - abs(_AA_rel_pos_start) )) # The last +1 is for "0" between -1 and 1.
    print("AA_rel_pos", AA_rel_pos)
    print(len(AA_rel_pos))

    ## AA residue
    df_AA_codon_table = pd.read_csv(join(_p_data, "AA_codon_table.txt"), sep='\s+', header=0, dtype=str)
    dict_AA_codon_table = {codon: residue for codon, residue in zip(df_AA_codon_table.iloc[:, 0], df_AA_codon_table.iloc[:, 3])}

    l_residue = [dict_AA_codon_table[codon] if codon in dict_AA_codon_table else codon for codon in _ref_codons_AA]

    ### as DataFrame
    df_Map0_AA = pd.DataFrame.from_dict({
        'AA_rel_pos_Nuc': AA_rel_pos,
        'Codon': _ref_codons_AA,
        'residue': l_residue,
        'Type': l_annots
    })

    print("df_Map0_AA:\n{}\n".format(df_Map0_AA))

    return df_Map0_AA


def makeMap0_SNPS(_ref_seq, _hla, _p_data):
    """

    (Note)
    `_ref_seq` contains '|'.
    """

    # check
    # for idx, base in enumerate(_ref_seq):
    #     print("{}: {}".format(idx, base))
    #     if idx > 200: break

    ## INS removal
    p_INS = re.compile(r'([ACGT])\.+([ACGT])') # Only the dots enclosed with bases. (<=> |...| : it can be saved. (e.g. DQB1's 5th exon)
    """
    Removing only dots enclosed with bases can save dots enclosed with '|', which can be not a INS.
    Practically, This exception is for just HLA-DQB1. (its 5th exon)
    """

    _ref_seq_INSremoved = p_INS.sub(r'\1\2', _ref_seq)
    # _ref_seq_INSremoved = p_INS.sub(r"\1Z\2", _ref_seq)
    print(_ref_seq_INSremoved)

    ## '|' revmoal
    _ref_seq_vbar_split = _ref_seq_INSremoved.split('|') # for `l_annots`
    _ref_seq_INSremoved_2 = ''.join(_ref_seq_vbar_split) # for 'l_SNPS_rel_pos` and `l_codons` (*** No INS and '|')
    print(_ref_seq_INSremoved_2)

    ## length check of `_ref_seq_INSremoved_2`(No INS and '|').
    if len(_ref_seq_INSremoved_2) % 3 != 0:
        print("Nuc reference seq for HLA-{} is not multiple of 3, which can cause frameshift.".format(_hla))

    ## Annots, Codons, SNPS_rel_pos
    l_annots = ['exon{}'.format(i+1) for i in range(len(_ref_seq_vbar_split)) for _ in range(len(_ref_seq_vbar_split[i]))]
    print("l_annots: ", l_annots)
    print(len(l_annots))

    # l_SNPS_rel_pos = list(range(1, len(_ref_seq_INSremoved_2)+1))
    l_SNPS_rel_pos = []
    idx_SNPS_rel_pos = 1
    for base in _ref_seq_INSremoved_2:
        if base == '.':
            l_SNPS_rel_pos.append(np.nan)
        else:
            l_SNPS_rel_pos.append(idx_SNPS_rel_pos)
            idx_SNPS_rel_pos += 1
    sr_SNPS_rel_pos = pd.Series(l_SNPS_rel_pos, name='SNPS_rel_pos_Nuc').interpolate()
    print("l_SNPS_rel_pos", l_SNPS_rel_pos)
    print(len(l_SNPS_rel_pos))


    l_codons = [_ref_seq_INSremoved_2[i:i+3] for i in range(0, len(_ref_seq_INSremoved_2), 3)]
    l_codons_3times = [codon for codon in l_codons for _ in range(3)]
    print("l_codons: ", l_codons)
    print(len(l_codons))
    print("l_codons_3times: ", l_codons_3times)
    print(len(l_codons_3times))

    ## AA residue
    df_AA_codon_table = pd.read_csv(join(_p_data, "AA_codon_table.txt"), sep='\s+', header=0, dtype=str)
    dict_AA_codon_table = {codon: residue for codon, residue in zip(df_AA_codon_table.iloc[:, 0], df_AA_codon_table.iloc[:, 3])}
    # printDict(dict_AA_codon_table, 10)

    l_residue_3times = [dict_AA_codon_table[codon] if codon in dict_AA_codon_table else codon for codon in l_codons_3times]

    ### as DataFrame

    df_Map0_SNPS = pd.DataFrame.from_dict({
        "SNPS_rel_pos_Nuc": sr_SNPS_rel_pos,
        "base": list(_ref_seq_INSremoved_2),
        "Codon": l_codons_3times,
        "residue": l_residue_3times,
        "Type": l_annots
    })

    print("df_Map0_SNPS:\n{}\n".format(df_Map0_SNPS))

    return df_Map0_SNPS


def mergeMap0s(_df_Map0_AA, _df_Map0_SNPS):

    print(_df_Map0_AA)
    print(_df_Map0_SNPS)

    ## df_Map0_AA_3times
    df_Map0_AA_3times = pd.DataFrame([
        [rel_pos, codon, residue, type] \

        for rel_pos, codon, residue, type in zip(
            _df_Map0_AA.iloc[:, 0],
            _df_Map0_AA.iloc[:, 1],
            _df_Map0_AA.iloc[:, 2],
            _df_Map0_AA.iloc[:, 3]
        ) \

        for _ in range(3) # repeat 3 times
    ], columns=_df_Map0_AA.columns)
    print("df_Map0_AA_3times:\n{}\n".format(df_Map0_AA_3times))

    df_merge0 = _df_Map0_SNPS.merge(df_Map0_AA_3times, on=['Codon', 'residue'])
    print("df_merge0:\n{}\n".format(df_merge0))

    df_merge0.to_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_AA+SNPS.nuc.txt", sep='\t', header=True, index=True)

    return 0

def BioSeq_test():

    df_Map0_AA = pd.read_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_AA.txt", sep='\s+', header=0, dtype=str)
    print(df_Map0_AA)
    df_Map0_AA_nuc = pd.read_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_AA.nuc.txt", sep='\s+', header=0, dtype=str)
    print(df_Map0_AA_nuc)

    seq_AA = Seq(''.join(df_Map0_AA.iloc[:, 1].tolist()))
    print(seq_AA)
    seq_AA_nuc = Seq(''.join(df_Map0_AA_nuc.iloc[:, 2].tolist()))
    print(seq_AA_nuc)

    alignments = pairwise2.align.globalxx(seq_AA, seq_AA_nuc)
    for match in alignments:
        print(match)
        print(pairwise2.format_alignment(*match))

    return 0

def BioSeq_test_SNPS():

    df_Map0_SNPS = pd.read_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_SNPS.txt", sep='\s+', header=0, dtype=str)
    df_Map0_SNPS = df_Map0_SNPS[df_Map0_SNPS['Type'].str.match('exon')]
    print(df_Map0_SNPS)

    df_Map0_SNPS_nuc = pd.read_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_SNPS.nuc.txt", sep='\s+', header=0, dtype=str)
    print(df_Map0_SNPS_nuc)

    seq_SNPS = Seq(''.join(df_Map0_SNPS.iloc[:, 3].tolist()))
    print(seq_SNPS)
    seq_SNPS_nuc = Seq(''.join(df_Map0_SNPS_nuc.iloc[:, 1].tolist()))
    print(seq_SNPS_nuc)

    alignments = pairwise2.align.globalmx(seq_SNPS, seq_SNPS_nuc, 100, 10)
    for match in alignments:
        print(match)
        print(pairwise2.format_alignment(*match))

    return 0

if __name__ == '__main__':

    # nuc = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3500/alignments/DQB1_nuc.txt"
    nuc = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3500/alignments/DQB1_nuc.txt"

    # processNuc(nuc, "DQB1", "05:01:01:01", _p_data="/home/wansonchoi/sf_VirtualBox_Share/HATK/IMGT2Seq/data")

    ### mergeMap0s
    # df_merged = mergeMap0s(
    #     pd.read_csv('/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_AA.nuc.txt', sep='\s+', header=0, dtype=str),
    #     pd.read_csv('/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_SNPS.nuc.txt', sep='\s+', header=0, dtype=str)
    # )


    ###
    # BioSeq_test()
    BioSeq_test_SNPS()
    pass