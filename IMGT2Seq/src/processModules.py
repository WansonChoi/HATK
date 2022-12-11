# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

from src.util import printDict


def chopAlignmentFile(_file:str) -> dict:

    """
    Chop given alignment file(ex. *.{prot,gen,nuc}.txt) by rows with '\n'.

    (ex)

    # file: A_prot.txt
    # date: 2022-10-12
    # version: IPD-IMGT/HLA 3.50.0
    # origin: http://hla.alleles.org/wmda/A_prot.txt
    # repository: https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/A_prot.txt
    # author: WHO, Steven G. E. Marsh (steven.marsh@ucl.ac.uk)
                       <= (***) SPLIT HERE.
     Prot              -30                                1
                       |                                  |
     A*01:01:01:01           MAVM APRTLLLLLS GALAL..TQT WAGSHSMRYF FTSVSRPGRG EPRFIAVGYV DDTQFVRFDS DAASQKMEPR APWIEQEGPE YWDQETRNMK
     A*01:01:01:02N          ---- ---------- -----..--- ---------- ---------- ---------- ---------- ---------- ---------- ----------
     A*01:01:01:03           ---- ---------- -----..--- ---------- ---------- ---------- ---------- ---------- ---------- ----------
     A*01:01:01:04           ---- ---------- -----..--- ---------- ---------- ---------- ---------- ---------- ---------- ----------

     (...)

     A*80:06                 ---- P--------- -----..--- ---------- ---------- ---------- --S---Q--- -----R---- -------E-- ---E----V-
     A*80:07                 ---- P--------- -----..--- ---------- ---------- ---------- --S---Q--- -----R---- -------E-- ---E----V-
     A*80:08N                ---- P--------- -----..--- ---------- ---------- ---------- --S---Q--- -----R---- -------E-- ---E----V-
     A*80:09N                **** ********** *****..*** ***------- ---------- ---------- --S---Q--- -----R---- -------E-- ---E----V-
                       <= (***) SPLIT HERE.
     Prot              69
                       |
     A*01:01:01:01     AHSQTDRANL GTLRGYYNQS EDGSHTIQIM YGCDVGPDGR FLRGY.RQDA YDGKDY.IAL NEDLRSWTAA DMAAQITKRK WEAVHAAE.. ..........
     A*01:01:01:02N    ---------- ---------- ---DPGPGRR SRPLIPHGRA RSPTV.SGSE IHPEAA.GLR DPCPGRGPG- FTRFHFQF-P KIPPGWSG.. ..........
     A*01:01:01:03     ---------- ---------- ---------- ---------- -----.---- ------.--- ---------- ---------- --------.. ..........
     A*01:01:01:04     ---------- ---------- ---------- ---------- -----.---- ------.--- ---------- ---------- --------.. ..........

    """

    p_sep = re.compile(r'\s*\n$')

    dict_chunks = {}

    with open(_file, 'r') as f_alignment:

        l_temp_chunk = []

        count_line = 0
        count_chunk = 0

        for line in f_alignment:

            # print(line)
            # print("[{}th]: {}".format(count_line, line))

            if bool(p_sep.match(line)): ## chop

                # print(line)
                # print(l_temp_chunk)

                dict_chunks[count_chunk] = l_temp_chunk.copy()
                l_temp_chunk.clear()
                count_chunk += 1

            else:
                l_temp_chunk.append(line)


            count_line += 1
            # if count_line >= 1000: break


    # print output
    # for idx, chunk in dict_chunks.items():
    #     print("\n{}: ".format(idx))
    #     print(pd.Series(chunk))

    return dict_chunks


"""
Compared to the above chopping job,

trimming should be done differently to prot, gen and nuc.

For example, In prot, there is no '|' between exon and intron, while gen and nuc have.

"""

def join_within_chunks_Prot(_dict_chunks:dict):

    """

    assumes
    (1) the 1st and last chunks are useless.
    (2) In the rest of chunks, the 1st and 2nd lines are useless.
        (e.g.)
    0        Prot              -30                        ...
    1                          |                          ...
    2        A*01:01:01:01           MAVM APRTLLLLLS GALAL...
    3        A*01:01:01:02N          ---- ---------- -----...
    4        A*01:01:01:03           ---- ---------- -----...
                                  ...
    7449     A*80:05                 ---- P--------- -----...
    7450     A*80:06                 ---- P--------- -----...

    => 0 and 1 rows are useless.

    """

    ## To return.
    l_joined_within = [] # a list of dictionaries.

    for i in range(1, len(_dict_chunks) - 1): # No need the 1st and last chunks

        # print("\n{}th chunk.".format(i))

        l_temp = _dict_chunks[i]
        dict_temp = {}


        count = 0

        for j in range(2, len(l_temp)): # No need 1st and 2nd lines.

            line_temp = l_temp[j].split()
            # print(line_temp)

            allele = line_temp[0]
            seq = ''.join(line_temp[1:]) # plain concatenation.

            dict_temp[allele] = seq

            count += 1
            # if count >= 5: break

        # printDict(dict_temp)
        l_joined_within.append(dict_temp.copy())

    """
    (before)
    ['A*01:01:01:01', 'MAVM', 'APRTLLLLLS', 'GALAL..TQT', 'WAGSHSMRYF', 'FTSVSRPGRG', 'EPRFIAVGYV', 'DDTQFVRFDS', 'DAASQKMEPR', 'APWIEQEGPE', 'YWDQETRNMK']
    ['A*01:01:01:02N', '----', '----------', '-----..---', '----------', '----------', '----------', '----------', '----------', '----------', '----------']
    ['A*01:01:01:03', '----', '----------', '-----..---', '----------', '----------', '----------', '----------', '----------', '----------', '----------']
    ['A*01:01:01:04', '----', '----------', '-----..---', '----------', '----------', '----------', '----------', '----------', '----------', '----------']
    ['A*01:01:01:05', '----', '----------', '-----..---', '----------', '----------', '----------', '----------', '----------', '----------', '----------']
    
    (after)
    A*01:01:01:01: MAVMAPRTLLLLLSGALAL..TQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMK
    A*01:01:01:02N: -------------------..-------------------------------------------------------------------------
    A*01:01:01:03: -------------------..-------------------------------------------------------------------------
    A*01:01:01:04: -------------------..-------------------------------------------------------------------------
    A*01:01:01:05: -------------------..-------------------------------------------------------------------------        
    
    """

    return l_joined_within



def join_between_chunks_Prot(_l_chunks:list) -> dict:

    dict_joined_between = {}

    whole_alleles = list(_l_chunks[0].keys())
    # print(whole_alleles)
    # print(len(whole_alleles))
    # print(bool("A*03:437Q" in whole_alleles))

    for allele in whole_alleles:

        str_temp = ''.join([chunk[allele] if allele in chunk else "" for chunk in _l_chunks])

        # print("{}: {}".format(allele, str_temp))
        dict_joined_between[allele] = str_temp


    return dict_joined_between



def ligateRightEnd(_dict_joined_between):

    dict_RETURN = {}

    ### Find the max length.
    allele_max, LEN_max = "", -1
    for allele, seq in _dict_joined_between.items():
        if len(seq) > LEN_max:
            allele_max = allele
            LEN_max = len(seq)

    # print("(Max Seq):\n{}: {}".format(allele_max, LEN_max))



    ### virtual sequences (seqs that were aligned to the ref seq.)
    count = 0
    for allele, seq in _dict_joined_between.items():

        if len(seq) < LEN_max:
            ## short seqs

            tmp_ligated = seq + 'x'*(LEN_max - len(seq))

            # print(
            #     "{}:\n"
            #     "{}\n"
            #     "{}".format(allele, seq, tmp_ligated)
            # )

            dict_RETURN[allele] = tmp_ligated # update with ligated one.

            count += 1
            # if count >= 10: break

        else:
            ## normal ones (as it is.)
            dict_RETURN[allele] = seq


    return dict_RETURN



def substituteBase(_dict_ligated):

    """

    Based on the official definition(https://www.ebi.ac.uk/ipd/imgt/hla/alignment/help/),

    '-': identical to that of the ref seq.
    '*': unknown at any point in the alignment.
    '.': an insertion or deletion
    blank: sequence following the termination codon. (Related to protein seq.)
    'X': the 'Stop' codons.

    In IMGT2Seq,
        - The blank was preprocessed to 'x' in the previous 'ligateRightEnd()' function.
        - '*' will be substituted by 'x'
        - '-' will be substituted by the base of the ref seq.
        - '.' will be processed after this function.


    Sometimes, the 1st seq is not ref seq. (ex. imgt3320, DQB1 prot: DQB1*05:01:01:01 / gen: DQB1*02:01:01)
    But, the 1st seq must be used as (a sort of) ref seq anyway. (i.e. the 1st seq != ref seq)

    In making map file, the exact reference seq will be used.

    """

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

    iter_dict_ligated = iter(_dict_ligated.items())

    ### reference seq. (the 1st seq, practically.)
    ref_allele, ref_seq = next(iter_dict_ligated)
    # print("ref allele: {}".format(ref_allele))
    # print("ref seq: {}".format(ref_seq))

    dict_RETURN[ref_allele] = ref_seq   # Don't forget the 1st seq!


    ### virtual sequences (seqs that were aligned to the ref seq.)
    count = 0
    for allele, seq in iter_dict_ligated: # the 2nd ~

        tmp_substituted = ''.join([SubstitutionRule(base_v, base_r) for base_v, base_r in zip(seq, ref_seq)]) # update with substituted one.

        # print(
        #     "{}:\n"
        #     "{}\n"
        #     "{}".format(allele, seq, tmp_substituted)
        # )

        dict_RETURN[allele] = tmp_substituted

        count += 1
        # if count >= 300: break


    return dict_RETURN


def splitByVerticalBar(_dict_substituted):

    """
    When 'gen' and 'nuc', the substituted dictionary will be split by '|' before processing Indel.

    Split parts represent functional difference such as '5_prime', 'intron', or 'exon'.

    This is mainly for HLA-DQB1's exon5.

    cDNA              745
    AA codon          217
                      |
    DQB1*05:01:01:01  CTT ATC ATC CGT CAA AGG AGT CGG AAA G|..... ..... ..... ..... ....|GG CTT CTG CAC TGA
    DQB1*05:01:01:02  --- --- --- --- --- --- --- --- --- -|..... ..... ..... ..... ....|-- --- --- --- ---
    DQB1*05:01:01:03  --- --- --- --- --- --- --- --- --- -|..... ..... ..... ..... ....|-- --- --- --- ---

    if the segment itself is either exon or insertion, it will be 'INDEL_{AA,SNPS}_...' in the map file.
    To avoid, split first then process indel afterward with exception handling in the next step.

    """

    dict_RETURN = {}

    ### reference seq. (the 1st seq, practically.)
    ref_allele, ref_seq = next(iter(_dict_substituted.items()))
    # print("ref allele: {}".format(ref_allele))
    # print("ref seq: {}".format(ref_seq))

    iter_vbars = list(re.finditer(r'\|', ref_seq))
    N_vbar = len(iter_vbars)


    """
    0: (Impossible)
    1: 5_prime | 3_prime (Impossible)
    2: 5_prime | exon1 | 3_prime
    3: 5_prime | exon1 | intron1 | 3_prime (Impossible)
    4: 5_prime | exon1 | intron1 | exon2 | 3_prime 
    5: 5_prime | exon1 | intron1 | exon2 | intron2 | 3_prime (Impossible)
    6: 5_prime | exon1 | intron1 | exon2 | intron2 | exon3 | 3_prime

    => # of '|' should be even because an exon must be covered by two '|'s.    
    
    """

    ### Setting annot.

    l_annots = []

    if len(iter_vbars) > 0 and len(iter_vbars) % 2 == 0: # should be even number.

        l_annots = ["exon{}".format(i+1) if j % 2 == 0 else "intron{}".format(i+1) \
                    for i in range(0, int(N_vbar / 2)) for j in range(2)]
        l_annots.pop() # abandon the last intronN. (Last exon must end with '3-prim'.)

        l_annots = ['5_prime'] + l_annots + ['3_prime']

        # print(l_annots)
    else:
        ### Possibly invalid seq. structure to process.
        print("Invalid # of vertical bar('|').")
        return -1

    ### Setting positions
    l_vbar_start = [item.start() for item in iter_vbars]
    l_vbar_end = [item.end() for item in iter_vbars]

    l_str_strat = [0] + l_vbar_end
    l_str_end = l_vbar_start + [len(ref_seq)]
    l_str_pos = list(zip(l_str_strat, l_str_end))
    # print(l_str_pos)
    # print(len(l_str_pos))

    if len(l_str_pos) != len(l_annots):
        print("Length mismatch")
        return -1

    ### split and gather
    for annot, s, e in zip(l_annots, l_str_strat, l_str_end):
        # print("{}: {}-{}".format(annot, s, e))

        dict_temp = {}

        for allele, seq in _dict_substituted.items():
            dict_temp[allele] = seq[s:e]

        dict_RETURN[annot] = dict_temp


    return l_annots, dict_RETURN


def processIndel(_dict_substituted):

    """
    Precisely, this funciton processes INSERTION.
    - r'\.+'(dots) spot in the ref. seq. represents an INSERTION locus.
    - In virtual sequences, r'\.+' spots in other than the above INSERTION spots represent DELITION.

    This function must process only INSERTIONS, so using 're.sub()' which can substitute both INS and DEL is not recommended.
    """

    dict_RETURN ={}

    iter_dict_ligated = iter(_dict_substituted.items())

    ### reference seq. (the 1st seq, practically.)
    ref_allele, ref_seq = next(iter_dict_ligated)
    # print("ref allele: {}".format(ref_allele))
    # print("ref seq: {}".format(ref_seq))

    """
    '.'s in ref seq. are INSERTIONs essentially because those '.'s are not allocated the relative position.
    
    '.'s in virtual seqs. are DELETIONs essentially. 
    """
    iter_indels = list(re.finditer(r'\.+', ref_seq))

    if len(iter_indels) > 0: # when there is at least one indel.

        if len(iter_indels) == 1 and (len(iter_indels[0].group()) == len(ref_seq)): # No need Indel processing.
            # Exception for DQB1's exon5.
            return _dict_substituted


        l_indels_start = [item.start() for item in iter_indels]
        l_indels_end = [item.end() for item in iter_indels]
        l_indels_pos = list(zip(l_indels_start, l_indels_end))

        l_str_start = [0] + l_indels_end
        l_str_end = l_indels_start + [len(ref_seq)]
        l_str_pos = list(zip(l_str_start, l_str_end))


        ### Replacing
        """
        `p_INDEL` is actually to represent "There is no Insertion." if matched.
        
        In virtual seqs, '*' or 'x' will be considered to "No Insertion".
        
        """
        p_INDEL = re.compile(r'\.+|\*+|x+')
        count = 0
        for allele, seq in _dict_substituted.items(): # iterate from the 1st(ref) seq.

            # 'z' or 'Z'. 'Z' means there is INSERTION in virtual sequences.
            l_z_or_Z = ['z' if bool(p_INDEL.match(seq[s : e])) else 'Z' for s, e in l_indels_pos]
            # print(l_z_or_Z)

            l_temp = [seq[s:e] for s, e in l_str_pos]
            l_temp2 = [l_temp[0]]
            l_temp2.extend([item for pair in zip(l_z_or_Z, l_temp[1:]) for item in pair])

            seq_indel_processed = ''.join(l_temp2)

            # print(
            #     "{}:\n"
            #     "{}\n"
            #     "(before): {}\n"
            #     "( after): {}".format(allele, l_z_or_Z, seq, seq_indel_processed)
            # )

            dict_RETURN[allele] = seq_indel_processed

            count += 1
            # if count >= 3000: break

        return dict_RETURN

    else:
        return _dict_substituted


def processIndel_gen(_l_annots, _dict_split):

    dict_RETURN = {annot: None for annot in _l_annots}

    for annot in _l_annots:
        dict_RETURN[annot] = processIndel(_dict_split[annot])

    return dict_RETURN


def join_between_annots(_l_annot, _dict_IndelProcessed):


    return 0



if __name__ == '__main__':

    """
    A module to chop *_{gen,nuc,prot}.txt file.
    """

    pass