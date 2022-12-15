# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

import IMGT2Seq.src.processModules as modules
from src.util import printDict


def processGen(_gen, _hla, _BP_start:int, _direction:str, _ref_allele:str,
               _f_prime_seq=False):

    ### Load gen raw seqs
    dict_chunks_gen = modules.chopAlignmentFile(_gen)

    ### Join seqs within chunks
    l_joined_within = modules.join_within_chunks_Prot(dict_chunks_gen)

    ### Join seqs between chunks
    dict_joined_between = modules.join_between_chunks_Prot(l_joined_within)
    # printDict(dict_joined_between, 10)

    ### Ligate the right end.
    dict_ligated = modules.ligateRightEnd(dict_joined_between)
    # printDict(dict_ligated, 10)

    ### Substitute bases.
    dict_substituted = modules.substituteBase(dict_ligated)
    # printDict(dict_substituted, 10)

    ### split by the vertical bar('|'). (for functional annot.)
    l_annots, dict_split = modules.splitByVerticalBar(dict_substituted)
    # for annot in l_annots:
    #     print("{}:".format(annot))
    #     printDict(dict_split[annot], 10)

    ### process Indel for each annotated seqs.
    dict_Indel_Processed = modules.processIndel_gen(l_annots, dict_split)
    # for annot in l_annots:
    #     print("[{}]:".format(annot))
    #     print("(before):")
    #     printDict(dict_split[annot], 10)
    #     print("(after):")
    #     printDict(dict_Indel_Processed[annot], 10)


    ### (Map) allocate BP
    _ref_allele = getRefSeq_gen(_hla, _ref_allele, dict_Indel_Processed)
    dict_ref_seq = {annot: dict_Indel_Processed[annot][_ref_allele] for annot in l_annots}
    # printDict(dict_ref_seq)

    df_Map0_SNPS = allocateBP_gen(dict_ref_seq, _BP_start, _direction, _hla) # will be used to make AA Map0 file.
    print(df_Map0_SNPS)

    # export for testing
    df_Map0_SNPS.to_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_SNPS.txt", sep='\t', header=True, index=False)


    ### Join betw. annots. (Final step for gen seq dictionary)
    dict_RETURN = modules.join_between_annots(l_annots, dict_Indel_Processed)
    printDict(dict_RETURN, 20)

    return dict_RETURN, df_Map0_SNPS


def getRefSeq_gen(_hla, _ref_allele, _dict_Indel_Processed):

    _ref_allele = "{}*{}".format(_hla, _ref_allele)

    _dict_5prime = next(iter(_dict_Indel_Processed.values())) # `_dict_Indel_Processed` of gen is a dictionary of dictionaries.

    if _ref_allele in _dict_5prime.keys():
        return _ref_allele
    else:
        return list(_dict_5prime.keys())[0]



def allocateBP_gen(_dict_ref_seq, _BP_start:int, _direction:str, _hla):

    l_rel_pos = []
    l_gen_pos = []
    l_annot = []

    iter_dict_ref_seq = iter(_dict_ref_seq.items())

    ### 5_prime part.
    annot, refseq = next(iter_dict_ref_seq)

    temp_rel_pos = -1
    temp_gen_pos = _BP_start + (-1 if _direction == '+' else 1)

    for base in reversed(refseq): # easier to process 5_prime seq when reversed.
        if base == 'z' or base == 'Z':

            # last rel_pos and gen_pos
            last_rel_pos = temp_rel_pos + 1
            last_gen_pos = temp_gen_pos + (1 if _direction == '+' else -1)

            INS_rel_pos = "{}x{}".format(last_rel_pos, temp_rel_pos)
            INS_gen_pos = int(round((last_gen_pos + temp_gen_pos) / 2))
            INS_label = "INS_SNPS_{HLA}_{REL_POS}_{GEN_POS}" \
                            .format(HLA=_hla, REL_POS=INS_rel_pos, GEN_POS=INS_gen_pos)

            l_rel_pos.append(INS_rel_pos)
            l_gen_pos.append(INS_gen_pos)
            l_annot.append(annot)

        else:

            label = "SNPS_{HLA}_{REL_POS}_{GEN_POS}_{ANNOT}" \
                    .format(HLA=_hla, REL_POS=temp_rel_pos, GEN_POS=temp_gen_pos, ANNOT=annot)

            l_rel_pos.append(temp_rel_pos)
            l_gen_pos.append(temp_gen_pos)
            l_annot.append(annot)


            temp_rel_pos += -1
            temp_gen_pos += (-1 if _direction == '+' else 1)

    l_rel_pos = list(reversed(l_rel_pos))
    l_gen_pos = list(reversed(l_gen_pos))
    l_annot = list(reversed(l_annot))


    ### Rest parts. (exon1 ~ 3_prime)

    temp_rel_pos = 1
    temp_gen_pos = _BP_start

    for annot, refseq in iter_dict_ref_seq:
        for base in refseq:
            if base == 'z' or base == 'Z':
                # l_rel_pos.append('i')
                # l_gen_pos.append('i')
                # l_annot.append('i')

                # last rel_pos and gen_pos
                last_rel_pos = temp_rel_pos - 1
                last_gen_pos = temp_gen_pos + (-1 if _direction == '+' else 1)

                INS_rel_pos = "{}x{}".format(last_rel_pos, temp_rel_pos)
                INS_gen_pos = int(round((last_gen_pos + temp_gen_pos) / 2))
                INS_label = "INS_SNPS_{HLA}_{REL_POS}_{GEN_POS}" \
                    .format(HLA=_hla, REL_POS=INS_rel_pos, GEN_POS=INS_gen_pos)

                l_rel_pos.append("{}x{}".format(last_rel_pos, temp_rel_pos))
                l_gen_pos.append(int(round((last_gen_pos + temp_gen_pos) / 2)))
                l_annot.append(annot)

            else:

                label = "SNPS_{HLA}_{REL_POS}_{GEN_POS}_{ANNOT}" \
                    .format(HLA=_hla, REL_POS=temp_rel_pos, GEN_POS=temp_gen_pos, ANNOT=annot)

                l_rel_pos.append(temp_rel_pos)
                l_gen_pos.append(temp_gen_pos)
                l_annot.append(annot)


                temp_rel_pos += 1
                temp_gen_pos += (1 if _direction == '+' else -1)


    # mid check
    # count = 0
    # for label, rel_pos, gen_pos, annot, base in zip(l_SNPS_labels, l_rel_pos, l_gen_pos, l_annot, ''.join(list(_dict_ref_seq.values()))):
    #     print("{} {} {} {} {}".format(rel_pos, gen_pos, annot, base, label))
    #
    #     count += 1
    #     # if count >= 5 : break


    df_RETURN = pd.DataFrame([
        [rel_pos, gen_pos, annot, base] for rel_pos, gen_pos, annot, base in  \
                zip(l_rel_pos, l_gen_pos, l_annot, ''.join(list(_dict_ref_seq.values())))
    ], columns=['SNPS_rel_pos', 'SNPS_gen_pos', 'Type', 'base'])

    df_RETURN.to_csv("/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/df_Map0_SNPS.txt", sep='\t', header=True, index=False)

    return df_RETURN



def export_dict_refseq(_dict_ref_seq, _out): # a function for manual testing.

    with open(_out, 'w') as f_out:
        for annot, seq in _dict_ref_seq.items():
            f_out.write("{}\t{}\n".format(annot, seq))

    return _out



if __name__ == '__main__':

    # A, gen
    # gen = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3500/alignments/A_gen.txt"
    # gen = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3500/alignments/A_gen.txt"

    # DQB1, gen
    gen = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3500/alignments/DQB1_gen.txt"
    # gen = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3500/alignments/DQB1_gen.txt"

    r = processGen(gen, 'DQB1', 32742362, '-', '05:01:01:01')


    # dict_ref_seq = {}
    # with open("/Users/wansonchoi/Git_Projects/HATK/tests/refseq_DQB1.txt", 'r') as f_refseq:
    #     for line in f_refseq:
    #         l = line.split()
    #         dict_ref_seq[l[0]] = l[1]
    #
    # # printDict(dict_ref_seq)
    # allocateBP_gen(dict_ref_seq, 32742362, '-', "DQB1")


    """
    (Gen)
    - Load raw sequences.
    - Get relative position start in the 1st chunk.

    간단하게 생각해서, finditer를 쓰려면 전에 ligateRightEnd를 해야함.
    만약 '|'로 partition을 한 seq들에 대해 prot처럼 처리하는걸 반복한다고 하더라도 ligateRightEnd는 앞서서 해야한다는 거지.
    만약 '|'로 partition을 한다고 하더라도 ref seq상의 '|'만 있으면 되긴 함. (과거에 df처럼 처리하던 것과 달리 ref seq에만 finditer하면 돼서.)
        - 다시 말해 substitute할때도 '|'이어야 할 부분을 'x'로 채워놔도, ref seq상의 '|'의 position으로 partition하면 노상관.

    사실 이런식면 ligateRightEnd랑 substituteBase까지 똑같이하고 ProcessIndel만 '|'로 잘라서 하면 되는거기도 함.



    """

    """
    (One exception to consider)
    seq1: ACGT|TCGGAA
    seq2: (blank)

    After ligation

    seq1: ACGT'|'TCGGAA
    seq2: xxxx'x'xxxxxx

    => 're.sub' function shouldn't be used.

    (deprecated)
    """