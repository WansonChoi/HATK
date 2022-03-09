# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import argparse, textwrap
import numpy as np

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


def HLAtoSequences(_chped, _dictionary, _type, _out_prefix, _f_asLump=False, _f_hasHeader=True,
                   _HLA_target=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

    ### Main Variables ###
    dict_seq = {_HLA_target[i]: {} for i in range(0, len(_HLA_target))}  # Initialization.
    dict_len = {_HLA_target[i]: -1 for i in range(0, len(_HLA_target))}


    ### Main Actions ###

    ##### < [1] Loading HLA Dictionary > #####
    with open(_dictionary, "r") as f_dictionary:

        count = 0

        for line in f_dictionary:

            t_line = re.split(r'\s+', line.rstrip('\n'))
            # ex1) (AA) ['B*58:01:01', 'MLVMAPRTVLLLLSAALALTETWAG...'] (len == 2)
            # ex2) (SNPS) ['B*58:01:01', 'CTAGTCCTGCTTCAGGGTCCGGGGCCCG...'] (len == 2)
            # print(t_line)

            for i in range(0, len(_HLA_target)):

                # if re.match(r'{}\*'.format(_HLA_target[i]), t_line[0]):
                if t_line[0].split('*')[0] == _HLA_target[i]:

                    dict_seq[_HLA_target[i]][t_line[0]] = t_line[1] # Sequence information.

                    if dict_len[_HLA_target[i]] == -1:
                        dict_len[_HLA_target[i]] = len(t_line[1])

                    break  # If a line of given dictionary belongs to either HLA, then checking whether it belongs to other HLAs is useless.

            count += 1
            # if count > 5 : break


    # Result check
    # for i in range(0, len(_HLA_target)):
    #     print("\n{} :".format(_HLA_target[i]))
    #
    #     idx = 0
    #     for k, v in dict_seq[_HLA_target[i]].items():
    #         print("{} : {}".format(k, v))
    #
    #         idx += 1
    #         if idx > 10 : break
    #
    # for k, v in dict_len.items():
    #     print("The length of HLA-{} : {}".format(k, v))

    # The insertion will be precessed in the dictionary generation in advance.
    # print("Insertion check : {}".format(haveInsertion))


    ##### < [2] Transforming each HLA alleles to corresponding sequences > #####

    with open(_out_prefix + ".{}.ped".format(_type), 'w') as f_output:
        f_output.writelines(GenerateLines(_chped, _type, dict_seq, dict_len, _HLA_target, _f_asLump=_f_asLump, _f_hasHeader=_f_hasHeader))



def GenerateLines(_chped, _type, _dict_seq, _dict_len, _HLA_target, _f_asLump=False, _f_hasHeader=True):

    with open(_chped, "r") as f_chped:

        if _f_hasHeader:
            l_header = f_chped.readline().split()
            # print(l_header)

        for line in f_chped:
            # t_line = re.split(r'\s+', line.rstrip('\n'))
            t_line = line.split()
            # print(t_line)

            """
            [0,1,2,3,4,5] := ped file information
            [6,7] := HLA-A,
            [8,9] := HLA-B,
            ...,
            [20, 21] := HLA-DRB1
            (... before update)
            """

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_part__ = '\t'.join([
                BringSequences(t_line[(2 * i + 6)], t_line[(2 * i + 7)], _type, _HLA_target[i],
                               _dict_seq[_HLA_target[i]], _dict_len[_HLA_target[i]], _f_asLump=_f_asLump) for i in range(0, len(_HLA_target))
            ])

            # mem_p2 = process.memory_info().rss / 1024 ** 2
            # print("{}(Mb)".format(mem_p2 - mem_p1))

            yield '\t'.join([__ped_info__, __genomic_part__]) + "\n"



def BringSequences(_HLA_allele1, _HLA_allele2, _type, _hla, _dict_seq, _len, _f_asLump=False):


    ### Dealing with generalized 4-field HLA alleles.

    if _HLA_allele1 in _dict_seq and _HLA_allele2 in _dict_seq:
        Seq1 = _dict_seq[_HLA_allele1]
        Seq2 = _dict_seq[_HLA_allele2]
        # No reversing HLA sequences.
    else:
        Seq1 = list(np.repeat('0', _len))
        Seq2 = list(np.repeat('0', _len))


    if _f_asLump:
        return str([[_HLA_allele1, ''.join(Seq1)], [_HLA_allele2, ''.join(Seq2)]])
    else:
        return '\t'.join(np.array(list(zip(Seq1, Seq2))).flatten())



def BringSequences_old(_HLA_allele1, _HLA_allele2, _type, _hla, _dict_seq, _len, _f_asLump=False):


    ### Dealing with generalized 4-field HLA alleles.

    # `_HLA_allele1`.
    try:
        Seq1 = _dict_seq[_HLA_allele1]
    except KeyError:
        Seq1 = "-1"  # fail

    # `_HLA_allele2`.
    try:
        Seq2 = _dict_seq[_HLA_allele2]
    except KeyError:
        Seq2 = "-1"  # fail

    # No reversing HLA sequences.


    if not _f_asLump:

        ### Main sequence information processing
        if Seq1 != "-1" and Seq2 != "-1":

            # Only when both HLA alleles can get the corresponding HLA sequence information.

            l_temp = []

            for i in range(0, len(Seq1)):
                l_temp.append(Seq1[i])
                l_temp.append(Seq2[i])

            # Reversing
            __return__ = '\t'.join(l_temp)

        else:
            __return__ = '\t'.join(["0" for z in range(0, 2 * _len)])


        return __return__


    else:

        # As each Strings.

        t_Seq1 = ""
        t_Seq2 = ""

        if Seq1 != "-1" and Seq2 != "-1":
            t_Seq1 = Seq1
            t_Seq2 = Seq2
        else:
            t_Seq1 = ''.join(["0" for z in range(0, _len)])
            t_Seq2 = ''.join(["0" for z in range(0, _len)])

        return str([[_HLA_allele1, t_Seq1], [_HLA_allele2, t_Seq2]])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        Original Author: Sherman Jia, 2012

        HLAtoSequences.py
        - This script Converts HLA alleles (in .ped file format) to amino acid or DNA sequences

        Input file should contain: FID, IID, pID, mID, sex, pheno, HLA-A (2), B (2), C (2),
        DPA1 (2), DPB1 (2), DQA1 (2), DQB1 (2), DRB1 (2) ... Broad Order

    #########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"
    # parser._optionals.description = "- Necessary main options.\n"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--chped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("--dict", help="\nHLA dictonary file name(ex. 'HLA_DICTIONARY_AA.txt')\n\n", required=True)
    parser.add_argument("--type", help="\nAA(for Amino Acid) or SNP(for SNPs)\n\n", choices=["AA", "SNPS"], required=True)
    parser.add_argument("--out", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n",
                        default=['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')

    # parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as original version.\n\n",
    #                     action='store_true')
    parser.add_argument("--asLump", help="\n(for Testing) Not zipped result to check strings.\n\n",
                        action='store_true')




    ##### <for Test> #####

    ## (2019. 01. 06.) Introducinig compatibility to work with old version of dictionary

    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference_old/HLA_DICTIONARY_AA.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/HAPMAP_CEU.old.enCODED",
    #                           "--previous-version"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference_old/HLA_DICTIONARY_SNPS.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/HAPMAP_CEU.old.enCODED",
    #                           "--previous-version"])


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    HLAtoSequences(args.chped, args.dict, args.type, args.out, _f_asLump=args.asLump, _HLA_target=args.HLA)
