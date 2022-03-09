# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join
import argparse, textwrap

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


def encodeVariants(_ped, _map, _out_prefix, _f_asSmallLetter=True):

    # Processing line by line with python Generators to save memory

    ########## < Control Flags > ##########

    _1_ALLELE_OVERLAPPING = 1
    _2_MAKING_NEW_PEDFILE = 1
    _3_MAKING_NEW_MAPFILE = 1
    _4_MAKING_ALLELELIST = 1


    if _1_ALLELE_OVERLAPPING:

        ########## < [1] Allele overlapping > ##########

        # Acquiring column number
        n_loci = 0
        n_row = 0

        with open(_ped, 'r') as f_ped:

            for line in f_ped:

                t_line = re.split(r'\s+', line.rstrip('\n'))
                # print(t_line[6:])

                if n_row == 0:

                    genomic_info = t_line[6:]
                    n_loci = int(len(genomic_info)/2) # == len(l_factors)

                    # Initializing the list containing factors which appear in each locus.
                    l_factors = [[] for i in range(0, n_loci)] # Initialization


                for i in range(0, n_loci):

                    idx1 = 2*i + 6 # index for `t_line`
                    idx2 = idx1 + 1

                    ##### Allele overlapping
                    if t_line[idx1] != "0" and (t_line[idx1] not in l_factors[i]):
                        l_factors[i].append(t_line[idx1])
                    if t_line[idx2] != "0" and (t_line[idx2] not in l_factors[i]):
                        l_factors[i].append(t_line[idx2])


                n_row += 1
                # if n_row > 5 : break


        # Sorting elements of each lists.
        for i in range(0, len(l_factors)):
            l_factors[i].sort()

        ### --- `l_factors` done.




    if _2_MAKING_NEW_PEDFILE:

        ########## < [2] Making new .ped file > ##########

        with open(_out_prefix + ".ped", 'w') as f_NewPed:
            f_NewPed.writelines(MakeNewPed(_ped, l_factors, _f_asSmallLetter))



    if _3_MAKING_NEW_MAPFILE:

        ########## < [3] Making new .map file > ##########

        with open(_out_prefix + ".map", 'w') as f_NewMap:
            f_NewMap.writelines(MakeNewMap(_map, l_factors))



    if _4_MAKING_ALLELELIST:

        ########## < [4] Making *.allelelist file > ##########

        with open(_out_prefix + ".factors", 'w') as f_allelelist:
            f_allelelist.writelines(MakeAlleleList(_map, l_factors))



def MakeNewPed(_ped, _l_factors, _f_asSmallLetter=True):

    count = 0

    with open(_ped, 'r') as f_ped:
        for line in f_ped:
            t_line = re.split(r'\s+', line.rstrip('\n'))

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                divideToBinaryMarkers(t_line[2 * i + 6], t_line[2 * i + 7], _l_factors[i], _f_asSmallLetter) for i in range(0, len(_l_factors))
            ])

            __return__ = '\t'.join([__ped_info__, __genomic_info__])

            yield __return__ + "\n"


def divideToBinaryMarkers(_SNP1, _SNP2, _factors, _f_asSmallLetter=True):

    _present_ = "p" if _f_asSmallLetter else "P"
    _absent_ = "a" if _f_asSmallLetter else "A"

    Seq = []

    if len(_factors) > 2:

        for j in range(0, len(_factors)):

            if _SNP1 == "0" or _SNP2 == "0":
                Seq.append("0"); Seq.append("0")

            else:

                if _factors[j] == _SNP1:
                    Seq.append(_present_)
                else:
                    Seq.append(_absent_)

                if _factors[j] == _SNP2:
                    Seq.append(_present_)
                else:
                    Seq.append(_absent_)

        # if len(_factors) > 3:
        #
        #     j_end = 1 if len(_factors) == 4 else len(_factors)
        #
        #     for j in range(0, j_end):
        #
        #         for k in range(j + 1, len(_factors)):
        #
        #             if _SNP1 == "0" or _SNP2 == "0":
        #                 Seq.append("0"); Seq.append("0")
        #
        #             else:
        #                 if _factors[j] == _SNP1 or _factors[k] == _SNP1:
        #                     Seq.append(_present_)
        #                 else:
        #                     Seq.append(_absent_)
        #
        #                 if _factors[j] == _SNP2 or _factors[k] == _SNP2:
        #                     Seq.append(_present_)
        #                 else:
        #                     Seq.append(_absent_)
        #
        #     if len(_factors) > 5:
        #
        #         j_end = 1 if len(_factors) == 6 else len(_factors)
        #
        #         for j in range(0, j_end):
        #             for k in range(j + 1, len(_factors)):
        #                 for l in range(k + 1, len(_factors)):
        #
        #                     if _SNP1 == "0" or _SNP2 == "0":
        #                         Seq.append("0"); Seq.append("0")
        #
        #                     else:
        #                         if _factors[j] == _SNP1 or _factors[k] == _SNP1 or _factors[l] == _SNP1:
        #                             Seq.append(_present_)
        #                         else:
        #                             Seq.append(_absent_)
        #
        #                         if _factors[j] == _SNP2 or _factors[k] == _SNP2 or _factors[l] == _SNP2:
        #                             Seq.append(_present_)
        #                         else:
        #                             Seq.append(_absent_)



    else:
        # Most of Cases have length less than equal 2, they will fall into this if-else block.
        if _SNP1 == "0" or _SNP2 == "0":
            Seq.append("0"); Seq.append("0")
        else:
            Seq.append(_SNP1); Seq.append(_SNP2)

    return '\t'.join(Seq)


def MakeNewMap(_p_map, _l_factors):

    count = 0

    with open(_p_map, 'r') as f_map:

        for l in f_map:

            idx = count # index for `l_factors`

            if len(_l_factors[idx]) > 2:

                t_line = re.split(r'\s+', l.rstrip('\n'))

                for j in range(0, len(_l_factors[idx])):
                    yield '\t'.join([t_line[0], t_line[1] + '_' + _l_factors[idx][j], t_line[2], t_line[3]]) + "\n"


                # if len(_l_factors[idx]) > 3:
                #
                #     j_end = 1 if len(_l_factors[idx]) == 4 else len(_l_factors[idx])
                #
                #     for j in range(0, j_end):
                #         for k in range(j+1, len(_l_factors[idx])):
                #             yield '\t'.join([t_line[0], t_line[1]+ '_' + _l_factors[idx][j]+_l_factors[idx][k], t_line[2], t_line[3]]) + "\n"
                #
                #
                #     if len(_l_factors[idx]) > 5:
                #
                #         j_end = 1 if len(_l_factors[idx]) == 6 else len(_l_factors[idx])
                #
                #         for j in range(0, j_end):
                #             for k in range(j+1, len(_l_factors[idx])):
                #                 for l in range(k+1, len(_l_factors[idx])):
                #                     yield '\t'.join([t_line[0], t_line[1]+ '_' + _l_factors[idx][j]+_l_factors[idx][k]+_l_factors[idx][l], t_line[2], t_line[3]]) + "\n"

            else:
                yield l # "\n" is included.

            count += 1



def MakeAlleleList(_p_map, _l_factors):

    count = 0

    with open(_p_map, 'r') as f_map:

        for l in f_map:

            idx = count

            t_line = re.split(r'\s+', l.rstrip('\n'))

            locus_label = t_line[1]
            alleleset = _l_factors[idx]

            yield '\t'.join([locus_label, str(alleleset)]) + "\n"

            count += 1




if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        Original Author : Sherman Jia, 2012
     
        encodeVariants.py
     
        - This script generates PLINK binary markers which encodes multi-alleles(factors) at 
            a locus (Amino acids and SNPs) 

    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--ped", help="\nThe *.ped file which is generated by 'HLAtoSequences.py'.\n\n", required=True)
    parser.add_argument("--map", help="\nThe *.map file which is generated by 'HLAtoSequences.py'.\n\n", required=True)
    parser.add_argument("--out", help="\nOutput file prefix.\n\n", required=True)


    ##### <for Test> #####

    # (2018. 7. 13.)
    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.enCODED"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.enCODED"])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    encodeVariants(args.ped, args.map, args.o)