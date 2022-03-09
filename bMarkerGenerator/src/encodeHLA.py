# -*- coding: utf-8 -*-

# (2017/11/27) recoded by Wanson Choi
import os, re
from os.path import basename, dirname, join
import argparse, textwrap
import numpy as np

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


# Patterns
p_1field = re.compile(r'\w+\*\d{2,3}')


def encodeHLA(_chped, _out_prefix, _hg="18", _f_asSmallLetter=True, _f_get_1field=False, _f_hasHeader=True,
              _HLA_target=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"), _p_data='IMGT2Seq/data'):

    # (2018. 9. 25.) Replaced by lift-over values.
    genepos_hg = LoadGenenomicPosition(_hg, _HLA_target, _p_data)

    ### Acquiring `HLA_allele_sets`.

    HLA_allele_sets = {_HLA_target[i]: [] for i in range(0, len(_HLA_target))}


    with open(_chped, 'r') as f_chped:

        if _f_hasHeader:
            l_header = f_chped.readline().split()
            # print(l_header)

        count = 0

        for l in f_chped:

            """
            l[:6] := ("FID", "IID", "PID", "MID", "Sex", "Phe")
            l[6:8] := HLA-A
            l[8:10] := HLA-B
            ...
            l[20:22] := HLA-DRB1
            """

            t_line = re.split(r'\s+', l.rstrip('\n'))
            # print(t_line)

            for i in range(0, len(_HLA_target)):

                idx1 = 2*i + 6
                idx2 = idx1 + 1

                al1 = t_line[idx1]
                al2 = t_line[idx2]


                # Allele 1
                if al1 != "0":

                    if al1 not in HLA_allele_sets[_HLA_target[i]]:
                        HLA_allele_sets[_HLA_target[i]].append(al1)


                    if _f_get_1field:

                        m = p_1field.match(al1)

                        if m:

                            al1_1field = m.group()

                            if al1_1field not in HLA_allele_sets[_HLA_target[i]]:
                                HLA_allele_sets[_HLA_target[i]].append(al1_1field)


                # Allele 2
                if al2 != "0":

                    if al2 not in HLA_allele_sets[_HLA_target[i]]:
                        HLA_allele_sets[_HLA_target[i]].append(al2)


                    if _f_get_1field:

                        m = p_1field.match(al2)

                        if m:

                            al2_1field = m.group()

                            if al2_1field not in HLA_allele_sets[_HLA_target[i]]:
                                HLA_allele_sets[_HLA_target[i]].append(al2_1field)


            count += 1
            # if count > 5 : break


    for i in range(0, len(_HLA_target)):
        HLA_allele_sets[_HLA_target[i]].sort()


    # # Result checking
    # print("\nHLA alleles.")
    # for k, v in HLA_allele_sets.items():
    #     print("{}: {}".format(k, v))



    ### Making a new *.HLA.map file.

    map_LABELS = ['_'.join(["HLA", HLA_allele_sets[_HLA_target[i]][j]]) for i in range(0, len(_HLA_target)) for j in range(0, len(HLA_allele_sets[_HLA_target[i]]))]
    # print(map_LABELS)

    map_POS = [str(genepos_hg[_HLA_target[i]]) for i in range(0, len(_HLA_target)) for _ in range(0, len(HLA_allele_sets[_HLA_target[i]]))]
    # print(map_POS)

    with open(_out_prefix + ".map", 'w') as f_HLA_map:
        f_HLA_map.writelines(('\t'.join(["6", map_LABELS[i], "0", map_POS[i]]) + "\n" for i in range(0, len(map_LABELS))))


    ### Making a new *.HLA.ped file.

    with open(_out_prefix + ".ped", 'w') as f_HLA_ped:
        f_HLA_ped.writelines(MakeHLAPed(_chped, HLA_allele_sets, _HLA_target, _f_asSmallLetter=_f_asSmallLetter, _f_hasHeader=_f_hasHeader))



    return _out_prefix + ".ped", _out_prefix + ".map"



def MakeHLAPed(_CHPED, _HLA_allele_sets, _HLA_target, _f_asSmallLetter=False, _f_hasHeader=True):

    with open(_CHPED, 'r') as f_chped:

        if _f_hasHeader:
            l_header = f_chped.readline().split()
            # print(l_header)

        count = 0

        for l in f_chped:

            t_line = re.split(r'\s+', l.rstrip('\n'))

            """
            t_line[:6] := ("FID", "IID", "PID", "MID", "Sex", "Phe")
            t_line[6:8] := HLA-A
            t_line[8:10] := HLA-B
            ...
            t_line[20:22] := HLA-DRB1
            """

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                PrintGenotypes4(t_line[2 * i + 6], t_line[2 * i + 7], _HLA_allele_sets[_HLA_target[i]],
                                _f_asSmallLetter=_f_asSmallLetter)
                for i in range(0, len(_HLA_target)) if len(_HLA_allele_sets[_HLA_target[i]]) > 0
            ])

            __return__ = '\t'.join([__ped_info__, __genomic_info__])


            yield __return__ + "\n"

            count += 1



# (2019. 1. 3.) Introduced for memory issues.
def PrintGenotypes4(_allele1, _allele2, _HLA_allele_sets, _f_asSmallLetter=False):

    l_output = []

    _present_ = "p" if _f_asSmallLetter else "P"
    _absent_ = "a" if _f_asSmallLetter else "A"



    for i in range(0, len(_HLA_allele_sets)):

        _ALLELE = _HLA_allele_sets[i]

        G1 = "-1"
        G2 = "-1"

        if _allele1 != "0" and _allele2 != "0":

            # Dealing with generalized 4-field HLA alleles.

            p_1field = re.compile(r'\w+\*\d{2,3}')


            # Allele 1
            m = p_1field.match(_allele1)
            _al1_1field = m.group() # Just assume that given alleles is definitely in the form of r'\w+\*(\d{2,3}:?)+'

            # Determining 'Present' or 'Absent'.
            if _allele1 == _ALLELE or _al1_1field == _ALLELE:
                G1 = _present_
            else:
                G1 = _absent_


            # Allele 2
            m = p_1field.match(_allele2)
            _al2_1field = m.group()
            if _allele2 == _ALLELE or _al2_1field == _ALLELE:
                G2 = _present_
            else:
                G2 = _absent_


        else:
            # If at least one HLA allele which is given in *.chped is "0", then consider both of them are "0"
            G1 = "0"
            G2 = "0"




        if G1 == "0" or G2 == "0":
            l_output.append("0")
            l_output.append("0")
        else:
            l_output.append(G1)
            l_output.append(G2)


    return '\t'.join(l_output)



def LoadGenenomicPosition(_hg, _HLA_target, _p_data):

    HLA_EXON1_START_POS = join(_p_data, "HLA_EXON1_START_CODON_POSITIONS_hg{}.txt".format(_hg))

    dict_HLA_EXON1_START_POS = {}

    with open(HLA_EXON1_START_POS, 'r') as f_pos:
        for line in f_pos:

            [hla, direction, exon1_pos] = line.split()

            if exon1_pos != 'NA':
                dict_HLA_EXON1_START_POS[hla] = int(exon1_pos)


    if not np.isin(_HLA_target, dict_HLA_EXON1_START_POS).all():
        for hla in _HLA_target:
            if hla not in dict_HLA_EXON1_START_POS:
                raise ValueError(std_ERROR_MAIN_PROCESS_NAME + "HLA-{}'s BP information is NA in '{}'.".format(hla, HLA_EXON1_START_POS))

    return dict_HLA_EXON1_START_POS



if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        encodeHLA.py

        This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
        The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
                                            2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

    ###########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--chped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("--out", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="hg", default="18")

    parser.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n",
                        default=['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')

    # parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work like original version.\n\n",
    #                     action='store_true')
    parser.add_argument("--asSmallLetter", help="\n'P'resent and 'A'bsent to 'p'resent and 'a'bsent.\n\n",
                        action='store_true')
    # parser.add_argument("--addDummyMarker", help="\nAdd dummy marker to prevent the glitch in work with plink(1.07).\n\n",
    #                     action='store_true')



    ##### <for Test> #####

    # (2019. 01. 06.)
    # # --previous-version
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-hg", "18",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190107/_3_encodeHLA/_Case_HAPMAP_CEU.HLA",
    #                           "--previous-version"])

    # # Generalized 4-field HLA alleles
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "18",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190107/_3_encodeHLA_4field/_Case_HAPMAP_CEU.HLA_4fieldTest",
    #                           "--asSmallLetter",
    #                           "--addDummyMarker"])





    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    encodeHLA(args.chped, args.out, args.hg, _f_asSmallLetter=(not args.asSmallLetter), _HLA_target=args.HLA,
              _p_data="IMGT2Seq/data")

    # LoadGenenomicPosition("18", ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"), "/home/wansonchoi/sf_VirtualBox_Share/HATK/IMGT2Seq/data/")