# -*- coding: utf-8 -*-

import os, sys, re
import subprocess
import pandas as pd
import argparse, textwrap


def main(_nuc, _gen = "Not_given", _prot = "Not_given", _out = "Not_given"):

    ### Test Stage

    df_raw_Markers_nuc = LoadRawSeq2(_nuc, "nuc")
    df_raw_Markers_nuc.to_csv(_out+".nuc.raw.markers.txt", sep='\t', header=False, index=True)
    print(df_raw_Markers_nuc.head())

    df_raw_Markers_gen = LoadRawSeq2(_gen, "gen")
    df_raw_Markers_nuc.to_csv(_out+".gen.raw.markers.txt", sep='\t', header=False, index=True)
    print(df_raw_Markers_gen.head())

    df_raw_Markers_prot = LoadRawSeq2(_prot, "prot")
    df_raw_Markers_nuc.to_csv(_out+".prot.raw.markers.txt", sep='\t', header=False, index=True)
    print(df_raw_Markers_prot.head())


    return 0


# def LoadRawSeq(_nuc)

def LoadRawSeq2(_nuc_filename,  _type):

    if not (_type == "nuc" or _type == "gen" or _type == "prot"):
        print("\nArugment \"_type\" must be either \"nuc\", \"gen\" or \"prot\"")
        return -1

    print("\n\nLoading \"{0}\" file".format(_nuc_filename))

    ##### < Core Local Variables > #####

    line_number_marks = None
    Sequence_Chunks = []

    f = None
    k_first = None

    dict_first = None


    ##### MAKRDING_CHUNKS

    # Classification based on file type
    if _type == "nuc":
        Mark_Char = "cDNA"
    elif _type == "gen":
        Mark_Char = "gDNA"
    elif _type == "prot":
        Mark_Char = "Prot"
    else:
        print("\nArugment \"_type\" must be either \"nuc\", \"gen\" or \"prot\"")
        return -1


    line_number_marks = subprocess.check_output(' '.join(["awk", "'{print $1}'", _nuc_filename, "|", "grep -n", Mark_Char]), shell=True)
    line_number_marks = line_number_marks.decode('UTF-8').rstrip('\n').split('\n')
    line_number_marks = [int(item.split(':')[0]) - 1 for item in line_number_marks]

    print("\nMarks is {0}".format(line_number_marks))

    # Classification based on file type 2
    if _type == "nuc":
        allele_start = line_number_marks[0] + 3
        allele_end = line_number_marks[1] - 2
    elif _type == "gen":
        allele_start = line_number_marks[0] + 2
        allele_end = line_number_marks[1] - 3
    elif _type == "prot":
        allele_start = line_number_marks[0] + 2
        allele_end = line_number_marks[1] - 2


    NumberOfAlleles = allele_end - (allele_start - 1)

    print("\nallele start and end are {0} and {1}\n".format(allele_start, allele_end))
    print("# of Alleles is {0}\n".format(NumberOfAlleles))


    ##### READING_RAW_SEQUENCES

    f = open(_nuc_filename)
    f_lines = f.readlines()
    f.close()

    # print(f_lines[-2]) # "Please see http://hla.alleles.org/terms.html for terms of use."
    line_number_marks.append(len(f_lines) - 2)


    ##### APPENDING_RAW_SEQUENCES

    # Preparing HLA allele names as a set of keys.
    k_first = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(allele_start, allele_end + 1)]
    # print("k_first are {0}".format(k_first))

    dict_first = {}

    for i in range(0, len(k_first)):
        s = re.split('\s+', f_lines[allele_start + i].lstrip(' ').rstrip(' \n'))[1:]
        dict_first[k_first[i]] = ''.join(s)

    # for k,v in dict_first.items():
    #     print("{0}: {1}".format(k, v))

    Sequence_Chunks.append(dict_first)

    ### Processing each Chunks

    for i in range(1, len(line_number_marks) - 1):
        t_keys = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(line_number_marks[i] + 3, line_number_marks[i + 1] - 2 + 1)]
        t_s = [''.join(re.split('\s+', f_lines[(line_number_marks[i] + 3) + k].lstrip(' ').rstrip(' \n'))[1:]) for k in range(0, len(t_keys))]

        Sequence_Chunks.append({k: v for k, v in zip(t_keys, t_s)})


    ### Merging respective Dictionaries

    for i in range(0, len(k_first)):
        t_string = ''.join([Sequence_Chunks[j][k_first[i]] for j in range(0, len(Sequence_Chunks)) if k_first[i] in Sequence_Chunks[j].keys()])
        dict_first[k_first[i]] = tuple(t_string)

    # Unprocessed DataFrame of codon information.
    df_CODON = pd.DataFrame.from_dict(dict_first, orient='index').fillna('x')





    # (2)

    conserved_1st_seq = df_CODON.iloc[0, :]
    # print("\nconserved_1st_seq\n")
    # print(conserved_1st_seq)

    # fordrop = conserved_1st_seq[conserved_1st_seq == '.'].index.tolist()
    #
    # df_CODON_trimmed = df_CODON.drop(fordrop, axis=1)
    # df_CODON_trimmed.columns = pd.Index(range(0, df_CODON_trimmed.shape[1]))
    # DataFrame where insertions('....') are removed.


    for i in range(0, len(conserved_1st_seq)):

        if conserved_1st_seq[i] != '.':

            col = df_CODON.iloc[1:, i]
            col[col == '-'] = conserved_1st_seq[i]
            col[col == '*'] = 'x'


    return df_CODON

    ##### < Trimming insertion > #####

    # In "*_nuc.txt" file, '....' and '|' should be removed.
    # After splitting by '|', the number of splitted parts must match the number of exon parts.

    """
    (1) Starting as tuples.
    (2) Removing insertions seqs in `conserved_1st_seq`.
    (3) Making a single string by aggregating  columns and splitting it by '|'.
    (4) Checking genomic_position information by exons.
    (5) Assuming the length is matched, make a header per base positions.
    (6) Finally, generate DataFrame by 3-base(codon).
    """


    # # (3)
    #
    # # DataFrame of condon sequence.
    # df_CODON_seq = df_CODON_trimmed.apply(lambda x: ''.join(x.astype(str)), axis=1)
    # df_CODON_seq.to_csv(_OUTPUT + '.codon.seq.txt', sep='\t', header=None, index=True)
    #
    # dict_CODON_seq = df_CODON_seq.apply(lambda x: x.split('|')).to_dict()
    #
    # # (4)
    #
    # # DataFrame splitted by exons.
    # df_CODON_seq_splited = pd.DataFrame.from_dict(dict_CODON_seq, orient='index')
    # df_CODON_seq_splited.to_csv(_OUTPUT + '.codon.exon.txt', sep='\t', header=False, index=True)
    #
    # # (5)
    #
    # # Genomic Positions information will be allocated to `df_CODON_seq_splitted` using "HLA_INTEGRATED_POSITIONS_hg.txt".
    #
    # # if len(name_check) == 1:
    # #     name_check = name_check.iat[-1]
    # # else:
    # #     print("HLA name got wrong")
    # #     sys.exit()
    #
    # print("\nCurrent file's HLA type is %s" % HLA_name)
    #
    # HLA_EXON_POS = HLA_INTEGRATED_POSITIONS_hg.loc[HLA_name].filter(regex='exon\d', axis=0)
    #
    # ### ***Most important condition***
    # if df_CODON_seq_splited.shape[1] == HLA_EXON_POS.shape[0]:
    #
    #     # This condition is the reason why "*_nuc.txt" file is used.
    #     # It is checking whether the number of exons is matched to that of the parts created by splitting sequence with '|'.
    #
    #     df_CODON_seq_trimmed = df_CODON_seq_splited.apply(lambda x: ''.join(x.astype(str)), axis=1)
    #     # df_CODON_seq_trimmed.to_csv(_OUTPUT+'.codon.seq.trim.txt', sep='\t', header=False, index=True)
    #     df_CODON_marker_trimmed = pd.concat(
    #         [pd.DataFrame([tuple(df_CODON_seq_trimmed.iat[i]) for i in range(0, df_CODON_seq_trimmed.shape[0])])],
    #         axis=1)
    #
    #     genuine_contents = re.findall('\w+[^x]', df_CODON_seq_trimmed.iat[0])
    #
    #     if len(genuine_contents) == 1:
    #
    #         # The length of `genuine_contents` must be 1.
    #         genuine_contents = genuine_contents[0]
    #
    #         if len(genuine_contents) % 3 != 0:
    #             # Strange cases.
    #             print('string of nuc file is wrong !')
    #             sys.exit()
    #
    #     else:
    #         print('string of nuc file is wrong !')
    #         sys.exit()
    #
    #     genuine_contents_codon = [''.join([genuine_contents[i], genuine_contents[i + 1], genuine_contents[i + 2]]) for i
    #                               in range(0, len(genuine_contents), 3)]
    #
    #     ### Extracting genomic positions corresponding to exon parts in "HLA_INTEGRATED_POSITIONS_hg.txt"
    #     l_header = []
    #
    #     if isREVERSE[HLA_name]:
    #
    #         for i in range(0, HLA_EXON_POS.shape[0]):
    #             l_header.extend([z for z in range(HLA_EXON_POS.iat[i, 0], HLA_EXON_POS.iat[i, 1] - 1, -1)])
    #
    #     else:
    #
    #         for i in range(0, HLA_EXON_POS.shape[0]):
    #             l_header.extend([z for z in range(HLA_EXON_POS.iat[i, 0], HLA_EXON_POS.iat[i, 1] + 1)])
    #
    #     genomic_position_by3 = [l_header[i] for i in range(1, len(l_header), 3)]
    #
    #     # AA_HEADER = pd.Series(genuine_contents_codon, index=genomic_position_by3)
    #     # AA_HEADER.to_csv(_OUTPUT+'.HEADER.txt', sep='\t', header=False, index=True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        LoadRawSeq.py


    ###########################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-nuc", help="\nInput *_nuc.txt file.\n\n", required=True)
    parser.add_argument("-gen", help="\nInput *_gen.txt file.\n\n")
    parser.add_argument("-prot", help="\nInput *_prot.txt file.\n\n")
    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)


    # for Publish
    # args = parser.parse_args()

    # for Testing
    args = parser.parse_args(["-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
                              "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
                              "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt",
                              "-o", "NEW_Dict"])

    print(args)

    main(args.nuc, args.gen, args.prot, args.o)