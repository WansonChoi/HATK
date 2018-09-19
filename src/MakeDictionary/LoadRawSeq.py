# -*- coding: utf-8 -*-

import os, sys, re
import subprocess
import pandas as pd
import argparse, textwrap


def main(_out, _hla, _type, _hg_Table, _imgt,
         _nuc, _gen = "Not_given", _prot = "Not_given"):

    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}


    """
    The list of DataFrames to return.
    
    1. df_raw_Makrers
    2. df_raw_Seqs
    3. df_raw_Seqs_splitted (for checking Exon region)
    """


    ### Test Stage

    print("\n=============<NUC>=============\n")

    # Raw Markers
    df_raw_Markers_nuc, start_offset_nuc = LoadRawSeq2(_nuc, "A", "nuc")
    df_raw_Markers_nuc.to_csv(_out+".HLA_{0}.nuc.raw.markers.txt".format("A"), sep='\t', header=False, index=True)
    print(df_raw_Markers_nuc.head())
    print(df_raw_Markers_nuc.tail())

    # Raw Seqeunces
    df_raw_Seqs_nuc = get_DF_rawSeqs(df_raw_Markers_nuc)
    df_raw_Seqs_nuc.to_csv(_out+".HLA_{0}.nuc.raw.seqs.txt".format("A"), sep='\t', header=False, index=True)
    print("\nRaw seqs of nuc\n")
    print(df_raw_Seqs_nuc.head())
    print(df_raw_Seqs_nuc.tail())

    # Raw Seqeunces splitted (for Exon checking)
    df_raw_Seqs_splitted_nuc = get_DF_rawSeqsSplitted(df_raw_Seqs_nuc)
    df_raw_Seqs_splitted_nuc.columns = pd.Index(["Exon"+str(i+1) for i in range(0, df_raw_Seqs_splitted_nuc.shape[1])])
    df_raw_Seqs_splitted_nuc.index.name = "Alleles"
    df_raw_Seqs_splitted_nuc.to_csv(_out+".HLA_{0}.nuc.raw.seqs.splitted.txt".format("A"), sep='\t', header=True, index=True)
    print("\nRaw seqs of nuc splitted.\n")
    print(df_raw_Seqs_splitted_nuc.head())

    # (Trimmed) Sequences splitted without any indel character(for Exon checking with gen file)
    df_Seqs_splited_noIndel_nuc = df_raw_Seqs_splitted_nuc.apply(lambda x : ProcessIndel(x, _remove_indel=True), axis = 0)
    print(df_Seqs_splited_noIndel_nuc.head())
    df_Seqs_splited_noIndel_nuc.to_csv(_out+".HLA_{0}.nuc.seqs.NoIndelChar.splitted.txt".format("A"), sep='\t', header=True, index=True)

    # No more processing to *_nuc.txt file is needed.


    print("\n=============<GEN>=============\n")

    # Raw Markers
    df_raw_Markers_gen, start_offset_gen = LoadRawSeq2(_gen, "A", "gen")
    df_raw_Markers_gen.to_csv(_out+".HLA_{0}.gen.raw.markers.txt".format("A"), sep='\t', header=False, index=True)
    print(df_raw_Markers_gen.head())
    print(df_raw_Markers_gen.tail())

    # Raw Sequences
    df_raw_Seqs_gen = get_DF_rawSeqs(df_raw_Markers_gen)
    df_raw_Seqs_gen.to_csv(_out+".HLA_{0}.gen.raw.seqs.txt".format("A"), sep='\t', header=False, index=True)
    print("\nRaw seqs of gen\n")
    print(df_raw_Seqs_gen.head())
    print(df_raw_Seqs_gen.tail())

    # Raw Seqeunces splitted (for Exon checking)
    df_raw_Seqs_splitted_gen = get_DF_rawSeqsSplitted(df_raw_Seqs_gen)


    # (Trimmed) Sequences splitted with processed indel.
    df_Seqs_splited_noIndel_gen = df_raw_Seqs_splitted_gen.apply(lambda x : ProcessIndel(x, _remove_indel=True), axis=0)
    print(df_Seqs_splited_noIndel_gen.head())


    # Making columns for gen file using nuc file.

    s_gen = df_Seqs_splited_noIndel_gen.iloc[0, :]
    s_nuc = df_Seqs_splited_noIndel_nuc.iloc[0, :]


    flag_Exon = s_gen.isin(s_nuc)
    # print(flag_Exon)


    l_idx_gen = []
    i_exon = 1
    i_intron = 1

    for i in range(0, len(flag_Exon)):

        if flag_Exon.iat[i]:
            l_idx_gen.append("Exon" + str(i_exon))
            i_exon += 1

        else:
            l_idx_gen.append("intron" + str(i_intron))
            i_intron += 1

    print(l_idx_gen)


    df_Seqs_splited_noIndel_gen.columns = pd.Index(l_idx_gen)
    df_Seqs_splited_noIndel_gen.index.name = "Alleles"

    print("\nRaw seqs of gen splitted.\n")
    print(df_Seqs_splited_noIndel_gen.head())
    df_Seqs_splited_noIndel_gen.to_csv(_out+".HLA_{0}.gen.raw.seqs.NoIndelChar.splitted.txt".format("A"), sep='\t', header=True, index=True)


    # (Trimmed) Sequences splitted with processed indel.
    df_Seqs_splited_Indel_gen = df_raw_Seqs_splitted_gen.apply(ProcessIndel, axis=0)
    df_Seqs_splited_Indel_gen.columns = df_Seqs_splited_noIndel_gen.columns
    print(df_Seqs_splited_Indel_gen.head())
    df_Seqs_splited_Indel_gen.to_csv(_out+".HLA_{0}.gen.seqs.IndelChar.splitted.txt".format("A"), sep='\t', header=True, index=True) # ***


    # Acquiring Genomic position of exon1 from "HLA_INTEGRATED_POSITIONS_hg18.txt"

    HLA_INTEGRATED_POS = pd.read_table('/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt',
                                       sep='\t', header=None, names=["HLA", "start", "end", "Type", "Direction"], index_col=0).loc["A", :]

    sr_temp = HLA_INTEGRATED_POS.loc[:, "Type"]
    exon1_offset = HLA_INTEGRATED_POS.loc[sr_temp == "exon1", "start"].iat[0]
    print(exon1_offset)

    l_genomic_positions = getPositionInfo("gen", df_Seqs_splited_Indel_gen, "A", isREVERSE["A"], exon1_offset)
    l_relative_positions = getPositionInfo("rel", df_Seqs_splited_Indel_gen, "A", isREVERSE["A"])

    df_Markers_gen = SeqsToMarkers(df_Seqs_splited_Indel_gen, l_genomic_positions, l_relative_positions) # ***
    print("\nFinal output as markers\n")
    print(df_Markers_gen.head())
    df_Markers_gen.to_csv(_out+".HLA_{0}.gen.MarkerTable.txt".format("A"), sep='\t', header=True, index=True)

    df_forMAP = df_Markers_gen.columns.to_frame(index=False)
    print("\nframes for .map file.\n")
    print(df_forMAP.head())
    df_forMAP.to_csv(_out+".HLA_{0}.gen.forMAP.txt".format("A"), sep='\t', header=True, index=False)


    df_Seqs_gen = df_Markers_gen.apply(lambda x : ''.join(x.astype(str)), axis=1) # *** Final dictionary
    df_Seqs_gen.to_csv(_out+".HLA_{0}.gen.seqs.txt".format("A"), sep='\t', header=False, index=True)

    print("\ndf_Markers_gen\n")
    print(df_Seqs_gen.head())


    # print("\n=============<PROT>=============\n")
    #
    # df_raw_Markers_prot = LoadRawSeq2(_prot, "A", "prot")
    # df_raw_Markers_prot.to_csv(_out+".HLA_{0}.prot.raw.markers.txt".format("A"), sep='\t', header=False, index=True)
    # print(df_raw_Markers_prot.head())
    # print(df_raw_Markers_prot.tail())


    return 0


# def LoadRawSeq(_nuc)

### Main Module 1
def LoadRawSeq2(_nuc_filename, _hla,  _type):


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
        offset_s = 3
        offser_e = -2
        allele_start = line_number_marks[0] + offset_s
        allele_end = line_number_marks[1] + offser_e
    elif _type == "gen":
        offset_s = 2
        offser_e = -3
        allele_start = line_number_marks[0] + offset_s
        allele_end = line_number_marks[1] + offser_e
    elif _type == "prot":
        offset_s = 2
        offser_e = -2
        allele_start = line_number_marks[0] + offset_s
        allele_end = line_number_marks[1] + offser_e


    NumberOfAlleles = allele_end - (allele_start - 1)

    print("\nallele start and end are {0} and {1}\n".format(allele_start, allele_end))
    print("# of Alleles is {0}\n".format(NumberOfAlleles))


    ##### READING_RAW_SEQUENCES

    f = open(_nuc_filename)
    f_lines = f.readlines()
    f.close()

    # print(f_lines[-2]) # "Please see http://hla.alleles.org/terms.html for terms of use."
    line_number_marks.append(len(f_lines) - 2)


    ### Acquiring starting offset

    # (2018. 9. 7.) Newly introduced.

    # Using "*_gen.txt" and "*_nuc.txt" files.

    if _type == "nuc":
        line_offset = re.split('\s+', f_lines[line_number_marks[0] + 1].lstrip(' ').rstrip(' \n')).pop()
    elif _type == "gen":
        line_offset = re.split('\s+', f_lines[line_number_marks[0]].lstrip(' ').rstrip(' \n')).pop()


    if _type == "nuc" or _type == "gen":
        print("\nStart Offset is : {0}\n".format(line_offset))


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

    # The rest of dictionaries are to be processed. (2nd chunks ~ The last chunks)

    for i in range(1, len(line_number_marks) - 1):
        t_keys = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(line_number_marks[i] + offset_s, line_number_marks[i + 1] + offser_e + 1)]
        t_s = [''.join(re.split('\s+', f_lines[(line_number_marks[i] + offset_s) + k].lstrip(' ').rstrip(' \n'))[1:]) for k in range(0, len(t_keys))]

        Sequence_Chunks.append({k: v for k, v in zip(t_keys, t_s)})




    ### Merging respective Dictionaries

    for i in range(0, len(k_first)):
        t_string = ''.join([Sequence_Chunks[j][k_first[i]] for j in range(0, len(Sequence_Chunks)) if k_first[i] in Sequence_Chunks[j].keys()])
        dict_first[k_first[i]] = tuple(t_string)

    # Unprocessed DataFrame of codon information.
    df_CODON = pd.DataFrame.from_dict(dict_first, orient='index')
    df_CODON.fillna(value='x', inplace=True)


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
            col[col == '-'] = conserved_1st_seq[i] # Processing '-' character.
            col[col == '*'] = 'x'


    return [df_CODON, line_offset]



### Main Module 2
def get_DF_rawSeqs(_df, _keep_index = True):


    if _keep_index:

        df_RETURN = _df.apply(lambda x : ''.join(x.astype(str)), axis = 1)
        df_RETURN.index = _df.index

        return df_RETURN

    else:

        return _df.apply(lambda x : ''.join(x.astype(str))).reset_index(drop=True)



### (Deprecate it later.)
def get_DF_rawSeqsSplitted(_sr, _keep_index=True):

    return _sr.str.split('\|', expand=True)



### Main Module 3
def getIndelSpots(_conserved_1st_seq):

    p = re.compile('\.+')

    iters = p.finditer(_conserved_1st_seq)
    l_INDELs = [i.span() for i in iters]

    if bool(l_INDELs):

        # First span
        l_notINDELs = [(0, l_INDELs[0][0])]

        for i in range(0, len(l_INDELs) - 1):
            l_notINDELs.append((l_INDELs[i][1], l_INDELs[i + 1][0]))

        # Last span
        l_notINDELs.append((l_INDELs[-1][1], len(_conserved_1st_seq)))

        return [l_INDELs, l_notINDELs]

    else:
        return l_INDELs



### Main Module 4
def getTrimmedSeqs(_string, _l_target_idx, _remove_indel=False):

    if len(_l_target_idx) == 2:

        IndelSeqs = pd.Series([_string[idx[0]:idx[1]] for idx in _l_target_idx[0]])

        if not _remove_indel:

            p_INSERTIONS = re.compile('\.+|\*+|x+')  # If Indel part corresponds to this pattern, then consider it as deletion('z')
            flag = IndelSeqs.str.match(p_INSERTIONS)

            IndelChar = pd.Series(['z' if flag[i] else 'Z' for i in range(0, IndelSeqs.shape[0])])

        else:

            IndelChar = pd.Series(['' for i in range(0, IndelSeqs.shape[0])])


        # =============

        t_string = _string[_l_target_idx[1][0][0]:_l_target_idx[1][0][1]]

        for i in range(1, len(_l_target_idx[1])):
            t_string = IndelChar[i - 1].join([t_string, _string[_l_target_idx[1][i][0]:_l_target_idx[1][i][1]]])

        return t_string

    else:

        return _string


def ProcessIndel(_sr, _remove_indel=False):

    # Acquiring span information using "conserved_1st_seq"
    l_spanInfo = getIndelSpots(_sr[0])

    return _sr.apply(lambda x: getTrimmedSeqs(x, l_spanInfo, _remove_indel))




### Main Module 5
def getPositionInfo(_type, _df, _hla, _isReverse, _exon1_offset=1, _as_flattened=False):

    if not (_type == "gen" or _type == "rel"):
        return -1

    if _type == "gen":
        offset_UTR = _exon1_offset
        offset_Rest = _exon1_offset
    elif _type == "rel":
        offset_UTR = 0  # Previously -1
        offset_Rest = 1

    # Preparing input dataframe
    df_UTR = _df.iloc[0, 0]  # Single String
    df_Rest = _df.iloc[0, 1:]  # Series

    # Final output to be returned.
    l_RETURN = []

    ### Processing UTR

    l_UTR = []
    N_indel = 0

    if not _isReverse:

        for i in range(0, len(df_UTR)):

            if (df_UTR[i] == 'A') or (df_UTR[i] == 'C') or (df_UTR[i] == 'G') or (df_UTR[i] == 'T'):
                #                 l_UTR.append((offset_UTR + N_indel) - i)
                l_UTR.append(offset_UTR - (len(df_UTR) - 1 - i) - N_indel)
            else:
                l_UTR.append('i')
                N_indel += 1

    else:

        print("A little bit later.")

    l_RETURN.append(l_UTR)


    ### Processing Rest of df

    curr_POS = offset_Rest

    for i in range(0, df_Rest.shape[0]):

        curr_string = df_Rest.iat[i]
        l_Rest = []
        N_indel = 0

        if not _isReverse:

            for j in range(0, len(curr_string)):

                if (curr_string[j] == 'A') or (curr_string[j] == 'C') or (curr_string[j] == 'G') or (
                        curr_string[j] == 'T'):
                    l_Rest.append((curr_POS - N_indel) + j)
                else:
                    l_Rest.append('i')
                    N_indel += 1

        else:

            print("A little bit later.")

        curr_POS = l_Rest[-1] + 1
        #         print(l_Rest)
        l_RETURN.append(l_Rest)

    if _as_flattened:
        return [element for l in l_RETURN for element in l]
    else:
        return l_RETURN


### Module 6
def SeqsToMarkers(_df, _l_gen_pos, _l_rel_pos):

    l_processed_sr = []

    for i in range(0, _df.shape[1]):

        sr_temp = _df.iloc[:, i].apply(lambda x: tuple(x))

        # Making it as markers
        df_temp = pd.DataFrame(sr_temp.tolist(), index=sr_temp.index)

        # Setting Columns(columns)
        idx_relative_pos = _l_rel_pos[i]
        idx_genomic_pos = _l_gen_pos[i]
        idx_TypeInfo = [sr_temp.name for j in range(0, len(sr_temp.iat[0]))]

        df_temp.columns = pd.MultiIndex.from_arrays([idx_TypeInfo, idx_relative_pos, idx_genomic_pos])

        l_processed_sr.append(df_temp)

    df_RETURN = pd.concat(l_processed_sr, axis=1)
    df_RETURN.index.name = None
    df_RETURN.columns.names = ["Type", "relative_POS", "genomic_POS"]

    return df_RETURN



### Sub-module 1
def ComplementChar(_character):
    if _character == "A":
        return "T"
    elif _character == "C":
        return "G"
    elif _character == "G":
        return "C"
    elif _character == "T":
        return "A"
    else:
        return _character



### Sub-module 2
def ComplementStr(_string):

    l_result = [ComplementChar(_string[i]) for i in range(0, len(_string))]

    return ''.join(l_result)







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

    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("-HLA", help="\nHLA gene name which you will process.\n\n", required=True, metavar='HLA',
                        choices=["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"])

    parser.add_argument("--type", "-t", help="\nSequence type to deal with(Amino Acids[AA] or SNPs[SNPS]\n\n", required=True, choices=["AA", "SNPS"], metavar='TYPE')
    parser.add_argument("-hg", help="\nHLA gene position information table.\n\n", required=True)
    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version", required=True)

    parser.add_argument("-nuc", help="\nInput *_nuc.txt file.\n\n", required=True)
    parser.add_argument("-gen", help="\nInput *_gen.txt file.\n\n", default="Not_given")
    parser.add_argument("-prot", help="\nInput *_prot.txt file.\n\n", default="Not_given")


    # for Publish
    # args = parser.parse_args()

    # for Testing
    # args = parser.parse_args(["-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt",
    #                           "-o", "NEW_Dict"])

    args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
                              "-HLA", "A",
                              "--type", "SNPS",
                              "-imgt", "370"
                              "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
                              "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
                              "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
                              "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt"
                              ])

    print(args)

    main(_out = args.o, _hla=args.HLA, _type=args.type, _hg_Table=args.hg, _imgt=args.imgt,
         _nuc=args.nuc, _gen=args.gen, _prot=args.prot)




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