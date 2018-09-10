# -*- coding: utf-8 -*-

import os, sys, subprocess
import pandas as pd
import argparse, textwrap
import re

def IMGTtoSequences(_out, _hla, _type, _hg_Table, _imgt,
                    _nuc, _gen = "Not_given", _prot = "Not_given",
                    _no_Indel_AA=False, _no_Indel_SNPS=False):


    ########## < Core Variables > ##########


    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting IMGTtoSequences.")


    ### General
    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}



    ########## < Argument Checking. > ##########

    # Preparing intermediate paths.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')

    INTERMEDIATE_PATH = os.path.dirname(_out)

    if not os.path.exists(INTERMEDIATE_PATH):
        os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))





    ########## < Main > ##########

    ##### (2018/1/17) Preparing HLA position information(exon, intron, etc.)

    # Preparing start offset of exon1
    HLA_INTEGRATED_POS = pd.read_table(_hg_Table, sep='\t', header=None, usecols = [0,1,2,3],
                                       names=['HLA', 'start', 'end', 'Type', "Direction"], index_col=0).loc["A", :]
    print("\nLoaded HLA information table.\n")
    print(HLA_INTEGRATED_POS.head())

    sr_temp = HLA_INTEGRATED_POS.loc[:, "Type"]
    exon1_offset = HLA_INTEGRATED_POS.loc[sr_temp == "exon1", "start"].iat[0]
    print(exon1_offset)



    ##### Loading *_nuc.txt file.

    print("\n\n=============<NUC>=============\n")
    print(std_MAIN_PROCESS_NAME + "Processing *_nuc.txt file.\n")

    # Raw Markers
    df_raw_Markers_nuc, start_offset_nuc = LoadRawSeq2(_nuc, _hla, "nuc")
    # df_raw_Markers_nuc.to_csv(_out+".HLA_{0}.nuc.raw.markers.txt".format(_hla), sep='\t', header=False, index=True)
    # print(df_raw_Markers_nuc.head())

    # Raw Seqeunces
    df_raw_Seqs_nuc = get_DF_rawSeqs(df_raw_Markers_nuc)
    # df_raw_Seqs_nuc.to_csv(_out+".HLA_{0}.nuc.raw.seqs.txt".format(_hla), sep='\t', header=False, index=True)
    # print("\nRaw seqs of nuc\n")
    # print(df_raw_Seqs_nuc.head())

    # Raw Seqeunces splitted (for Exon checking)
    df_raw_Seqs_splitted_nuc = df_raw_Seqs_nuc.str.split('\|', expand=True)

    df_raw_Seqs_splitted_nuc.columns = pd.Index(["Exon"+str(i+1) for i in range(0, df_raw_Seqs_splitted_nuc.shape[1])])
    df_raw_Seqs_splitted_nuc.index.name = "Alleles"
    df_raw_Seqs_splitted_nuc.to_csv(_out+".HLA_{0}.nuc.raw.seqs.splitted.txt".format(_hla), sep='\t', header=True, index=True) # Exporting 1 for nuc

    print("\nRaw seqs(with Indels) of nuc splitted.\n")
    print(df_raw_Seqs_splitted_nuc.head())

    # (Trimmed) Sequences splitted without any indel character(for Exon checking with gen file)
    df_Seqs_splited_noIndel_nuc = df_raw_Seqs_splitted_nuc.apply(lambda x : ProcessIndel(x, _remove_indel=True), axis = 0)

    print("\n\nSeqeunces(without Indels) of nuc splitted.\n")
    print(df_Seqs_splited_noIndel_nuc.head())

    df_Seqs_splited_noIndel_nuc.to_csv(_out+".HLA_{0}.nuc.seqs.NoIndelChar.splitted.txt".format(_hla), sep='\t', header=True, index=True) # Exporting 2 for nuc

    # # No more processing to *_nuc.txt file is needed.




    if _type == "AA":

        print(std_MAIN_PROCESS_NAME + "Processing Prot files for AA.")


    elif _type == "SNPS":

        print("\n")
        print("\n=============<GEN>=============\n")
        print(std_MAIN_PROCESS_NAME + "Processing Prot files for SNPS.\n")

        # Raw Markers
        df_raw_Markers_gen, start_offset_gen = LoadRawSeq2(_gen, _hla, "gen")
        # df_raw_Markers_gen.to_csv(_out+".HLA_{0}.gen.raw.markers.txt".format("A"), sep='\t', header=False, index=True)
        print(df_raw_Markers_gen.head())
        # print(df_raw_Markers_gen.tail())

        # Raw Sequences
        df_raw_Seqs_gen = get_DF_rawSeqs(df_raw_Markers_gen)
        # df_raw_Seqs_gen.to_csv(_out+".HLA_{0}.gen.raw.seqs.txt".format(_hla), sep='\t', header=False, index=True)
        # print("\n\nRaw seqs(with Indels) of gen splitted.\n")
        # print(df_raw_Seqs_gen.head())
        # print(df_raw_Seqs_gen.tail())

        # Raw Seqeunces splitted (for Exon checking)
        df_raw_Seqs_splitted_gen = df_raw_Seqs_gen.str.split('\|', expand=True)

        # # (Trimmed) Sequences splitted without any indel character(for Exon checking with nuc file)
        # df_Seqs_splited_noIndel_gen = df_raw_Seqs_splitted_gen.apply(lambda x : ProcessIndel(x, _remove_indel=True), axis=0)
        # # print(df_Seqs_splited_noIndel_gen.head())


        ##### Exon checking

        # Making header for gen file using nuc file.

        df_Seqs_splited_noIndel_gen = df_raw_Seqs_splitted_gen.apply(lambda x : ProcessIndel(x, _remove_indel=True), axis=0)

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

        print("\n{0}\n".format(l_idx_gen))


        if _no_Indel_SNPS:

            df_Seqs_splited_noIndel_gen.columns = pd.Index(l_idx_gen)
            df_Seqs_splited_noIndel_gen.index.name = "Alleles"

            print("\nSeqs without Indel characters of gen splitted.\n")
            print(df_Seqs_splited_noIndel_gen.head())

            df_Seqs_splited_noIndel_gen.to_csv(_out+".HLA_{0}.gen.seqs.NoIndelChar.splitted.txt".format(_hla), sep='\t', header=True, index=True)

        else:

            # (Trimmed) Sequences splitted with processed indel.
            df_Seqs_splited_Indel_gen = df_raw_Seqs_splitted_gen.apply(ProcessIndel, axis=0)

            df_Seqs_splited_Indel_gen.columns = pd.Index(l_idx_gen)
            df_Seqs_splited_Indel_gen.index.name = "Alleles"

            print("\nSeqs with Indel characters of gen splitted.\n")
            print(df_Seqs_splited_Indel_gen.head())

            df_Seqs_splited_Indel_gen.to_csv(_out+".HLA_{0}.gen.seqs.IndelChar.splitted.txt".format(_hla), sep='\t', header=True, index=True) # ***



        ##### Making final output as markers and seqs.

        df_temp = df_Seqs_splited_Indel_gen if not _no_Indel_SNPS else df_Seqs_splited_noIndel_gen

        l_genomic_positions = getPositionInfo("gen", df_temp, _hla, isREVERSE[_hla], exon1_offset)
        l_relative_positions = getPositionInfo("rel", df_temp, _hla, isREVERSE[_hla])

        # Final output as markers
        df_Markers_gen = SeqsToMarkers(df_temp, l_genomic_positions, l_relative_positions) # ***

        print("\nFinal output as markers\n")
        print(df_Markers_gen.head())
        df_Markers_gen.to_csv(_out+".HLA_{0}.gen.MarkerTable.txt".format(_hla), sep='\t', header=True, index=True)


        df_forMAP = df_Markers_gen.columns.to_frame(index=False)
        print("\nframes for .map file.\n")
        print(df_forMAP.head())
        df_forMAP.to_csv(_out+".HLA_{0}.gen.forMAP.txt".format("A"), sep='\t', header=True, index=False)


        # Final output as Seqs
        df_Seqs_gen = df_Markers_gen.apply(lambda x : ''.join(x.astype(str)), axis=1) # *** Final dictionary
        df_Seqs_gen.to_csv(_out+".HLA_{0}.gen.seqs.txt".format("A"), sep='\t', header=False, index=True)

        print("\ndf_Markers_gen\n")
        print(df_Seqs_gen.head())


    else:

        print(std_ERROR_MAIN_PROCESS_NAME + "You gave wrong value to `_type` argument.")







#################### < Core Functions > ####################

### Main Module 1
def LoadRawSeq2(_nuc_filename, _hla,  _type):


    print("\nLoading \"{0}\" file".format(_nuc_filename))

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

        IMGTtoSequences.v2.py

        : Processing(Parsing HLA sequence information distributed by IMGT-HLA.
        
        Renewed in 2018.09.10


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("-HLA", help="\nHLA gene name which you will process.\n\n", required=True, metavar='HLA',
                        choices=["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"])

    parser.add_argument("--type", "-t", help="\nSequence type to deal with(Amino Acids[AA] or SNPs[SNPS]\n\n", required=True, choices=["AA", "SNPS"], metavar='TYPE')
    parser.add_argument("--hg-table", "-hg", help="\nHLA gene position information table.\n\n", required=True)
    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="imgt_version", required=True)

    parser.add_argument("-nuc", help="\nInput *_nuc.txt file.\n\n", required=True)
    parser.add_argument("-gen", help="\nInput *_gen.txt file.\n\n", default="Not_given")
    parser.add_argument("-prot", help="\nInput *_prot.txt file.\n\n", default="Not_given")

    # Optional arguments
    parser.add_argument("-nI-SNPS", help="\nNo Indel Markers in SNPS.\n\n", action="store_true")
    parser.add_argument("-nI-AA", help="\nNo Indel Markers in AA\n\n", action="store_true")




    ##### < for Test > #####

    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "A",
    #                           "--type", "SNPS",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt"
    #                           ])

    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "A",
    #                           "--type", "SNPS",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt",
    #                           "-nI-SNPS"
    #                           ])



    ##### < for Publish > #####
    args = parser.parse_args()

    print(args)


    # main function execution
    IMGTtoSequences(_out = args.o, _hla=args.HLA, _type=args.type, _hg_Table=args.hg_table, _imgt=args.imgt,
                    _nuc=args.nuc, _gen=args.gen, _prot=args.prot, _no_Indel_SNPS=args.nI_SNPS, _no_Indel_AA=args.nI_AA)
