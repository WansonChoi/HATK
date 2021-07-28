# -*- coding: utf-8 -*-

import os, sys, re, subprocess
from os.path import basename, dirname
import pandas as pd
import argparse, textwrap




########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (basename(__file__))




def ProcessIMGTv2(_out, _hla, _imgt, _BP_start_codon, isREVERSE, _gen, _nuc, _prot, _no_Ins=False, _include_UTR=False,
                  _save_intermediates=False):

    """

    # SNPS dictionary
    # AA dictionary

    """


    ### < (1) SNPS dictionary > ###

    # Raw sequence as each position
    df_gen_raw_markers, start_offset_gen = LoadRawSeq2(_gen, _hla, "gen")
    # print("df_gen_raw_markers:\n{}\n".format(df_gen_raw_markers))
    # print("start_offset_gen:{}".format(start_offset_gen))

    # Raw Seqeunces splitted by type (ex. 5_prime, exon1, ...)
    df_gen_raw_seqs_splitted = df_gen_raw_markers \
        .apply(lambda x: ''.join(x.astype(str)), axis=1) \
        .str.split('\|', expand=True)
    setTypeHeader_gen(df_gen_raw_seqs_splitted)
    #     print("df_gen_raw_seqs_splitted:\n{}\n".format(df_gen_raw_seqs_splitted))

    # Processing Insertion spots
    df_gen_seqs_splitted_InsAsZ = ProcessInsertion(df_gen_raw_seqs_splitted, "gen", _hla)
    #     print("df_gen_seqs_splitted_InsAsZ:\n{}\n".format(df_gen_seqs_splitted_InsAsZ))

    if _save_intermediates:
        df_gen_raw_seqs_splitted.to_csv(_out + '.HLA-{}.df_gen_raw_seqs_splitted.txt'.format(_hla), sep='\t',
                                        header=True, index=True)
        df_gen_seqs_splitted_InsAsZ.to_csv(_out + '.HLA-{}.df_gen_seqs_splitted_InsAsZ.txt'.format(_hla),
                                           sep='\t', header=True, index=True)



    df_gen_markers = get_PositionTable_SNPS(df_gen_seqs_splitted_InsAsZ, int(_BP_start_codon), start_offset_gen, isREVERSE)
    print("df_gen_markers:\n{}\n".format(df_gen_markers))
    """
    The multi-index of the `df_gen_markers` is needed in AA.
    """

    df_gen_PositionTable = df_gen_markers.iloc[[0], :].transpose()
    print("df_gen_PositionTable:\n{}\n".format(df_gen_PositionTable))

    if _save_intermediates:
        df_gen_markers.to_csv(_out + '.HLA-{}.df_gen_markers.txt'.format(_hla),
                              sep='\t', header=True, index=True)

        df_gen_PositionTable.to_csv(_out+'.HLA-{}.df_gen_PositionTable.txt'.format(_hla),
                                    sep='\t', header=True, index=True)





    ### < (2) AA dictionary > ###

    df_prot_raw_markers, _ = LoadRawSeq2(_prot, _hla, "prot")
    # print("df_prot_raw_markers:\n{}\n".format(df_prot_raw_markers))

    df_prot_seqs_splitted = df_prot_raw_markers \
        .apply(lambda x: ''.join(x.astype(str)), axis=1) \
        .str.split('\|', expand=True)  # Redundant but effective.
    df_prot_seqs_splitted.columns = ['AA_seq']
    # print("df_prot_seqs_splitted:\n{}\n".format(df_prot_seqs_splitted))

    df_prot_seqs_splitted_InsAsZ = ProcessInsertion(df_prot_seqs_splitted, "prot", _hla)
    # print("df_prot_seqs_splitted_InsAsZ:\n{}\n".format(df_prot_seqs_splitted_InsAsZ))

    if _save_intermediates:
        df_prot_raw_markers.to_csv(_out + '.HLA-{}.df_prot_raw_markers.txt'.format(_hla), sep='\t', header=False, index=True)
        df_prot_seqs_splitted.to_csv(_out + '.HLA-{}.df_prot_seqs_splitted.txt'.format(_hla), sep='\t', header=False, index=True)
        df_prot_seqs_splitted_InsAsZ.to_csv(_out + '.HLA-{}.df_prot_seqs_splitted_InsAsZ.txt'.format(_hla), sep='\t',
                                            header=False, index=True)



    # Loading *.nuc file for AA rel_pos start offset.
    """
    AA dictionary is dependent on the next two other things.
    
    (1) *.nuc.txt file
    : The start of relative position.
    
    (2) `df_SNPS_forMAP`
    : BP of AA map file.
    
    """
    _, AA_rel_pos_start = LoadRawSeq2(_nuc, _hla, "nuc")
    print("AA_rel_pos_start: {}".format(AA_rel_pos_start))

    df_SNPS_forMAP = df_gen_markers.columns.to_frame(index=False).reindex(["SNP_rel_pos", "SNP_gen_pos", "Type"], axis=1)
    print("df_SNPS_forMAP:\n{}\n".format(df_SNPS_forMAP))



    df_AA_PositionTable = get_PositionTable_AA(df_prot_seqs_splitted_InsAsZ, int(AA_rel_pos_start), df_SNPS_forMAP)
    print("df_AA_PositionTable:\n{}\n".format(df_AA_PositionTable))

    df_AA_forMAP = df_AA_PositionTable[['AA_rel_pos', 'SNP_gen_pos', 'Type']]
    df_AA_forMAP.columns = ['AA_rel_pos', 'AA_gen_pos', 'Type']
    print("df_AA_forMAP:\n{}\n".format(df_AA_forMAP))

    df_prot_markers = \
        pd.DataFrame([list(aa_seq) for aa_seq in df_prot_seqs_splitted_InsAsZ['AA_seq']],
                     index=df_prot_seqs_splitted_InsAsZ.index,
                     columns=pd.MultiIndex.from_frame(df_AA_forMAP))
    print("df_prot_markers:\n{}\n".format(df_prot_markers)) # `df_prot_markers` is the '__MAPTABLE__`.

    if _save_intermediates:
        df_AA_PositionTable.to_csv(_out+'.HLA-{}.df_prot_PositionTable.txt'.format(_hla), sep='\t', header=True, index=False)
        df_prot_markers.to_csv(_out+'.HLA-{}.df_prot_markers.txt'.format(_hla), sep='\t', header=True, index=True)

        df_SNPS_forMAP.to_csv(_out+'.HLA-{}.df_gen_SNPS_forMAP.txt'.format(_hla), sep='\t', header=True, index=False)
        df_AA_forMAP.to_csv(_out+'.HLA-{}.df_prot_AA_forMAP.txt'.format(_hla), sep='\t', header=True, index=False)





    ### < (4) Filter > ###

    # Filtering condition to `df_gen_markers`, `df_prot_markers`, `df_SNPS_MAP`, and `df_SNPS_MAP`
    if _no_Ins:

        ## Filtering dictionary
        # SNPS
        df_gen_markers, df_SNPS_forMAP = FilterInsertion(df_gen_markers)
        print("df_gen_markers(No Insertion):\n{}\n".format(df_gen_markers))
        print("df_SNPS_forMAP(No Insertion):\n{}\n".format(df_SNPS_forMAP))

        # AA
        df_prot_markers, df_AA_forMAP = FilterInsertion(df_prot_markers)
        print("df_prot_markers(No Insertion):\n{}\n".format(df_prot_markers))
        print("df_AA_forMAP(No Insertion):\n{}\n".format(df_AA_forMAP))


        if _save_intermediates:
            df_gen_markers.to_csv(_out + '.HLA-{}.df_gen_markers.noIns.txt'.format(_hla),
                                  sep='\t', header=True, index=True)
            df_SNPS_forMAP.to_csv(_out+'.HLA-{}.df_gen_SNPS_forMAP.noIns.txt'.format(_hla), sep='\t', header=True, index=False)

            df_prot_markers.to_csv(_out+'.HLA-{}.df_prot_markers.noIns.txt'.format(_hla), sep='\t', header=True, index=True)
            df_AA_forMAP.to_csv(_out+'.HLA-{}.df_prot_AA_forMAP.noIns.txt'.format(_hla), sep='\t', header=True, index=False)



    if not _include_UTR:

        ## Filtering dictionary (Only to SNPS)
        df_gen_markers, df_SNPS_forMAP = FilterUTR(df_gen_markers) # New 'df_SNPS_forMAP' is generated based on the filtered 'df_gen_markers'.
        print("df_gen_markers(No UTRs):\n{}\n".format(df_gen_markers))
        print("df_SNPS_forMAP(No UTRs):\n{}\n".format(df_SNPS_forMAP))





    ### < (5) Complement > ###
    if isREVERSE:
        df_gen_markers = df_gen_markers.applymap(lambda x : ComplementBase(x))
        print("df_gen_markers(Complemented):\n{}\n".format(df_gen_markers))





    ### < (6) as Dictionary > ###
    df_gen_dictionary = df_gen_markers.apply(lambda x: ''.join(x.astype(str)), axis=1)
    print("df_gen_dictionary:\n{}\n".format(df_gen_dictionary))

    df_prot_dictionary = df_prot_markers.apply(lambda x: ''.join(x.astype(str)), axis=1)
    print("df_prot_dictionary:\n{}\n".format(df_prot_dictionary))





    ### < (3) forMAP > ###

    # make map file for `df_SNPS_forMAP` and `df_AA_forMAP`.
    df_SNPS_MAP = MakeMap(_hla, "SNPS", df_SNPS_forMAP)
    print("df_SNPS_MAP:\n{}\n".format(df_SNPS_MAP))

    df_AA_MAP = MakeMap(_hla, "AA", df_AA_forMAP)
    print("df_AA_MAP:\n{}\n".format(df_AA_MAP))





    ### < (7) dimension check > ###
    print("<SNPS (HLA-{})>".format(_hla))
    print("reference sequence length : ", len(df_gen_dictionary.iat[0]))
    print("# of markers : ", df_SNPS_MAP.shape[0])
    if len(df_gen_dictionary.iat[0]) != df_SNPS_MAP.shape[0]: print(std_WARNING_MAIN_PROCESS_NAME + "Not same.")

    print("<AA (HLA-{})>".format(_hla))
    print("reference sequence length : ", len(df_prot_dictionary.iat[0]))
    print("# of markers : ", df_AA_MAP.shape[0])
    if len(df_prot_dictionary.iat[0]) != df_AA_MAP.shape[0]: print(std_WARNING_MAIN_PROCESS_NAME + "Not same.")





    # [df_Seqs_gen, df_Seqs_IndelProcessed_prot, final_SNPS_forMAP, final_AA_forMAP, __MAPTABLE__] ... Previously.
    return (df_gen_dictionary, df_prot_dictionary, df_SNPS_MAP, df_AA_MAP, df_prot_markers)





def LoadRawSeq2(_nuc_filename, _hla,  _type):


    # print("\nLoading \"{0}\" file".format(_nuc_filename))

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

    # print("\nMarks is {0}".format(line_number_marks))

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

    # print("\nallele start and end are {0} and {1}\n".format(allele_start, allele_end))
    # print("# of Alleles is {0}\n".format(NumberOfAlleles))


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
    else:
        line_offset = -1


    # if _type == "nuc" or _type == "gen":
    #     print("\nStart Offset is : {0}\n".format(line_offset))


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


    for i in range(0, len(conserved_1st_seq)):

        if conserved_1st_seq[i] != '.':

            col = df_CODON.iloc[1:, i]
            col[col == '-'] = conserved_1st_seq[i] # Processing '-' character.
            col[col == '*'] = 'x'


    return [df_CODON, line_offset]



def setTypeHeader_gen(_df):

    ### Set Header(_df.columns) inplace.

    """
    ex. N_columns = 17

    17 - 2(5-prime UTR, 3-prime UTR) = 15 (*** The reason why 2 is subtracted from `N_columns`)

    15 / 2 = 7 (=> 7 exon-intron pairs + 1 last exon)

    1 : exon1, intron1
    2 : exon2, intron2
    ...
    7 : exon7, intron7
    8 : exon8.

    """

    N_columns = _df.shape[1]
    N_exonN_intronN_pairs = int((N_columns - 2) / 2)

    l_exonN_intronN = [item for i in range(N_exonN_intronN_pairs) for item in
                       ('exon{N}'.format(N=i + 1), 'intron{N}'.format(N=i + 1))]
    #     print("l_exonN_intronN:{}".format(l_exonN_intronN))

    l_Header = ['5_prime'] + l_exonN_intronN + ['exon{}'.format(N_exonN_intronN_pairs + 1)] + ['3_prime']
    print("l_Header:{}".format(l_Header))

    _df.columns = l_Header

    ### set index name(_df.index.name)
    _df.index.name = 'Allele'

    return 0



def getInsertionSpots(_sequence, _seq_type, _hla):

    if re.match(r'^\.+$', _sequence):
        """
        # Exception. (ex. Amino acid seq. DQB1 exon5.)

        (ex.)
        `_sequence` := '....................'

        """
        #
        return []

    p_insertion = re.compile(r'\.+')

    """
    (ex.)
    '.ACG.TAC...GTACGT...' (len == 20)

    (Insertion) => [(0, 1), (4, 5), (8, 11), (17, 20)]

    (Return) => [(0', 0)] + [(1, 4), (5, 8), (11, 17)] + [(20, len(_sequence))]


    """

    insertions = [item.span() for item in p_insertion.finditer(_sequence)]
    #     print(insertions)

    if _seq_type == "prot" and _hla == 'DQB1':
        """
        (2018. 9. 16.)
        This exception handling was introduced becuase the 5th exon in DQB1 is either exon or indel itself.
        Also, this exceoption need to be gone through only when dealing with AA not SNPS(I hope so.).
        """

        #         print("Skip the last insertion considering it as the 5th exon.")
        insertions.pop()  # Abandon last indel assuming this will correspond to exon5.

    if len(insertions) > 0:
        pass
    else:
        ### No insertion found.
        #         print("No insertion.")
        return ([], [])

    ## First item.
    sequence_toSave = [(0, insertions[0][0])]

    ## Main items.
    for i in range(1, len(insertions)):
        #         print("{}-th : {}".format(i, insertions[i]))

        sequence_toSave.append((insertions[i - 1][1], insertions[i][0]))

    ## Last item.
    sequence_toSave.append((insertions[-1][1], len(_sequence)))

    #     print(sequence_toSave)
    return (insertions, sequence_toSave)



def getTrimmedSeqs(_sequence, _l_insertion_spots):
    (insertions, sequence_toSave) = _l_insertion_spots

    #     print(insertions)
    #     print(sequence_toSave)

    ### z or Z (Insertion check.)
    l_z_or_Z = []
    for ins in insertions:

        t_string = _sequence[ins[0]: ins[1]]
        #         print(t_string)

        if re.match(r'^(\.+|\*+|x+)$', t_string):
            # Normal (Same to the reference seq.)
            l_z_or_Z.append('z')
        else:
            # Insertion (Different to the reference seq.)
            l_z_or_Z.append('Z')

    #     print(l_z_or_Z)

    ### Concatenation
    to_concat = _sequence[sequence_toSave[0][0]: sequence_toSave[0][1]]
    for i in range(1, len(sequence_toSave)):
        t_string = _sequence[sequence_toSave[i][0]: sequence_toSave[i][1]]

        to_concat = l_z_or_Z[i - 1].join([to_concat, t_string])
    #     print(to_concat)

    return to_concat



def ProcessInsertion(_df, _seq_type, _hla, _out=None):
    """
    (1) Find insertion spots in the reference sequence. (`getInsertionSpots()`)
    (2) Trim insertion spots in virtual sequences (`getTrimmedSeqs()`)

    This iterates over each column of the `_df`.

    """

    l_sr_col_trimmed = []

    for _type, _sr_col in _df.iteritems():  # Iterates over each column of the `_df`.

        #         print("============================================")
        #         print("_sr_col:\n{}\n".format(_sr_col))
        #         print("Length of each string : ", len(_sr_col.iat[0]))

        ### Acquire insertion spots.
        t_reference_seq = _sr_col.iat[0]
        #         print("t_reference_seq:\n", t_reference_seq)
        #         print("_seq_type: ", _seq_type)
        #         print("_hla: ", _hla)

        l_insertion_spots = getInsertionSpots(t_reference_seq, _seq_type, _hla)
        #         print("l_insertion_spots:\n", l_insertion_spots)

        if len(l_insertion_spots[0]) > 0:

            sr_col_trimmed = _sr_col.map(lambda x: getTrimmedSeqs(x, l_insertion_spots))
        #             print("sr_col_trimmed:\n{}\n".format(sr_col_trimmed))

        else:
            # No insertion in the reference seq. (`t_reference_seq`)
            #             print("No insertion in '{}'.".format(_type))

            sr_col_trimmed = _sr_col

        l_sr_col_trimmed.append(sr_col_trimmed)

    df_RETURN = pd.concat(l_sr_col_trimmed, axis=1)
    #     print("df_RETURN:\n{}\n".format(df_RETURN))

    if _out:
        df_RETURN.to_csv(_out, sep='\t', header=True, index=True)
        return _out
    else:
        return df_RETURN



def get_PositionTable_SNPS(_df, _BP_start_codon, _rel_start=-1, isREVERSE=False):
    """

    Get
    (1) genomic position,
    (2) relative position,
    (3) Type

    => These 3 information columns will be pd.MultiIndex object.

    """

    l_rel_pos = []
    l_gen_pos = []

    sr_reference_seq = _df.iloc[0, :]
    #     print("sr_reference_seq:\n{}\n".format(sr_reference_seq))

    ### Divide (1) UTR(5-prime) and (2) the rest(exon1 ~)
    str_5_prime_UTR = sr_reference_seq.iloc[0]
    #     print("str_5_prime_UTR:\n{}".format(str_5_prime_UTR))

    sr_Rest = sr_reference_seq.iloc[1:]
    #     print("sr_Rest:\n{}".format(sr_Rest))



    ### (1) UTR(5-prime) processing.

    ## Relative position (Forward/Backward direction doesn't matter.)

    l_rel_pos_5prime_UTR = []
    l_gen_pos_5prime_UTR = []

    i_rel_pos = -1
    i_gen_pos = _BP_start_codon + (+1 if isREVERSE else -1)

    for base in reversed(str_5_prime_UTR):

        if base == 'z' or base == 'Z':
            #             print("{}: {}".format(i+0.5, base))

            # rel_pos
            l_rel_pos_5prime_UTR.append('i')

            # gen_pos
            l_gen_pos_5prime_UTR.append('i')

        else:
            #             print("{}: {}".format(i, base))

            # rel_pos
            l_rel_pos_5prime_UTR.append(i_rel_pos)
            i_rel_pos -= 1

            # gen_pos
            l_gen_pos_5prime_UTR.append(i_gen_pos)
            i_gen_pos = (i_gen_pos + 1) if isREVERSE else (
                        i_gen_pos - 1)  # (***) depending on forward/backward direction.

    #     if l_rel_pos_5prime_UTR[-1] != _rel_start:
    #         print("[Warning] Start relative position doesn't match!") # Add it later.

    #     print("l_rel_pos_5prime_UTR(reversed):\n{}".format(l_rel_pos_5prime_UTR[::-1]))
    #     print("l_gen_pos_5prime_UTR(reversed):\n{}".format(l_gen_pos_5prime_UTR[::-1]))

    l_rel_pos.append(l_rel_pos_5prime_UTR[::-1])
    l_gen_pos.append(l_gen_pos_5prime_UTR[::-1])



    ### (2) the rest(exon1 ~)

    i_rel_pos = 1  # (***) starts from 1 not 0.
    i_gen_pos = _BP_start_codon

    for (_type, _sequence) in sr_Rest.iteritems():

        t_l_rel_pos_Rest = []
        t_l_gen_pos_Rest = []

        # print("\n{} : {}".format(_type, _sequence))

        for base in _sequence:
            if base == 'z' or base == 'Z':
                # rel_pos
                t_l_rel_pos_Rest.append('i')

                # gen_pos
                t_l_gen_pos_Rest.append('i')

            else:
                # rel_pos
                t_l_rel_pos_Rest.append(i_rel_pos)
                i_rel_pos += 1

                # gen_pos
                t_l_gen_pos_Rest.append(i_gen_pos)
                i_gen_pos += (-1 if isREVERSE else +1)

        # print("t_l_rel_pos_Rest:\n{}".format(t_l_rel_pos_Rest))
        # print("t_l_gen_pos_Rest:\n{}".format(t_l_gen_pos_Rest))

        l_rel_pos.append(t_l_rel_pos_Rest)
        l_gen_pos.append(t_l_gen_pos_Rest)

    ### Return

    #     print("l_rel_pos:\n{}".format(l_rel_pos))
    #     print("l_gen_pos:\n{}".format(l_gen_pos))

    """
    The part that was `SeqsToMarkers()`
    """

    l_temp = []

    count = 0

    for i in range(_df.shape[1]):  # iter over columns(5_prime, exon1, ...)

        sr_temp = _df.iloc[:, i].apply(lambda x: tuple(x))
        #         print(sr_temp)

        df_temp = pd.DataFrame(
            sr_temp.tolist(),
            index=sr_temp.index,
            columns=pd.MultiIndex.from_arrays([l_rel_pos[i], l_gen_pos[i], [sr_temp.name] * len(sr_temp.iat[0])])
        )
        #         print(df_temp)

        l_temp.append(df_temp)

        count += 1
    #         if count >= 1 : break

    df_RETURN = pd.concat(l_temp, axis=1)
    df_RETURN.columns.names = ["SNP_rel_pos", "SNP_gen_pos", "Type"]

    #     print("df_RETURN:\n{}\n".format(df_RETURN))

    return df_RETURN



def get_PositionTable_AA(_df, _AA_rel_pos_start, _df_SNPS_forMAP):

    """

    'getPositionInfo_AA' in the former version.

    """

    AA_reference_seq = _df.iat[0, 0]
    print("AA_reference_seq:\n{}\n".format(AA_reference_seq))

    f_insertion = _df_SNPS_forMAP['SNP_rel_pos'] == 'i'  # (***) exclude insertion spots
    print(f_insertion.any())
    f_exon = _df_SNPS_forMAP['Type'].str.match('exon')  # Only exons.
    df_SNPS_forMAP_exons = _df_SNPS_forMAP[(~f_insertion) & f_exon]
    print("df_SNPS_forMAP_exons:\n{}\n".format(df_SNPS_forMAP_exons))




    ### Relative position.

    l_AA_rel_pos = []
    t_rel_pos = _AA_rel_pos_start

    for base in AA_reference_seq:
        if base == 'z' or base == 'Z':
            l_AA_rel_pos.append('i')
        else:
            l_AA_rel_pos.append(t_rel_pos)
            t_rel_pos += 1

            if t_rel_pos == 0:
                t_rel_pos += 1

    print("l_AA_rel_pos:\n{}".format(l_AA_rel_pos))
    print("len(l_AA_rel_pos): ", len(l_AA_rel_pos))

    df_AA_rel_pos = pd.DataFrame({'AA_rel_pos': l_AA_rel_pos, 'AA_seq': list(AA_reference_seq)})
    print("df_AA_rel_pos:\n{}\n".format(df_AA_rel_pos))




    ### Genomic position

    # Mid point of BPs of exon parts.
    # df_SNPS_forMAP_exons_codon_mid = pd.concat(
    #     [df_SNPS_forMAP_exons.iloc[[i], 1:] for i in range(1, df_SNPS_forMAP_exons.shape[0], 3)])
    df_SNPS_forMAP_exons_codon_mid = \
        df_SNPS_forMAP_exons.iloc[list(range(1, df_SNPS_forMAP_exons.shape[0], 3)), :]
    print("df_SNPS_forMAP_exons_codon_mid:\n{}\n".format(df_SNPS_forMAP_exons_codon_mid))



    # considers insertion spots again.
    l_temp = []
    index_codon_mid = 0  # index for `df_SNPS_forMAP_exons_codon_mid`.

    for i in range(df_AA_rel_pos.shape[0]):

        base = df_AA_rel_pos.iat[i, 1]

        if base == 'z' or base == 'Z':
            l_temp.append(['i', 'i', 'i'])
        else:
            SNP_rel_pos = df_SNPS_forMAP_exons_codon_mid.iat[index_codon_mid, 0]
            SNP_gen_pos = df_SNPS_forMAP_exons_codon_mid.iat[index_codon_mid, 1]
            Type = df_SNPS_forMAP_exons_codon_mid.iat[index_codon_mid, 2]

            l_temp.append([SNP_rel_pos, SNP_gen_pos, Type])
            index_codon_mid += 1

    df_SNPS_forMAP_exons_codon_mid_INS = pd.DataFrame(l_temp,
                                                      columns=['SNP_rel_pos', 'SNP_gen_pos', 'Type'],
                                                      index=df_AA_rel_pos.index)
    #         print("df_SNPS_forMAP_exons_codon_mid_INS:\n{}\n".format(df_SNPS_forMAP_exons_codon_mid_INS))

    df_AA_PositionTable = pd.concat([df_AA_rel_pos, df_SNPS_forMAP_exons_codon_mid_INS], axis=1)
    # print("df_AA_PositionTable:\n{}\n".format(df_AA_PositionTable))

    return df_AA_PositionTable



def MakeMap(_hla, _type, _df_forMAP):

    l_Label = []
    l_BP = []

    for i in range(_df_forMAP.shape[0]):

        Rel_pos = _df_forMAP.iat[i, 0]
        Gen_pos = _df_forMAP.iat[i, 1]
        Func_Type = _df_forMAP.iat[i, 2]

        if Rel_pos == 'i':

            # Rel_pos = 'x'.join([_df_forMAP.iat[i - 1, 0], _df_forMAP.iat[i + 1, 0]])
            Rel_pos = '{}x{}'.format(_df_forMAP.iat[i - 1, 0], _df_forMAP.iat[i + 1, 0])
            Gen_pos = int(round((int(_df_forMAP.iat[i - 1, 1]) + int(_df_forMAP.iat[i - 1, 1])) / 2))

            t_Label = '_'.join(['INS', _type, _hla, str(Rel_pos), str(Gen_pos), Func_Type])

        else:

            t_Label = '_'.join([_type, _hla, str(Rel_pos), str(Gen_pos), Func_Type])

        l_Label.append(t_Label)
        l_BP.append(Gen_pos)



    # Return
    df_MAP = pd.concat([
        pd.Series(['6'] * _df_forMAP.shape[0], name='Chr'),
        pd.Series(l_Label, name='Label'),
        pd.Series(['0'] * _df_forMAP.shape[0], name='GD'),
        pd.Series(l_BP, name='BP')
    ], axis=1)
    #     print("df_MAP:\n{}\n".format(df_MAP))

    return df_MAP



def FilterInsertion(_df):

    df_columns = _df.columns.to_frame()
    #     print(df_columns)

    ## Filtering the 'df_gen_markers' (or 'df_prot_markers').
    f_INS = df_columns.iloc[:, 0] == 'i'

    df_markers_NoINS = _df.loc[:, ~f_INS]
    #     print("df_markers_NoINS:\n{}\n".format(df_markers_NoINS))


    ## Filtering the 'df_gen_forMAP' (or 'df_prot_forMAP'). (Generating a new one based on the above filtered 'df_gen_markers' (or 'df_prot_markers')
    df_forMAP_NoINS = df_markers_NoINS.columns.to_frame(index=False)


    return df_markers_NoINS, df_forMAP_NoINS



def FilterUTR(_df_gen_markers):

    df_columns = _df_gen_markers.columns.to_frame()
    # print("df_columns:\n{}\n".format(df_columns))

    """
    (ex.)
    
    df_columns:
                                    SNP_rel_pos SNP_gen_pos     Type
    SNP_rel_pos SNP_gen_pos Type                                    
    -300        30457009    5_prime        -300    30457009  5_prime
    -299        30457010    5_prime        -299    30457010  5_prime
    -298        30457011    5_prime        -298    30457011  5_prime
    ...                                     ...         ...      ...
     3520       30460828    3_prime        3520    30460828  3_prime
     3521       30460829    3_prime        3521    30460829  3_prime
     3522       30460830    3_prime        3522    30460830  3_prime
    
    
    (cf. Insertion)
    -13	    30457296	5_prime
    i	    i	        5_prime # Insertion is marked with type.
    -12	    30457297	5_prime
    
    
    """

    ## Filtering the 'df_gen_markers' (No '5_prime' and '3_prime')

    f_UTR = df_columns['Type'].str.match(r'[53]_prime')

    df_gen_markers_NoUTR = _df_gen_markers.loc[:, ~f_UTR]
    # print("df_gen_markers_NoUTR:\n{}\n".format(df_gen_markers_NoUTR))


    ## Filtering map file (Generating a new 'df_gen_forMAP' based on the above filtered 'df_gen_markers'.
    df_SNPS_forMAP_NoUTR = df_gen_markers_NoUTR.columns.to_frame(index=False)
    # print("df_SNPS_forMAP_NoUTR:\n{}\n".format(df_SNPS_forMAP_NoUTR))

    """
    New 'forMAP' file is generated from the filtered 'df_gen_markers_NoUTR' DataFrame.
    i.e. Marker numbers will be consistent between 'df_gen_markers_NoUTR' and 'forMAP'.
    """

    return df_gen_markers_NoUTR, df_SNPS_forMAP_NoUTR



def ComplementBase(_base):

    if _base == "A":
        return "T"
    elif _base == "C":
        return "G"
    elif _base == "G":
        return "C"
    elif _base == "T":
        return "A"
    else:
        return _base





if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        ProcessIMGTv2.py

        : Processing(Parsing HLA sequence information distributed by IMGT-HLA.

        Renewed in 2018.09.10
        Renewed in 2021.07.22


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--out", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("--HLA", help="\nHLA gene name which you will process.\n\n", required=True, metavar='HLA')
    parser.add_argument("--imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="imgt_version",
                        required=True)
    parser.add_argument("--BP-start-codon", help="\nThe Base position of the start codon of the HLA gene.\n\n",
                        required=True)
    parser.add_argument("--isReverse", help="\nIs Backward(Reverse) strand? (<-> Forward strand).\n\n",
                        action="store_true")

    parser.add_argument("--nuc", help="\nInput *_nuc.txt file.\n\n", required=True)
    parser.add_argument("--gen", help="\nInput *_gen.txt file.\n\n", required=True)
    parser.add_argument("--prot", help="\nInput *_prot.txt file.\n\n", required=True)

    # Optional arguments
    parser.add_argument("--no-Ins", help="\nNo Insertion Markers.\n\n", action="store_true")
    parser.add_argument("--include-UTR", help="\ninclude UTR parts(5-prime, 3-prime) in Markers.\n\n", action="store_true")
    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files. (DEBUG)\n\n", action='store_true')





    ##### < for Test > #####

    # HLA-E / imgt3440
    args = parser.parse_args(["--out", "/Users/wansonchoi/Git_Projects/HATK/tests/20210728_IMGT2Seq/Dict.HLA-E",
                              "--HLA", "E",
                              "--imgt", "3320",
                              "--BP-start-codon", "30457309",
                              "--nuc", "/Users/wansonchoi/Git_Projects/HATK/example/IMGTHLA3320/alignments/E_nuc.txt",
                              "--gen", "/Users/wansonchoi/Git_Projects/HATK/example/IMGTHLA3320/alignments/E_gen.txt",
                              "--prot", "/Users/wansonchoi/Git_Projects/HATK/example/IMGTHLA3320/alignments/E_prot.txt",
                              "--save-intermediates",
                              "--no-Ins"
                              ])

    # HLA-E / imgt3440
    # args = parser.parse_args(["--out", "/Users/wansonchoi/Git_Projects/HATK/tests/20210706_IMGT2Seq/Dict.HLA-E",
    #                           "--HLA", "E",
    #                           "--imgt", "3440",
    #                           "--BP-start-codon", "30457309",
    #                           "--nuc", "/Users/wansonchoi/Git_Projects/HATK/example/IMGTHLA3320/alignments/E_nuc.txt",
    #                           "--gen", "/Users/wansonchoi/Git_Projects/HATK/example/IMGTHLA3320/alignments/E_gen.txt",
    #                           "--prot", "/Users/wansonchoi/Git_Projects/HATK/example/IMGTHLA3320/alignments/E_prot.txt",
    #                           "--save-intermediates"
    #                           ])



    ##### < for Publish > #####
    # args = parser.parse_args()

    print(args)

    # main function execution
    ProcessIMGTv2(_out=args.out, _hla=args.HLA, _imgt=args.imgt, _BP_start_codon=args.BP_start_codon, isREVERSE=args.isReverse,
                  _gen=args.gen, _nuc=args.nuc, _prot=args.prot,
                  _no_Ins=args.no_Ins, _include_UTR=args.include_UTR, _save_intermediates=args.save_intermediates)
