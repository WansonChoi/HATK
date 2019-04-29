# -*- coding: utf-8 -*-

import os, sys, re, subprocess
import pandas as pd
import argparse, textwrap


########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}





def ProcessIMGT(_out, _hla, _hg, _imgt, _nuc=None, _gen=None, _prot=None, _no_Indel=False, _save_intermediates=False,
                _no_prime=True):


    ########## < Core Variables > ##########

    ### Paths
    p_data = "./data/IMGT2Sequence"


    ########## < Argument Checking. > ##########

    # Preparing intermediate paths.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)):
        INTERMEDIATE_PATH = os.path.dirname(_out)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)



    # if _nuc == "Not_given" or _prot == "Not_given" or _gen == "Not_given":
    if not (_nuc and _prot and _gen):
        print(std_ERROR_MAIN_PROCESS_NAME + "{0}_prot.txt, {0}_gen.txt or {0}_nuc.txt files can't be found. Please check it again.\n".format(_hla))
        sys.exit()



    ### "HLA_INTEGRATED_POSITIONS" file

    HLA_INTEGRATED_POSITIONS_filename = os.path.join(p_data, "HLA_INTEGRATED_POSITIONS_hg{0}.txt".format(_hg))

    if not os.path.exists(HLA_INTEGRATED_POSITIONS_filename):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}\" not found!".format(HLA_INTEGRATED_POSITIONS_filename))
        sys.exit()

    # print(HLA_INTEGRATED_POSITIONS_filename)





    ########## < Loading necessary Data. > ##########

    # (2018/1/17) Preparing HLA position information(exon, intron, etc.)

    HLA_INTEGRATED_POS = pd.read_table(HLA_INTEGRATED_POSITIONS_filename, sep='\t', header=None, usecols = [0, 1, 2, 3],
                                       names=['HLA', 'start', 'end', 'Type', "Direction"], index_col=0).loc[_hla, :]
    # print("\nLoaded HLA information table.\n")
    # print(HLA_INTEGRATED_POS.head())

    sr_temp = HLA_INTEGRATED_POS.loc[:, "Type"]
    exon1_offset = HLA_INTEGRATED_POS.loc[sr_temp == "exon1", "start"].iat[0] # `exon1_offset` will be used in generating SNPS dictionary.
    # print("\nStarting point of genomic position of HLA-{0} is {1}\n".format(_hla, exon1_offset))



    ########## < Control Flags. > ##########

    _1_LOADING_GEN = 1
    _2_GENERATING_SNPS_DICT = 1
    _3_LOADING_NUC = 1
    _4_LOADING_PROT = 1
    _5_GENERATING_AA_DICT = 1
    _6_GENERATING_forMAP = 1




    ########## < Main > ##########



    if _1_LOADING_GEN:


        # print("\n\n=============<GEN>=============\n")
        # print(std_MAIN_PROCESS_NAME + "Processing gen files for SNPS.\n")

        # Raw Markers
        df_raw_Markers_gen, start_offset_gen = LoadRawSeq2(_gen, _hla, "gen")
        # df_raw_Markers_gen.to_csv(_out+".HLA_{0}.gen.raw.markers.txt".format("A"), sep='\t', header=False, index=True)
        # print("1")
        # print(df_raw_Markers_gen.head())

        # Raw Sequences
        df_raw_Seqs_gen = get_DF_rawSeqs(df_raw_Markers_gen)
        # df_raw_Seqs_gen.to_csv(_out+".HLA_{0}.gen.raw.seqs.txt".format(_hla), sep='\t', header=False, index=True)
        # print("2")
        # print(df_raw_Seqs_gen.head())

        # Raw Seqeunces splitted (for Exon checking)
        df_raw_Seqs_splitted_gen = df_raw_Seqs_gen.str.split('\|', expand=True)
        # print("3")
        # print(df_raw_Seqs_splitted_gen.head())




    if _2_GENERATING_SNPS_DICT:


        ##### Exon/Intron Labeling

        # Making header for gen file using nuc file.
        # (2018. 9. 11.) Previous way to compare sequences directly was too rigid to be applied to all HLAs.

        df_Seqs_splited_noIndel_gen = df_raw_Seqs_splitted_gen.apply(lambda x : ProcessIndel(x, _remove_indel=True), axis=0)

        s_gen = df_Seqs_splited_noIndel_gen.iloc[0, :]

        l_idx_gen = ["5_prime"]
        i_exon = 1
        i_intron = 1


        for i in range(1, s_gen.shape[0]):

            if i % 2 == 0:
                l_idx_gen.append("intron" + str(i_intron))
                i_intron += 1
            else:
                l_idx_gen.append("exon" + str(i_exon))
                i_exon += 1

        l_idx_gen.pop()
        l_idx_gen.append("3_prime")

        # print("\n{0}\n".format(l_idx_gen))


        # (2018. 9. 13.) Processing NoIndel dataframe became necessary due to making map file of AA.
        df_Seqs_splited_noIndel_gen.columns = pd.Index(l_idx_gen)
        df_Seqs_splited_noIndel_gen.index.name = "Alleles"

        # print("\nSeqs without Indel characters of gen splitted.\n")
        # print(df_Seqs_splited_noIndel_gen.head())


        # File Writing
        if _save_intermediates:
            df_Seqs_splited_noIndel_gen.to_csv(_out+".HLA_{0}.gen.seqs.noindel.splitted.txt".format(_hla), sep='\t', header=True, index=True) # **




        ### Consenting with AA map file. ###

        l_genomic_positions = getPositionInfo_SNPS("gen", df_Seqs_splited_noIndel_gen, _hla, isREVERSE[_hla], exon1_offset,
                                                   _has_Indel=(not _no_Indel))
        l_relative_positions = getPositionInfo_SNPS("rel", df_Seqs_splited_noIndel_gen, _hla, isREVERSE[_hla],
                                                    _has_Indel=(not _no_Indel))


        df_Markers_NoIndel_gen = SeqsToMarkers(df_Seqs_splited_noIndel_gen, l_genomic_positions, l_relative_positions)


        if _no_Indel:

            print("\n\nMarkers without Indel characters.\n")
            print(df_Markers_NoIndel_gen.head())

            # File Writing
            if _save_intermediates:
                df_Markers_NoIndel_gen.to_csv(_out+".HLA_{0}.gen.markers.noindel.txt".format(_hla), sep='\t', header=True, index=True)


        # Filtering only exon columns in `df_Markers_NoIndel_gen`

        p = re.compile("(e|E)xon\d?")

        df_Markers_onlyExons =  df_Markers_NoIndel_gen.filter(regex=p, axis=1)
        # print("Only exons")
        # print(df_Markers_onlyExons.head())

        # df_Markers_onlyExons.to_csv(_out+".HLA_{0}.gen.onlyexons.markers.txt".format(_hla), sep='\t', header=True, index=True)

        precursor_AA_forMAP = df_Markers_onlyExons.columns.to_frame(index=False)
        precursor_AA_forMAP.columns = pd.Index(["SNP_rel_pos", "SNP_gen_pos", "Type"])

        # File Writing
        if _save_intermediates:
            precursor_AA_forMAP.to_csv(_out+".HLA_{0}.gen.precursor_AA_forMAP.noindel.txt".format(_hla), sep='\t', header=True, index=False)




        if not _no_Indel:

            # (Trimmed) Sequences splitted with processed indel.
            df_Seqs_splited_Indel_gen = df_raw_Seqs_splitted_gen.apply(ProcessIndel, axis=0)

            df_Seqs_splited_Indel_gen.columns = pd.Index(l_idx_gen)
            df_Seqs_splited_Indel_gen.index.name = "Alleles"

            print("\nSeqs with Indel characters of gen splitted.\n")
            print(df_Seqs_splited_Indel_gen.head())

            # File Writing
            if _save_intermediates:
                df_Seqs_splited_Indel_gen.to_csv(_out+".HLA_{0}.gen.seqs.indel.splitted.txt".format(_hla), sep='\t', header=True, index=True) # **

                """
                `df_Seqs_splited_Indel_gen`
                
                Alleles	        / 5_prime  / exon1    / intron1  / ... / 3_prime
                A*01:01:01:01   / CAGGA... / ATGGC... / GTGAG... / ... / GACAGCTGCCTT...
                A*01:01:01:02N  / CAGGA... / ATGGC... / GTGAG... / ... / GACAGCTGCCTT...
                ...
                
                """



        ##### Making final output as markers and seqs.

        df_temp = df_Seqs_splited_Indel_gen if not _no_Indel else df_Seqs_splited_noIndel_gen



        # ### (2019. 01. 21.) Excluding "5-prime" and "3-prime" sequence parts.
        # p = re.compile(r'exon\d|intron\d')
        # df_temp = df_temp.filter(regex=p, axis=1)
        #
        # print("df_temp")
        # print(df_temp.head())



        l_genomic_positions = getPositionInfo_SNPS("gen", df_temp, _hla, isREVERSE[_hla], exon1_offset,
                                                   _has_Indel=(not _no_Indel))
        l_relative_positions = getPositionInfo_SNPS("rel", df_temp, _hla, isREVERSE[_hla],
                                                    _has_Indel=(not _no_Indel))


        # for i in range(0, df_temp.shape[1]):
        #     print(len(l_genomic_positions[i]))
        #     print(len(l_relative_positions[i]))
        #     print(len(df_temp.iloc[0, i]))

        # Final output as markers
        df_Markers_gen = SeqsToMarkers(df_temp, l_genomic_positions, l_relative_positions)

        print("\nFinal output as markers\n")
        print(df_Markers_gen.head())



        if _no_prime:
            ### (2019. 01. 22.) Excluding prime sequences.
            p = re.compile(r'exon\d|intron\d')
            df_Markers_gen = df_Markers_gen.filter(regex=p, axis=1)

            # print("Prime filetered")
            # print(df_Markers_gen.head())

            if _save_intermediates:
                df_Markers_gen.to_csv(_out + ".HLA_{0}.gen.markers.filt.indel.txt".format(_hla), sep='\t', header=True, index=True)




        # File Writing
        if _save_intermediates:

            if not _no_Indel:
                df_Markers_gen.to_csv(_out + ".HLA_{0}.gen.markers.indel.txt".format(_hla), sep='\t', header=True, index=True) # ***
            else:
                df_Markers_gen.to_csv(_out + ".HLA_{0}.gen.markers.noindel.txt".format(_hla), sep='\t', header=True, index=True)  # ***


        ### precursor which will be used in making map file.
        precursor_SNPS_forMAP = df_Markers_gen.columns.to_frame(index=False)
        print("\nframes for .map file.\n")
        print(precursor_SNPS_forMAP.head())


        # File Writing
        if _save_intermediates:

            if not _no_Indel:
                precursor_SNPS_forMAP.to_csv(_out + ".HLA_{0}.gen.precursor_SNPS_forMAP.indel.txt".format(_hla), sep='\t', header=True, index=False) # ***
            else:
                precursor_SNPS_forMAP.to_csv(_out + ".HLA_{0}.gen.precursor_SNPS_forMAP.noindel.txt".format(_hla), sep='\t', header=True, index=False)  # ***



        # Final output as Seqs
        df_Seqs_gen = df_Markers_gen.apply(lambda x : ''.join(x.astype(str)), axis=1) # *** Final dictionary


        # File Writing
        if _save_intermediates:

            if not _no_Indel:
                df_Seqs_gen.to_csv(_out+".HLA_{0}.gen.seqs.indel.txt".format(_hla), sep='\t', header=False, index=True) # ***
            else:
                df_Seqs_gen.to_csv(_out+".HLA_{0}.gen.seqs.noindel.txt".format(_hla), sep='\t', header=False, index=True) # ***


        print("\ndf_Markers_gen\n")
        print(df_Seqs_gen.head())






    if _3_LOADING_NUC:


        ##### Loading *_nuc.txt file.

        print("\n\n=============<NUC>=============\n")
        print(std_MAIN_PROCESS_NAME + "Processing *_nuc.txt file.\n")

        # Raw Markers
        df_raw_Markers_nuc, start_offset_nuc = LoadRawSeq2(_nuc, _hla, "nuc")
        # df_raw_Markers_nuc.to_csv(_out + ".HLA_{0}.nuc.raw.markers.txt".format(_hla), sep='\t', header=False, index=True)
        # print(df_raw_Markers_nuc.head())

        # Raw Seqeunces
        df_raw_Seqs_nuc = get_DF_rawSeqs(df_raw_Markers_nuc)
        # df_raw_Seqs_nuc.to_csv(_out + ".HLA_{0}.nuc.raw.seqs.txt".format(_hla), sep='\t', header=False, index=True)
        # print(df_raw_Seqs_nuc.head())

        # Raw Seqeunces splitted (for Exon checking)
        df_raw_Seqs_splitted_nuc = df_raw_Seqs_nuc.str.split('\|', expand=True)

        df_raw_Seqs_splitted_nuc.columns = pd.Index(["exon" + str(i + 1) for i in range(0, df_raw_Seqs_splitted_nuc.shape[1])])
        df_raw_Seqs_splitted_nuc.index.name = "Alleles"
        # df_raw_Seqs_splitted_nuc.to_csv(_out + ".HLA_{0}.nuc.raw.seqs.splitted.txt".format(_hla), sep='\t', header=True, index=True)

        print("\nRaw seqs(with Indels) of nuc splitted.\n")
        print(df_raw_Seqs_splitted_nuc.head())

        # (Trimmed) Sequences splitted without any indel character(for Exon checking with gen file)
        df_Seqs_splited_noIndel_nuc = df_raw_Seqs_splitted_nuc.apply(lambda x: ProcessIndel(x, _remove_indel=True), axis=0)

        if _hla == "C":
            # Hardcoding only for HLA-C due to exceptional form.
            df_Seqs_splited_noIndel_nuc = HardCodingforHLA_C(df_Seqs_splited_noIndel_nuc)

        print("\n\nSeqeunces(without Indels) of nuc splitted.\n")
        print(df_Seqs_splited_noIndel_nuc.head())


        # File Writing
        if _save_intermediates:
            df_Seqs_splited_noIndel_nuc.to_csv(_out + ".HLA_{0}.nuc.seqs.noindel.splitted.txt".format(_hla), sep='\t', header=True, index=True)  # Exporting 2 for nuc

        # No more processing to *_nuc.txt file is needed.




    if _4_LOADING_PROT:


        print("\n\n=============<Prot>=============\n")
        print(std_MAIN_PROCESS_NAME + "Processing Prot files for AA.\n")

        # Raw Markers
        df_raw_Markers_prot, dummy = LoadRawSeq2(_prot, _hla, "prot")
        print(df_raw_Markers_prot.head())

        # Raw Sequences
        sr_raw_Seqs_prot = get_DF_rawSeqs(df_raw_Markers_prot)
        print(sr_raw_Seqs_prot.head())




    if _5_GENERATING_AA_DICT:


        if _no_Indel:

            df_Seqs_IndelProcessed_prot = ProcessIndel(sr_raw_Seqs_prot, _remove_indel=True, _hla=("DQB1" if _hla == "DQB1" else "Not_given"))

            # File Writing
            if _save_intermediates:
                df_Seqs_IndelProcessed_prot.to_csv(_out + ".HLA_{0}.prot.seqs.noindel.txt".format(_hla), sep='\t', header=False, index=True)

        else:

            df_Seqs_IndelProcessed_prot = ProcessIndel(sr_raw_Seqs_prot, _remove_indel=False, _hla=("DQB1" if _hla == "DQB1" else "Not_given"))

            # File Writing
            if _save_intermediates:
                df_Seqs_IndelProcessed_prot.to_csv(_out + ".HLA_{0}.prot.seqs.indel.txt".format(_hla), sep='\t', header=False, index=True)


        ### Making final outputs.

        # Making columns for final outputs.
        l_relative_positions, l_genomic_positions, l_type = getPositionInfo_AA(df_Seqs_IndelProcessed_prot.iat[0], isREVERSE[_hla], int(start_offset_nuc),
                                                                               _has_Indel=(not _no_Indel), _df_precursor_map=precursor_AA_forMAP)

        print("\ngenomic positions")
        print(l_genomic_positions)
        print("\n\nrelative positions")
        print(l_relative_positions)
        print("\n\ntype information.")
        print(l_type)


        l_Markers = df_Seqs_IndelProcessed_prot.apply(lambda x : tuple(x)).tolist()

        df_Markers_prot = pd.DataFrame(l_Markers, index=df_Seqs_IndelProcessed_prot.index,
                                       columns=pd.MultiIndex.from_arrays([l_relative_positions, l_genomic_positions, l_type],
                                                                         names=["AA_rel_pos", "AA_gen_pos", "Type"]))

        precursor2_AA_forMAP = df_Markers_prot.columns.to_frame(index=False)

        print("\nFinal output as markers.\n")
        print(df_Markers_prot.head())
        print(precursor2_AA_forMAP.head())

        # File Writing
        if _save_intermediates:
            if _no_Indel:
                df_Markers_prot.to_csv(_out + ".HLA_{0}.prot.markers.noindel.txt".format(_hla), sep='\t', header=True, index=True)
                precursor2_AA_forMAP.to_csv(_out + ".HLA_{0}.prot.precursor2_AA_forMAP.noindel.txt".format(_hla), sep='\t', header=True, index=False)
            else:
                df_Markers_prot.to_csv(_out + ".HLA_{0}.prot.markers.indel.txt".format(_hla), sep='\t', header=True, index=True)
                precursor2_AA_forMAP.to_csv(_out + ".HLA_{0}.prot.precursor2_AA_forMAP.indel.txt".format(_hla), sep='\t', header=True, index=False)



        # Maptable in Heatmap. (2018. 10. 26.)
        __MAPTABLE__ = os.path.join(INTERMEDIATE_PATH, "HLA_MAPTABLE_{0}.hg{1}.imgt{2}.txt".format(_hla, _hg, _imgt))
        df_Markers_prot.to_csv(__MAPTABLE__, sep='\t', header=True, index=True)






    if _6_GENERATING_forMAP:

        print(std_MAIN_PROCESS_NAME + "Generating precursor file for *.map file.")

        """
        
        After INDEL('i") markers are processed.
        
        1. AA_forMAP
        2. precursor_AA_forMAP
        3. precurosr_SNPS_forMAP
        
        `SNPS_forMAP` will be made by merging those DataFrames sequentially.
        """

        # (2019. 01. 22.)
        if not _no_prime:

            final_AA_forMAP, final_SNPS_forMAP = coatingPrecursor(precursor_AA_forMAP.astype(str), precursor2_AA_forMAP.astype(str), precursor_SNPS_forMAP.astype(str))

        else:

            setPositionOfIndel(precursor2_AA_forMAP)
            setPositionOfIndel(precursor_SNPS_forMAP)

            final_AA_forMAP = precursor2_AA_forMAP
            final_SNPS_forMAP = precursor_SNPS_forMAP


        print("\nFinally created \"forMAP\" of AA.\n")
        print(final_AA_forMAP.head())

        print("\nFinally created \"forMAP\" of SNPS.\n")
        print(final_SNPS_forMAP.head())


        # File Writing
        if _save_intermediates:
            final_AA_forMAP.to_csv(_out + ".HLA_{0}.prot.forMAP.txt".format(_hla), sep='\t', header=True, index=False)
            final_SNPS_forMAP.to_csv(_out + ".HLA_{0}.gen.forMAP.txt".format(_hla), sep='\t', header=True, index=False)




    ### Complementing reverse strand.

    if isREVERSE[_hla]:

        """
        (2018. 9. 17.)
        I decided to complement reverse strands only in final output "HLA_DICTIONARY_{AA,SNPS}_hg{18,19,38}.{txt,map}".
        """

        df_Seqs_gen = df_Seqs_gen.apply(lambda x : ComplementStr(x))


        # File Writing
        if _save_intermediates:
            df_Seqs_gen.to_csv(_out + ".HLA_{0}.gen.complementSeqs.txt".format(_hla))



    return [df_Seqs_gen, df_Seqs_IndelProcessed_prot, final_SNPS_forMAP, final_AA_forMAP, __MAPTABLE__]

# end - ProcessIMGT()






#################### < Core Functions > ####################

### Main Module 1
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
def getIndelSpots(_conserved_1st_seq, _hla = "Not_given"):

    p = re.compile('\.+')

    iters = p.finditer(_conserved_1st_seq)
    l_INDELs = [i.span() for i in iters]

    if _hla == "DQB1":

        """
        (2018. 9. 16.)
        This exception handling was introduced becuase the 5th exon in DQB1 is either exon or indel itself.
        Also, this exceoption need to be gone through only when dealing with AA not SNPS(I hope so.).
        """

        l_INDELs.pop() # Abandon last indel assuming this will correspond to exon5.


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



def ProcessIndel(_sr, _remove_indel=False, _hla="Not_given"):

    # Acquiring span information using "conserved_1st_seq"
    l_spanInfo = getIndelSpots(_sr[0], _hla)

    return _sr.apply(lambda x: getTrimmedSeqs(x, l_spanInfo, _remove_indel))



### Main Module 5
def getPositionInfo_SNPS(_type, _df, _hla, _isReverse, _exon1_offset=1, _as_flattened=False, _has_Indel=False):

    if not (_type == "gen" or _type == "rel"):
        return -1

    if _type == "gen":

        offset_UTR = (_exon1_offset - 1) if not _isReverse else (_exon1_offset + 1)
        offset_Rest = _exon1_offset

    elif _type == "rel":

        offset_UTR = -1
        offset_Rest = 1


    # Preparing input dataframe
    df_UTR = _df.iloc[0, 0]  # Single String
    df_Rest = _df.iloc[0, 1:]  # Series

    # Final output to be returned.
    l_RETURN = []


    ### Processing UTR

    l_UTR = []
    N_indel = 0

    if _has_Indel:
        df_UTR = df_UTR[::-1]

        # Making Genomic Positions for UTR part.

    for i in range(0, len(df_UTR)):

        if (df_UTR[i] == 'z') or (df_UTR[i] == 'Z'):

            l_UTR.append('i')
            N_indel += 1

        else:

            if _type == "gen":
                l_UTR.append(offset_UTR + (-(i - N_indel) if not _isReverse else (i - N_indel)))

            elif _type == "rel":
                l_UTR.append(offset_UTR + (-(i - N_indel)))


    l_RETURN.append(l_UTR[::-1])


    ### Processing Rest of df

    curr_POS = offset_Rest

    for i in range(0, df_Rest.shape[0]):

        curr_string = df_Rest.iat[i]
        l_Rest = []
        N_indel = 0

        for j in range(0, len(curr_string)):

            if (curr_string[j] == 'z') or (curr_string[j] == 'Z'):

                l_Rest.append('i')
                N_indel += 1

            else:

                if _type == "gen":
                    l_Rest.append(curr_POS + ((j - N_indel) if not _isReverse else (-(j - N_indel))))

                elif _type == "rel":
                    l_Rest.append(curr_POS + (j - N_indel))


        # curr_POS = l_Rest[-1] + (+1 if _type=="rel" else -1)

        if _type == "rel":
            curr_POS = l_Rest[-1] +1

        elif _type == "gen":
            curr_POS = l_Rest[-1] + (+1 if not _isReverse else -1)

        l_RETURN.append(l_Rest)

    if _as_flattened:
        return [element for l in l_RETURN for element in l]
    else:
        return l_RETURN



### Main Module 5.5
def getPositionInfo_AA(_1st_string, _is_Reverse, _rel_s_offset, _has_Indel=False, _df_precursor_map="Not_given"):

    ### Checking necessary conditions. ###

    if not isinstance(_df_precursor_map, pd.DataFrame):
        return -1

    print(_df_precursor_map.head(10))
    print(_df_precursor_map.shape[0])



    ### Core variables. ###

    L = len(_1st_string)

    l_RETURN_rel_pos = []
    l_RETURN_gen_pos = []
    l_RETURN_type = []



    ##### Making relative positions. #####

    _idx_count = _rel_s_offset

    for i in range(0, L):

        if (_has_Indel) and ((_1st_string[i] == 'z') or (_1st_string[i] == 'Z')):
            l_RETURN_rel_pos.append('i')

        else:

            if _idx_count == 0:
                _idx_count += 1

            l_RETURN_rel_pos.append(_idx_count)
            _idx_count += 1



    ##### Making genomic positions and type information. #####

    N_indel = 0
    flag_aborted = False

    #     print("\nLength of first string is {0}\n".format(L))
    #     print(_1st_string)

    for i in range(0, L):

        if (_has_Indel) and ((_1st_string[i] == 'z') or (_1st_string[i] == 'Z')):

            l_RETURN_gen_pos.append('i')
            l_RETURN_type.append('i')
            N_indel += 1

        else:

            j = (3 * (i - N_indel) + 1)

            if j < _df_precursor_map.shape[0]:
                l_RETURN_gen_pos.append(int(_df_precursor_map.iat[j, 1]))
                l_RETURN_type.append(_df_precursor_map.iat[j, 2])
            else:

                L_l_RETURN = len(l_RETURN_gen_pos)

                Last_Element_gen = l_RETURN_gen_pos[-1]
                Last_Element_type = l_RETURN_type[-1]

                flag_aborted = True

                break

    # appending difference parts
    if flag_aborted:

        for i in range(0, (L - L_l_RETURN)):
            l_RETURN_gen_pos.append(Last_Element_gen + (3 * (i + 1) if not _is_Reverse else -3 * (i + 1)))
            l_RETURN_type.append(Last_Element_type)

    #         print("\nLength of l_RETURN is {0}\n".format(len(l_RETURN)))


    return [l_RETURN_rel_pos, l_RETURN_gen_pos, l_RETURN_type]



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

        df_temp.columns = pd.MultiIndex.from_arrays([idx_relative_pos, idx_genomic_pos, idx_TypeInfo])

        l_processed_sr.append(df_temp)

    df_RETURN = pd.concat(l_processed_sr, axis=1)
    df_RETURN.index.name = None
    df_RETURN.columns.names = ["SNP_rel_pos", "SNP_gen_pos", "Type"]

    return df_RETURN



### Module 7
def HardCodingforHLA_C(_df_C):

    # Additional Processing for nuc file of HLA-C.
    # _df_C := (trimmed) seqs / NoIndelChar / Splitted

    sr_temp = _df_C.iloc[:, -1]
    s_1st_seq = sr_temp.iloc[0]

    _exon_name = sr_temp.name

    p = re.compile('^[^x]+')
    s = p.match(s_1st_seq)


    my_span = s.span()

    sr_temp2 = sr_temp.apply(lambda x: x[my_span[0]:my_span[1]])
    #     print(sr_temp2.head())

    df_temp2 = pd.DataFrame(sr_temp2.tolist(), index=sr_temp2.index, columns=[_exon_name])
    #     print(df_RETURN.head())

    df_RETURN = pd.concat([_df_C.iloc[:, :-1], df_temp2], axis=1)

    return df_RETURN



### Module 8
def coatingPrecursor(_df_precursor1_AA, _df_precursor2_AA, _df_precursor1_SNPS):

    df_RETURN_AA_forMAP = None
    df_RETURN_SNPS_forMAP = None


    ##### (1) Marking "Signal Peptide" (only for `_df_precursor2_AA`)
    flag1 = _df_precursor2_AA.iloc[:, 0].str.match('^-')

    if flag1.any():

        sr_temp = _df_precursor2_AA.loc[flag1, "Type"].apply(lambda x: '_'.join([x, "Signal_Peptide"]))

        sr_Type = sr_temp.append(_df_precursor2_AA.loc[~flag1, "Type"]).reset_index(drop=True)
        sr_Type.name = "Type"

        t_df_precursor2_AA = pd.concat([_df_precursor2_AA.iloc[:, [0, 1]], sr_Type], axis=1)

    else:
        t_df_precursor2_AA = _df_precursor2_AA



    ##### (2) Re-labeling indel("i").

    # cf) `_df_precursor1_AA` doesn't have indels.

    ### <AA>

    flag_indel_AA = t_df_precursor2_AA.iloc[:, 0] == 'i'

    if flag_indel_AA.any():

        idx_indel_AA = t_df_precursor2_AA.loc[flag_indel_AA, :].index.tolist()

        for i in idx_indel_AA:
            t_rel_pos = 'x'.join([t_df_precursor2_AA.iat[i - 1, 0], t_df_precursor2_AA.iat[i + 1, 0]])
            t_gen_pos = round((int(t_df_precursor2_AA.iat[i - 1, 1]) + int(t_df_precursor2_AA.iat[i + 1, 1])) / 2)

            t_df_precursor2_AA.iat[i, 0] = t_rel_pos
            t_df_precursor2_AA.iat[i, 1] = t_gen_pos


    df_RETURN_AA_forMAP = t_df_precursor2_AA  # (***) "forMAP" file for AA is done.




    ### <SNPS>

    t_df_precursor1_SNPS = _df_precursor1_SNPS.copy()

    flag_indel_SNPS = t_df_precursor1_SNPS.iloc[:, 0] == 'i'

    if flag_indel_SNPS.any():

        idx_indel_SNPS = t_df_precursor1_SNPS.loc[flag_indel_SNPS, :].index.tolist()

        for i in idx_indel_SNPS:
            t_rel_pos = 'x'.join([t_df_precursor1_SNPS.iat[i - 1, 0], t_df_precursor1_SNPS.iat[i + 1, 0]])
            t_gen_pos = round((int(t_df_precursor1_SNPS.iat[i - 1, 1]) + int(t_df_precursor1_SNPS.iat[i + 1, 1])) / 2)

            t_df_precursor1_SNPS.iat[i, 0] = t_rel_pos
            t_df_precursor1_SNPS.iat[i, 1] = t_gen_pos


    ##### (3) Merging

    ### 1st Merging : `_df_precursor1_AA` (No indel) vs `t_df_precursor2_AA`

    df_merged1 = pd.merge(_df_precursor1_AA.astype(str), t_df_precursor2_AA.loc[~flag_indel_AA, :],
                          how='outer', left_on="SNP_gen_pos", right_on="AA_gen_pos")


    ## dividing `df_merged1` into "Matched" and "Unmatched" parts.
    flag_Unmatched = df_merged1.iloc[:, 0].isna()

    df_merged1_matched = df_merged1.loc[~flag_Unmatched, :].astype(str)
    df_merged1_unmatched = df_merged1.loc[flag_Unmatched, :].astype(str)

    ## Processing matched part.

    if df_merged1_matched.shape[0] % 3 != 0:
        print("\nThe # of rows of `df_merged1_matched` is not 3 times.\n")
        return -1

    for i in range(1, df_merged1_matched.shape[0], 3):
        t_AA_rel = df_merged1_matched.iat[i, 3]
        t_AA_gen = df_merged1_matched.iat[i, 4]
        t_AA_type = df_merged1_matched.iat[i, 5]

        df_merged1_matched.iat[i - 1, 3] = t_AA_rel
        df_merged1_matched.iat[i - 1, 4] = t_AA_gen + "_1st"
        df_merged1_matched.iat[i - 1, 5] = t_AA_type

        df_merged1_matched.iat[i, 4] = t_AA_gen + "_2nd"

        df_merged1_matched.iat[i + 1, 3] = t_AA_rel
        df_merged1_matched.iat[i + 1, 4] = t_AA_gen + "_3rd"
        df_merged1_matched.iat[i + 1, 5] = t_AA_type

    ## Processing Unmatched part.

    flag_3_prime = _df_precursor1_SNPS.iloc[:, 2].str.match("3_prime")
    sr_temp = df_merged1_unmatched.iloc[:, 5].apply(lambda x: '_or_'.join([x, "3_prime"]))
    sr_temp.name = "Type"
    df_merged1_unmatched = pd.concat([df_merged1_unmatched.iloc[:, [3, 4]], sr_temp], axis=1)


    ### 2nd Merging : `df_merged1_unmatched` vs. `_df_precursor1_SNPS`("3-prime" parts).

    df_merged2 = pd.merge(_df_precursor1_SNPS.loc[flag_3_prime, :].astype(str), df_merged1_unmatched,
                          how='outer', left_on="SNP_gen_pos", right_on="AA_gen_pos").iloc[
                 :df_merged1_unmatched.shape[0] * 3, :]

    for i in range(1, df_merged2.shape[0], 3):
        t_AA_rel = df_merged2.iat[i, 3]
        t_AA_gen = df_merged2.iat[i, 4]
        t_AA_type = df_merged2.iat[i, 5]

        df_merged2.iat[i - 1, 3] = t_AA_rel
        df_merged2.iat[i - 1, 4] = t_AA_gen + "_1st"
        df_merged2.iat[i - 1, 5] = t_AA_type

        df_merged2.iat[i, 4] = t_AA_gen + "_2nd"

        df_merged2.iat[i + 1, 3] = t_AA_rel
        df_merged2.iat[i + 1, 4] = t_AA_gen + "_3rd"
        df_merged2.iat[i + 1, 5] = t_AA_type


    ### 3rd(Last) Merging

    df_merged1_matched = pd.concat([df_merged1_matched, df_merged2]).reset_index(drop=True).iloc[:, [0, 3, 4, 5]]

    df_merged3 = pd.merge(t_df_precursor1_SNPS, df_merged1_matched,
                          how='left', left_on="SNP_rel_pos", right_on="SNP_rel_pos")

    df_RETURN_SNPS_forMAP = df_merged3  # (***) forMAP file for SNPS is done.


    return [df_RETURN_AA_forMAP, df_RETURN_SNPS_forMAP]




def setPositionOfIndel(_df_precursor):

    """

    """

    for i, value in _df_precursor.iloc[:, 0].iteritems():

        if value == 'i':

            t_rel_pos = 'x'.join([str(_df_precursor.iat[i-1, 0]), str(_df_precursor.iat[i+1, 0])])
            t_gen_pos = int(round((_df_precursor.iat[i-1, 1] + _df_precursor.iat[i+1, 1]) / 2))

            _df_precursor.iat[i, 0] = t_rel_pos
            _df_precursor.iat[i, 1] = t_gen_pos

    return 0



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

        ProcessIMGT.py

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

    parser.add_argument("-hg", help="\nHuman Genome Version.\n\n", required=True, choices=["18", "19", "38"])
    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="imgt_version", required=True)

    parser.add_argument("-nuc", help="\nInput *_nuc.txt file.\n\n", required=True)
    parser.add_argument("-gen", help="\nInput *_gen.txt file.\n\n", default="Not_given")
    parser.add_argument("-prot", help="\nInput *_prot.txt file.\n\n", default="Not_given")


    # Optional arguments
    parser.add_argument("--no-Indel", "-nI", help="\nNo Indels in  Markers.\n\n", action="store_true")




    ##### < for Test > #####

    # HLA-A / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "A",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt"
    #                           ])

    # HLA-A / No Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "A",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/A_prot.txt",
    #                           "-nI"
    #                           ])

    # HLA-C / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "C",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_prot.txt",
    #                           ])

    # HLA-C / No Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "C",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_prot.txt",
    #                           "-nI"
    #                           ])

    # HLA-B / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "B",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/B_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/B_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/B_prot.txt",
    #                           ])

    # HLA-B / No Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "B",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/C_prot.txt",
    #                           "-nI"
    #                           ])

    # HLA-DRB1 / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "DRB1",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DRB_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DRB_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DRB_prot.txt"
    #                           ])

    # HLA-DRB1 / No Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                           "-HLA", "DRB1",
    #                           "-imgt", "370",
    #                           "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DRB_nuc.txt",
    #                           "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DRB_gen.txt",
    #                           "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DRB_prot.txt",
    #                           "-nI"
    #                           ])

    # HLA-DQA1 / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                            "-HLA", "DQA1",
    #                            "-imgt", "370",
    #                            "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                            "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DQA_nuc.txt",
    #                            "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DQA_gen.txt",
    #                            "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DQA_prot.txt"
    #                           ])

    # HLA-DQB1 / Indel / imgt370
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                            "-HLA", "DQB1",
    #                            "-imgt", "370",
    #                            "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                            "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DQB_nuc.txt",
    #                            "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DQB_gen.txt",
    #                            "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DQB_prot.txt"
    #                           ])

    # HLA-DQB1 / Indel / imgt3320
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2.imgt3320",
    #                            "-HLA", "DQB1",
    #                            "-imgt", "3320",
    #                            "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                            "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA3320/alignments/DQB1_nuc.txt",
    #                            "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA3320/alignments/DQB1_gen.txt",
    #                            "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA3320/alignments/DQB1_prot.txt"
    #                           ])

    # HLA-DPA1 / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                            "-HLA", "DPA1",
    #                            "-imgt", "370",
    #                            "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                            "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DPA_nuc.txt",
    #                            "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DPA_gen.txt",
    #                            "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DPA_prot.txt"
    #                           ])

    # HLA-DPB1 / Indel
    # args = parser.parse_args(["-o", "/Users/wansun/Git_Projects/HATK/MkDict_v2/Dictionary.v2",
    #                            "-HLA", "DPB1",
    #                            "-imgt", "370",
    #                            "-hg", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/MakeDictionary/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                            "-nuc", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DPB_nuc.txt",
    #                            "-gen", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DPB_gen.txt",
    #                            "-prot", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/IMGTHLA/IMGTHLA370/alignments/DPB_prot.txt"
    #                           ])


    ##### < for Publish > #####
    args = parser.parse_args()

    # print(args)


    # main function execution
    ProcessIMGT(_out=args.o, _hla=args.HLA, _hg=args.hg, _imgt=args.imgt, _nuc=args.nuc, _gen=args.gen, _prot=args.prot,
                _no_Indel=args.no_Indel)
