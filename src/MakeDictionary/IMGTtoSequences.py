# -*- coding: utf-8 -*-

import os, sys, subprocess
import pandas as pd
import argparse, textwrap
import re

def IMGTtoSequences(_inputfile, _OUTPUT, _TYPE, _HLA, _HG=19,
                    _p_HLA_INTEGRATED_POSITIONS="data/HLA_INTEGRATED_POSITIONS_hg19.txt", _return_as_dataframe=False):

    """

    """

    ########## < Core Variables > ##########


    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting IMGTtoSequences.")


    ### General
    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}


    ### MARKING_CHUNKS & Insertions.
    line_number_marks = []
    Sequence_Chunks = [] # list of Dictionarys (which contains chunks of raw sequences).

    start_p_AA = {'A': -24, 'C': -24, 'B': -24, 'DRB1': -29, 'DQA1': -23, 'DQB1': -32, 'DPA1': -31, 'DPB1': -29}
    start_p_SNPS = {'A': -300, 'C': -283, 'B': -284, 'DRB1': -599, 'DQA1': -746, 'DQB1': -525, 'DPA1': -523, 'DPB1': -366}

    insertions = []
    insertions_idx = []
    insertions_len = []

    # insertions_idx_shifted = [] # (2018.1.22)

    trimmed_len = 0


    ### Standard Sequence.
    conserved_1st_seq = ''              # This includes substring '...'(insertions) and '|'.
    trimmed_conserved_1st_seq = ''      # This includes only '|'.


    ### (2018/1/17) HLA position information(exon, intron, etc.)
    HLA_INTEGRATED_POSITIONS_hg = pd.read_table(_p_HLA_INTEGRATED_POSITIONS, sep='\t', header=None,
                                                usecols = [0,1,2,3],
                                                names=['HLA_name', 's_pos', 'e_pos', 'Class'],
                                                index_col=[0, 3])
    # print(HLA_INTEGRATED_POSITIONS_hg.head())


    ### HLA name
    HLA_name = _HLA
    # print(std_MAIN_PROCESS_NAME + "Given HLA gene : %s" % (HLA_name))


    ### DataFrame for final output.
    df_forRETURN = pd.DataFrame()



    ########## < Argument Checking. > ##########

    # Checking "_p_HLA_INTEGRATED_POSITIONS"

    # Preparing intermediate paths.
    _OUTPUT = _OUTPUT if not _OUTPUT.endswith('/') else _OUTPUT.rstrip('/')

    INTERMEDIATE_PATH = os.path.dirname(_OUTPUT)

    if not os.path.exists(INTERMEDIATE_PATH):
        os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))



    ########## < Control Flags > ##########

    # (2018.2.11) Segmentation of Main processes.

    MARKING_CHUNKS = 1
    READING_RAW_SEQUENCES = 1
    APPENDING_RAW_SEQUENCES = 1
    UNTRIMMED_MARKERS = 1
    UNTRIMMED_SEQUENCES = 1
    PROCESSING_INSERTIONS = 1
    TRIMMED_SEQUENCES = 1
    LOAD_NUC_FILE = 1  # (2018.2.13)
    TRIMMED_MARKERS = 1
    ENCODE_INSERTIONS = 1





    if MARKING_CHUNKS:

        ########## < 1. Finding Segments beforehand > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Finding Segments beforehand.")

        # Before reading "A_prot.txt" or "A_gen.txt" file, get the line numbers where chunks are made by using bash "awk" and "grep".
        # These line numbers will be stored in `line_number_marks'

        if _TYPE == 'AA':
            line_number_marks = subprocess.check_output(' '.join(["awk", "'{print $1}'", _inputfile, "|", "grep -n Prot"]), shell=True)
        elif _TYPE == 'SNPS':
            line_number_marks = subprocess.check_output(' '.join(["awk", "'{print $1}'", _inputfile, "|", "grep -n gDNA"]), shell=True)

        line_number_marks = line_number_marks.decode('UTF-8').rstrip('\n').split('\n')
        line_number_marks = [int(item.split(':')[0])-1 for item in line_number_marks]
        # Unix uses 1-base while python uses 0-base. That's why I subtract 1 for each acquired line numbers.

        print("\nAcquired line numbers :")
        print(line_number_marks)


        # Counting how many HLA alleles are there.
        try:
            # ex. [6, 3613, 7220, 10827, 14434] => ((3613-2)-(6+1) == ((7220-2)-(3613+1)) == ... == 3604
            # Actually, this approach is very fragile because slightly changed contents of the file could cause a disaster.
            # But, I'll assum that basic structure of the files won't change.
            # (2017.12.29) In "*_gen.txt" file, offset of -3 is the last line.

            f_keys_start = line_number_marks[0]+2
            f_keys_end = line_number_marks[1]-2 if _TYPE == 'AA' else line_number_marks[1] - 3 if _TYPE == 'SNPS' else None

            NumberofKeys = f_keys_end - (f_keys_start - 1)
            print("\n# of keys(== # of HLA alleles) is : %d" % (NumberofKeys))

        except KeyError:
            print(std_ERROR_MAIN_PROCESS_NAME + "Length of marks is less than 2. Someting Wrong!")
            sys.exit()



    if READING_RAW_SEQUENCES:

        ########## < 2. Reading input_file > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Reading Input files.")

        # As the information of the raw files are prepared, It's time to read the files through the filestream.
        f = open(_inputfile)
        f_lines = f.readlines()
        f.close()

        # "Please see http://hla.alleles.org/terms.html for terms of use." => This line is the cap for the end.
        # print(f_lines[len(f_lines)-2])
        line_number_marks.append(len(f_lines)-2)


    if APPENDING_RAW_SEQUENCES:

        ########## < 3. Making respective Dictionarys > ##########

        print(std_MAIN_PROCESS_NAME + "[3] Making respective Dictionaries.")

        # Chunks that start with "Prot" or "gDNA" will be contained in `Sequence_Chunks` dictionary.

        # I found that the number of HLA alleles could be different for each chunks.
        # So, I decided to introduce list of dictionary.

        print("\n\n")
        print("First line index is {0}\n\t{1}".format(f_keys_start, f_lines[f_keys_start]))
        print("Last line index is {0}\n\t{1}".format(f_keys_end, f_lines[f_keys_end]))

        print("\nLength of keys is %d" % (f_keys_end - (f_keys_start - 1)))  # length of keys



        # The number of HLA alleles in first chunks must be the maximum number of the keys.
        # Ex) [6, 3613, 7220, 10827, 14434] => (3613-2)-(6+1) == maximum # of keys

        # So, Make 1st key of dictionary with (marks[1] - marks[0]).
        # After that, keep checking whether the number of next key of dictionary is smaller than that of 1st dictionary.
        # Assumption that the number of chunks is bigger than 2.

        # f_keys_start : 1st line of sequence chunk contents(A*01:01:01:01) / f_keys_end : last line(A*80:02).


        ### Process 1st chunk first, then store it in `dict_first`.

        # k_first := ['A*01:01:01:01', 'A*01:01:01:02N', ... , 'A*80:02']
        k_first = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(f_keys_start, f_keys_end+1)]

        print("Length of k_first is {0}\n\n".format(len(k_first)))
        # print(k_first)

        dict_first = {}

        for i in range(0, len(k_first)):

            s = re.split('\s+', f_lines[f_keys_start + i].lstrip(' ').rstrip(' \n'))[1:]
            dict_first[k_first[i]] = ''.join(s)

        # for k,v in dict_first.items():
        #     print("key and value are : %s\t\t%s" % (k, v))

        Sequence_Chunks.append(dict_first)

        # So far, the 1st chunk is made and prepared in `dict_fisrt`.



        ### From here the rest of dictionaries are to be processed. (2nd dictionary ~ The last dictionary)
        for i in range(1, len(line_number_marks)-1): # Note that the start of the for loop is 1 not 0.

            if _TYPE == 'AA':
                t_keys = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(line_number_marks[i] + 2, line_number_marks[i + 1] - 2 + 1)]
            elif _TYPE == 'SNPS':
                t_keys = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(line_number_marks[i] + 2, line_number_marks[i + 1] - 3 + 1)] # 이거 -3 하나때문에 if-else추가해서 예외처리함

            t_s = [''.join(re.split('\s+', f_lines[(line_number_marks[i]+2) + k].lstrip(' ').rstrip(' \n'))[1:]) for k in range(0, len(t_keys))]

            Sequence_Chunks.append({k : v for k,v in zip(t_keys, t_s)})


        # All dictionaries have been processed.


        ########## < 4. Merging respective Dictionarys > ##########

        print(std_MAIN_PROCESS_NAME + "[4] Merging respective Dictionaries")

        # Merge all dictionaries based on the 1st dictionary.

        for i in range(0, len(k_first)):

            # print("key %d: %s" % (i, k_first[i]))

            # t_string = ""
            #
            # for j in range(0, len(Prots)):
            #
            #     if k_first[i] in Prots[j].keys():
            #         t_string += Prots[j][k_first[i]]


            # string concatenation
            t_string = ''.join([Sequence_Chunks[j][k_first[i]] for j in range(0, len(Sequence_Chunks)) if k_first[i] in Sequence_Chunks[j].keys()])

            dict_first[k_first[i]] = tuple(t_string)

        print("\nConcatenation is Done.")


    if UNTRIMMED_MARKERS:

        ########## < 5. Dictionary to DataFrame(untrimmed_markers) > ##########

        print(std_MAIN_PROCESS_NAME + "[5] Dictionary to DataFrame(untrimmed_markers)\n")

        df_MARKERS = pd.DataFrame.from_dict(dict_first, orient='index')
        df_MARKERS.fillna(value='x', inplace=True)


        ##### < Removing vertical bars('|') > #####

        conserved_1st_seq = df_MARKERS.iloc[0, :]                       # Bringing standard seq(1st seq.)
        vb_here = conserved_1st_seq[conserved_1st_seq == '|'].index     # Acquiring the indexes of '|'
        df_MARKERS = df_MARKERS.drop(vb_here, axis=1)                   # Drop those columns
        df_MARKERS.columns = pd.Index(range(0, df_MARKERS.shape[1]))    # Initialize the columns information.


        ##### < Replacing '-' characters with the equivalent of 1st sequence> #####

        conserved_1st_seq = df_MARKERS.iloc[0, :]

        print(df_MARKERS.head())

        for i in range(0, len(conserved_1st_seq)):

            if conserved_1st_seq[i] != '.':

                col = df_MARKERS.iloc[1:, i]
                col[col == '-'] = conserved_1st_seq[i]
                col[col == '*'] = 'x'

        if _TYPE == 'AA':
            df_MARKERS.to_csv(_OUTPUT + '.AA.markers.txt', sep='\t', header=False, index=True)
        elif _TYPE == 'SNPS':
            df_MARKERS.to_csv(_OUTPUT + '.SNPS.markers.txt', sep='\t', header=False, index=True)



    if UNTRIMMED_SEQUENCES:

        ########## < Untrimmed Sequences > ##########

        df_Seqs = pd.DataFrame(df_MARKERS.apply(lambda x : ''.join(x.astype(str)), axis=1)) # Aggregating columns to a string.

        # HLA Sequences in `df_Seqs` still have insertions.
        # These intertions will be processed by regular expressions.


        if _TYPE == 'AA':
            df_Seqs.to_csv(_OUTPUT + '.AA.seqs.txt', sep='\t', header=False)
        elif _TYPE == 'SNPS':
            df_Seqs.to_csv(_OUTPUT + '.SNPS.seqs.txt', sep='\t', header=False)


    if PROCESSING_INSERTIONS:

        """
        (2018.2.12) (1) PROCESSING_INSERTIONS, (2) TRIMMED_SEQUENCES => These two parts will be worked on the git branches.
        
        - The substring '....' will be substituted by character 'i'.
        - On the second thought, Just substitue '....' part and Do the rest of processing in 'TRIMMED_SEQUENCES' block.
                
        """

        ########## < 6. Processing Insertions > ##########

        print(std_MAIN_PROCESS_NAME + "[6] Processing Insertions.")

        insertions_idx = [(obj.start(0), obj.end(0)) for obj in re.finditer('\.+', ''.join(conserved_1st_seq))]

        if len(insertions_idx) > 0:

            insertions_len = [(idx[1] - idx[0]) for idx in insertions_idx]

            # Insertion processing on `conserved_1st_seq` is done.

            # Extracting values of insertion markers.
            insertions = [[df_Seqs.iat[i, 0][idx[0] : idx[1]] for i in range(0, len(df_Seqs))] for idx in insertions_idx]

            df_Insertions = pd.DataFrame(insertions).transpose() # This DataFrame will be used with final output.
            df_Insertions.index = pd.Index(df_Seqs.index)

            if _TYPE == 'AA':
                df_Insertions.to_csv(_OUTPUT+'.AA.insertions.txt', sep='\t', header=False, index=True)
            elif _TYPE == 'SNPS':
                df_Insertions.to_csv(_OUTPUT+'.SNPS.insertions.txt', sep='\t', header=False, index=True)

            # (2018. 8. 13.) Thsese ".insertions.txt" files are not necessarily useful. (Maybe should be excluded.).


    if TRIMMED_SEQUENCES:

        ########## < 7. Trimming out insertion sequences. > ##########

        print(std_MAIN_PROCESS_NAME + "[7] Trimming out insertion sequences.")

        # Substitute insertions sequences (ex. "......") to character 'i'.

        df_Seqs_trimmed = []

        for i in range(0, len(df_Seqs)):

            t_string = df_Seqs.iat[i, 0]

            shifted = 0

            for j in range(0, len(insertions_idx)):

                if j > 0:
                    # Applying shifted situation.
                    shifted += insertions_len[j-1] -1

                _idx_start = insertions_idx[j][0] - shifted
                _idx_end = insertions_idx[j][1] - shifted

                # Substituting part. (ex. "......." => 'i')
                t_string = t_string[0: _idx_start] + 'i' + t_string[_idx_end: ]


            df_Seqs_trimmed.append(t_string)

        df_Seqs_trimmed = pd.DataFrame(df_Seqs_trimmed)
        df_Seqs_trimmed.index = df_Seqs.index


        # (2018.3.15) Gather the dataframes which will be returned as `df_forRETURN`.
        df_forRETURN = df_Seqs_trimmed

        if _TYPE == 'AA':
            df_Seqs_trimmed.to_csv(_OUTPUT+'.AA.seqs.trim.txt', sep='\t', header=False)

        elif _TYPE == 'SNPS':
            df_Seqs_trimmed.to_csv(_OUTPUT+'.SNPS.seqs.trim.txt', sep='\t', header=False)

        # The product "*.AA.seqs.trim.txt" also could be unuseful.



    if LOAD_NUC_FILE:

        ########## < 7.5. Loading and Processing *_nuc.txt file > ##########

        if _TYPE == 'AA':

            print(std_MAIN_PROCESS_NAME + "[7.5] Loading and Processing *_nuc.txt file.")

            # Going to use "*_nuc.txt" file too.
            # Conduct the jobs done so far to "*_nuc.txt" file one more. (Only in case of AA.)

            ##### > ====

            # Initializing the variable used before.
            line_number_marks = None
            Sequence_Chunks = []

            f = None
            k_first = None

            dict_first = None

            ##### < Chunk indexing > #####

            _nuc_filename = re.sub(pattern='_prot', repl='_nuc', string=_inputfile)
            print("\n\nLoading \"{0}\" file".format(_nuc_filename))

            line_number_marks = subprocess.check_output(' '.join(["awk", "'{print $1}'", _nuc_filename, "|", "grep -n cDNA"]), shell=True)
            line_number_marks = line_number_marks.decode('UTF-8').rstrip('\n').split('\n')
            line_number_marks = [int(item.split(':')[0]) - 1 for item in line_number_marks]

            print("\nMarks is {0}".format(line_number_marks))

            allele_start = line_number_marks[0] + 3
            allele_end = line_number_marks[1] - 2
            NumberOfAlleles = allele_end - (allele_start - 1)

            print("\nallele start and end are {0} and {1}\n".format(allele_start, allele_end))
            print("# of Alleles is {0}\n".format(NumberOfAlleles))


            ##### < Reading Input files > #####
            f = open(_nuc_filename)
            f_lines = f.readlines()
            f.close()

            # print(f_lines[-2]) # "Please see http://hla.alleles.org/terms.html for terms of use."
            line_number_marks.append(len(f_lines) - 2)


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


            ##### < Processing each Chunks > #####

            for i in range(1, len(line_number_marks) - 1):
                t_keys = [re.split('\s+', f_lines[j].lstrip(' ').rstrip(' \n'))[0] for j in range(line_number_marks[i] + 3, line_number_marks[i + 1] - 2 + 1)]
                t_s = [''.join(re.split('\s+', f_lines[(line_number_marks[i] + 3) + k].lstrip(' ').rstrip(' \n'))[1:]) for k in range(0, len(t_keys))]

                Sequence_Chunks.append({k: v for k, v in zip(t_keys, t_s)})


            ##### < Merging all chunks > #####

            for i in range(0, len(k_first)):
                t_string = ''.join([Sequence_Chunks[j][k_first[i]] for j in range(0, len(Sequence_Chunks)) if k_first[i] in Sequence_Chunks[j].keys()])
                dict_first[k_first[i]] = tuple(t_string)

            # Unprocessed DataFrame of codon information.
            df_CODON = pd.DataFrame.from_dict(dict_first, orient='index').fillna('x')
            df_CODON.to_csv(_OUTPUT+'.codon.txt', sep='\t', header=False, index=True)
            print(df_CODON.head())


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

            # (2)

            conserved_1st_seq = df_CODON.iloc[0, :]

            fordrop = conserved_1st_seq[conserved_1st_seq == '.'].index.tolist()

            df_CODON_trimmed = df_CODON.drop(fordrop, axis=1)
            df_CODON_trimmed.columns = pd.Index(range(0, df_CODON_trimmed.shape[1]))

            # for j in range(0, df_CODON_trimmed.shape[1]):
            #     col = df_CODON_trimmed.iloc[1:, j]
            #     col[col == '-'] = df_CODON_trimmed.iat[0, j]
            #     col[col == '*'] = 'x'

            # DataFrame where insertions('....') are removed.
            df_CODON_trimmed.to_csv(_OUTPUT+'.codon.trim.txt', sep='\t', header=None, index=True)

            # (3)

            # DataFrame of condon sequence.
            df_CODON_seq = df_CODON_trimmed.apply(lambda x: ''.join(x.astype(str)), axis=1)
            df_CODON_seq.to_csv(_OUTPUT+'.codon.seq.txt', sep='\t', header=None, index=True)


            dict_CODON_seq = df_CODON_seq.apply(lambda x : x.split('|')).to_dict()

            # (4)

            # DataFrame splitted by exons.
            df_CODON_seq_splited = pd.DataFrame.from_dict(dict_CODON_seq, orient='index')
            df_CODON_seq_splited.to_csv(_OUTPUT+'.codon.exon.txt', sep='\t', header=False, index=True)


            # (5)

            # Genomic Positions information will be allocated to `df_CODON_seq_splitted` using "HLA_INTEGRATED_POSITIONS_hg.txt".

            # if len(name_check) == 1:
            #     name_check = name_check.iat[-1]
            # else:
            #     print("HLA name got wrong")
            #     sys.exit()

            print("\nCurrent file's HLA type is %s" % HLA_name)

            HLA_EXON_POS = HLA_INTEGRATED_POSITIONS_hg.loc[HLA_name].filter(regex='exon\d', axis=0)

            ### ***Most important condition***
            if df_CODON_seq_splited.shape[1] == HLA_EXON_POS.shape[0]:

                # This condition is the reason why "*_nuc.txt" file is used.
                # It is checking whether the number of exons is matched to that of the parts created by splitting sequence with '|'.

                df_CODON_seq_trimmed = df_CODON_seq_splited.apply(lambda x: ''.join(x.astype(str)), axis=1)
                # df_CODON_seq_trimmed.to_csv(_OUTPUT+'.codon.seq.trim.txt', sep='\t', header=False, index=True)
                df_CODON_marker_trimmed = pd.concat([pd.DataFrame([tuple(df_CODON_seq_trimmed.iat[i]) for i in range(0, df_CODON_seq_trimmed.shape[0])])], axis=1)


                genuine_contents = re.findall('\w+[^x]', df_CODON_seq_trimmed.iat[0])

                if len(genuine_contents) == 1:

                    # The length of `genuine_contents` must be 1.
                    genuine_contents = genuine_contents[0]

                    if len(genuine_contents)%3 != 0:

                        # Strange cases.
                        print('string of nuc file is wrong !')
                        sys.exit()

                else:
                    print('string of nuc file is wrong !')
                    sys.exit()


                genuine_contents_codon = [''.join([genuine_contents[i], genuine_contents[i+1], genuine_contents[i+2]]) for i in range(0, len(genuine_contents), 3)]

                ### Extracting genomic positions corresponding to exon parts in "HLA_INTEGRATED_POSITIONS_hg.txt"
                l_header = []

                if isREVERSE[HLA_name]:

                    for i in range(0, HLA_EXON_POS.shape[0]):
                        l_header.extend([z for z in range(HLA_EXON_POS.iat[i, 0], HLA_EXON_POS.iat[i, 1] - 1, -1)])

                else:

                    for i in range(0, HLA_EXON_POS.shape[0]):
                        l_header.extend([z for z in range(HLA_EXON_POS.iat[i, 0], HLA_EXON_POS.iat[i, 1] + 1)])


                genomic_position_by3 = [l_header[i] for i in range(1, len(l_header), 3)]


                # AA_HEADER = pd.Series(genuine_contents_codon, index=genomic_position_by3)
                # AA_HEADER.to_csv(_OUTPUT+'.HEADER.txt', sep='\t', header=False, index=True)

            else:
                print()


    if TRIMMED_MARKERS:

        ########## < 8. Trimmed Markers > ##########

        print(std_MAIN_PROCESS_NAME + "[8] Generating trimmed Markers.")

        # Make a new DataFrame where `df_Seqs` is divided by markers.

        df_MARKERS_trimmed = pd.concat([pd.DataFrame([tuple(df_Seqs_trimmed.iat[i, 0]) for i in range(0, len(df_Seqs_trimmed))])], axis=1)
        df_MARKERS_trimmed.index = df_MARKERS.index


        if _TYPE == 'AA':

            """
            I'll make a MultiIndex object with next three information.
            
            (1) genomic_position 
            (2) relative_position 
            (3) 3-base codon
            
            
            (2018.2.15)             
            I'll make a header with only codon contents except for the part where 'xxx..' appears in "*_nuc.txt" file.
            
            The problem related to the length will appear, so the next three jobs should be done as interating over `conserved_1st_seq` in "*.trim.marker".
            
                (1) Exclude insertion parts.
                (2) relative_position / Make relative_position.
                (3) Add a padding for the part where 'xxx' appears.
                
            
            """

            # Utilize `genuine_contents_codon` and `genomic_position_by3` made in previous "LOAD_NUC_FILE" code block.

            conserved_1st_seq = df_Seqs_trimmed.iat[0, 0]

            # Make a new list of "genomic_position" and "codon" due to the spot where 'i' makes a shift.

            LABEL_genomic_position = []
            LABEL_relative_position = []
            LABEL_codons = []

            count_indel = 0
            rel_start = start_p_AA[HLA_name]
            gen_last_position = genomic_position_by3[-1]

            for i in range(0, len(conserved_1st_seq)):

                if conserved_1st_seq[i] == 'i':
                    count_indel += 1

                    # Just append 'i' to the three lists when the value is 'i' in `convserved_1st_seq`.
                    LABEL_relative_position.append('i')
                    LABEL_genomic_position.append('i')
                    LABEL_codons.append('i')

                elif conserved_1st_seq[i] != 'i':

                    # relative_position
                    rel_pos = rel_start + (i - count_indel)
                    LABEL_relative_position.append(rel_pos if rel_pos < 0 else rel_pos + 1)

                    if (i - count_indel) < len(genuine_contents_codon):

                        # genomic_position
                        LABEL_genomic_position.append(genomic_position_by3[i - count_indel])

                        # codon
                        LABEL_codons.append(genuine_contents_codon[i - count_indel])

                    else:
                        # codon
                        LABEL_codons.append('xxx')

                        # genomic_position
                        gen_last_position += -3 if isREVERSE[HLA_name] else 3
                        LABEL_genomic_position.append(gen_last_position)



            df_MARKERS_trimmed.columns = pd.MultiIndex.from_arrays([LABEL_genomic_position, LABEL_relative_position, LABEL_codons], names=['genomic_position', 'relative_position' ,'codon'])
            print(df_MARKERS_trimmed.head())

            df_MARKERS_trimmed.to_csv(_OUTPUT+'.AA.markers.trim.labeled.txt', sep='\t', header=True, index=True)
            df_MARKERS_trimmed.iloc[0, :].transpose().to_csv(_OUTPUT+'.AA.forMAP.txt', sep='\t', index=True)


        elif _TYPE == 'SNPS':

            """
            In case of SNPS, a header will be made only with next two things. 
            
                (1) genomic_position
                (2) reltive_position
                
            """

            HLA_POS = HLA_INTEGRATED_POSITIONS_hg.loc[HLA_name]

            # print(HLA_POS.head())
            #
            # POS_start = HLA_POS.iat[0,0]
            # POS_end = HLA_POS.iat[HLA_POS.shape[0]-1, HLA_POS.shape[1]-1]
            #
            # print("POS_start : {0}\nPOS_end : {1}".format(POS_start, POS_end))

            conserved_1st_seq = df_Seqs_trimmed.iat[0, 0]

            # I need to make a new list for "genomic position" and "codon" due to the 'i'.

            LABEL_genomic_position = []
            LABEL_relative_position = []

            count_indel = 0
            gen_start = HLA_POS.iat[0, 0]
            rel_start = start_p_SNPS[HLA_name]

            for i in range(0, len(conserved_1st_seq)):

                if conserved_1st_seq[i] == 'i':
                    count_indel += 1

                    # Just append 'i' to the lists when the value of `conserved_1st_seq` is 'i'.
                    LABEL_relative_position.append('i')
                    LABEL_genomic_position.append('i')

                elif conserved_1st_seq[i] != 'i':

                    # relative_position
                    rel_pos = rel_start + (i - count_indel)
                    LABEL_relative_position.append(rel_pos if rel_pos < 0 else rel_pos + 1)

                    # genomic_position
                    gen_pos = gen_start + ((i - count_indel) if not isREVERSE[HLA_name] else (count_indel - i)) # -(i - count_indel) 을 더해주는거
                    LABEL_genomic_position.append(gen_pos)

            df_MARKERS_trimmed.columns = pd.MultiIndex.from_arrays([LABEL_genomic_position, LABEL_relative_position], names=['genomic_position', 'relative_position'])
            print(df_MARKERS_trimmed.head())

            df_MARKERS_trimmed.to_csv(_OUTPUT+'.SNPS.markers.trim.labeled.txt', sep='\t', header=True, index=True)
            df_MARKERS_trimmed.iloc[0, :].transpose().to_csv(_OUTPUT+'.SNPS.forMAP.txt', sep='\t', index=True)

    ### DataFrame which has marker information has been made.


    if ENCODE_INSERTIONS:

        ########## < 9. Encoding Insertions. > ##########

        print(std_MAIN_PROCESS_NAME + "[9] Encoding Insertions.")

        """
        (2018.2.18)        
        So far, Insertion or Deletion parts on the sequence were marked as only a single character 'i'.
        (Anyway, the marker labeled with 'i' is INDEL, not either Insertion or Deletion.)
        
        But, Professor Han requested me to introduce the way to distinguish given 'i' is insertion or deletion.
        So, I'll encode the substring of INDEL part to 'Z' or 'z'.  
        
        Z := "Insertion".
        z := "Deletion".
        
        """

        pattern_INSERTIONS = re.compile('\.+|\*+|x+')
        # Substrings which consists of these patterns will be classified to "Deletion" (not "Insertion").

        try:

            df_Insertions_encoded = df_Insertions.applymap(lambda x : 'z' if pattern_INSERTIONS.match(x) else 'Z')
            df_Insertions_encoded2 = df_Insertions.applymap(lambda x : '#' if pattern_INSERTIONS.match(x) else '@')

        except:

            """
            (2018.3.5) 
            If no INDEL exists in `conserved_1st_seq` of some HLA gene(ex. HLA-A), then just return `df_forRETURN`.
            Introduced `__return_as_dataframe` flag for this.
            
            """

            return df_forRETURN if _return_as_dataframe else 0

        print(df_Insertions_encoded.head(30))

        if _TYPE == 'AA':
            df_Insertions_encoded.to_csv(_OUTPUT+'.AA.insertions.encode.txt', sep='\t', header=False)
        elif _TYPE == 'SNPS':
            df_Insertions_encoded.to_csv(_OUTPUT+'.SNPS.insertions.encode.txt', sep='\t', header=False)



        # Processing the index of `df_MARKERS_trimmed` to extract the index of 'i' by using filter() and Regular Expressions.

        sr_temp = df_MARKERS_trimmed.columns.get_level_values(level='genomic_position').to_series().reset_index(drop=True)
        # No special reason to use `genomic_position`. `relative_position` also can be utilized to do this.

        HereINSERTIONS = sr_temp.apply(lambda x : True if re.match(pattern='i', string=str(x)) else False)

        idx_INSERTIONS_COLUMNS = HereINSERTIONS[HereINSERTIONS].index.tolist()
        print(idx_INSERTIONS_COLUMNS)

        print(df_MARKERS_trimmed.iloc[:, idx_INSERTIONS_COLUMNS].head())



        for i in range(0, df_Insertions_encoded.shape[1]):
            df_MARKERS_trimmed.iloc[:, idx_INSERTIONS_COLUMNS[i]] = df_Insertions_encoded.iloc[:, i]


        df_forRETURN = df_MARKERS_trimmed.apply(lambda x : ''.join(x.astype(str)), axis=1)


        if _TYPE == 'AA':
            df_MARKERS_trimmed.to_csv(_OUTPUT+'.AA.markers.trim.labeled.txt', sep='\t', header=True)
            df_forRETURN.to_csv(_OUTPUT+'.AA.seqs.trim.txt', sep='\t', header=False, index=True)
        elif _TYPE == 'SNPS':
            df_MARKERS_trimmed.to_csv(_OUTPUT+'.SNPS.markers.trim.labeled.txt', sep='\t', header=True)
            df_forRETURN.to_csv(_OUTPUT+'.SNPS.seqs.trim.txt', sep='\t', header=False, index=True)




        # If INDEL parts are marked as "P" or "A", it's quite hard to differentiate by eyes.
        # so, try marking them with 'B' and 'b' combination.
        for i in range(0, df_Insertions_encoded2.shape[1]):
            df_MARKERS_trimmed.iloc[:, idx_INSERTIONS_COLUMNS[i]] = df_Insertions_encoded2.iloc[:, i]

        if _TYPE == 'AA':
            df_MARKERS_trimmed.to_csv(_OUTPUT+'.AA.markers.trim.labeled.visual.txt', sep='\t', header=True)
            df_MARKERS_trimmed.apply(lambda x : ''.join(x.astype(str)), axis=1).to_csv(_OUTPUT+'.AA.seqs.trim.visual.txt', sep='\t', header=False, index=True)
        elif _TYPE == 'SNPS':
            df_MARKERS_trimmed.to_csv(_OUTPUT+'.SNPS.markers.trim.labeled.visual.txt', sep='\t', header=True)
            df_MARKERS_trimmed.apply(lambda x : ''.join(x.astype(str)), axis=1).to_csv(_OUTPUT+'.SNPS.seqs.trim.visual.txt', sep='\t', header=False, index=True)


        return df_forRETURN if _return_as_dataframe else 0



if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        IMGTtoSequences.py

        : Processing(Parsing HLA sequence information distributed by IMGT-HLA.


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
    parser.add_argument("-i", help="\nInput file name\n\n", required=True, metavar='INPUT')
    parser.add_argument("-o", help="\nOutput file name prefix\n\n", required=True, metavar='OUTPUT')

    parser.add_argument("-HLA", help="\nHLA gene name which you will process.\n\n", required=True, metavar='HLA',
                        choices=["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"])

    # parser.add_argument("-hg", help="Human Genome Version[18, 19 or 38]\n\n", required=True, metavar='hg', choices=["18", "19", "38"])
    parser.add_argument("--type", "-t", help="\nSequence type to deal with(Amino Acids[AA] or SNPs[SNPS]\n\n", required=True, choices=["AA", "SNPS"], metavar='TYPE')
    parser.add_argument("--HLA-pos", "-Hp", help="\nHLA gene position information table.\n\n", required=True)



    ##### < for Test > #####

    # # AA / imgt370 / hg18 / HLA-A
    # args = parser.parse_args(["-i", "../data/IMGTHLA370/alignments/A_prot.txt",
    #                           "-t", "AA",
    #                           "-HLA", "A",
    #                           "-Hp", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-o", "../TEST20180812/IMGT_A"
    #                           ])

    # # AA / imgt370 / hg18 / HLA-C
    # args = parser.parse_args(["-i", "../data/IMGTHLA370/alignments/C_prot.txt",
    #                           "-t", "AA",
    #                           "-HLA", "C",
    #                           "-Hp", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/HLA_INTEGRATED_POSITIONS_hg18.txt",
    #                           "-o", "../TEST20180812/IMGT_C"
    #                           ])

    # # AA / imgt3320 / hg19
    # args = parser.parse_args(["-i", "./data/IMGTHLA3320/alignments/A_prot.txt", "-o", "IMGT_A", "-t", "AA"])


    # SNPS / imgt370
    # args = parser.parse_args(["-i", "./data/IMGTHLA370/alignments/A_gen.txt", "-o", "IMGT_A", "--type", "SNPS"])



    ##### < for Publish > #####
    args = parser.parse_args()
    # print(args)


    # main function execution
    IMGTtoSequences(_inputfile=args.i, _OUTPUT=args.o, _TYPE=args.type, _HLA=args.HLA,
                    _p_HLA_INTEGRATED_POSITIONS=args.HLA_pos)
