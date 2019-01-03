# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd
import argparse, textwrap
import psutil
import time


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def HLAtoSequences(_chped, _dictionary, _type, _out, __use_pandas=False, __return_as_DataFrame=False):


    ########## < Core Variables > ##########

    # HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    # isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}

    HLA_DICTIONARY = pd.DataFrame()
    HLA_DICTIONARY_byHLA = {}
    HLA_SEQ_LENGTH = {}


    ########## < Argument checking > ##########

    # (1) ped file existence
    if not os.path.exists(_chped):
        print(_chped)
        print(std_ERROR_MAIN_PROCESS_NAME + "Given ped file dosen't exist. Please check it again.\n")
        sys.exit()

    # (2) HLA DICTIONARY file
    if not os.path.exists(_dictionary):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given dictionary file dosen't exist. Please check it again.\n")
        sys.exit()

    # (3) Chekcing `_type`
    if not (_type == "AA" or _type == "SNPS"):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given value for argument `_type` has wrong value. Please check it again.\n")
        sys.exit()




    """
    (2018. 12. 27.)
    For the memory issue, previous codes and logic bsaed on pandas will be used only when `__return_as_DataFrame` is True.
    In other cases, generting sequence informtion will be done with just plain for loops and file writing(I/O).
    """

    if __use_pandas:


        ########## < Control Flags > ##########

        LOADING_DICTIONARY = 1
        LOADING_CHPED = 1
        BRINGING_SEQUENCES = 1
        EXPORTING_OUTPUT_PED = 1



        if LOADING_DICTIONARY:

            ########## <1. Dictionary Preparation> ##########

            # print("\n[1] Loading Dictionary Data.\n")

            if os.path.isfile(_dictionary):
                HLA_DICTIONARY = pd.read_table(_dictionary, sep='\t', header=None, names=["Alleles", "Seqs"], index_col=0)
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Dictionary file doesn't exit!\n")
                sys.exit()

            # print(HLA_DICTIONARY.head())


            ### Dividing `HLA_DICTIONARY` in advance.

            # Dividing `HLA_DICTIONARY` by HLA. (`HLA_DICTIONARY_byHLA`)
            HLA_DICTIONARY_byHLA = {HLA_names[i]: HLA_DICTIONARY.filter(regex= ''.join(["^", HLA_names[i], "\*"]), axis=0) for i in range(0, len(HLA_names))}

            # for i in range(0, len(HLA_names)):
            #     print("\nSequence Information of %s" % HLA_names[i])
            #     print(HLA_DICTIONARY_byHLA[HLA_names[i]].head())

            HLA_SEQ_LENGTH = {HLA_names[i] : len(HLA_DICTIONARY_byHLA[HLA_names[i]].iat[0, 0]) for i in range(0, len(HLA_names))}

            # for i in range(0, len(HLA_names)):
            #     print("\nSequence Length of %s" % HLA_names[i])
            #     print(HLA_SEQ_LENGTH[HLA_names[i]])



        if LOADING_CHPED:

            ########## <2. Loading Cleaned HPED(*.chped) file> ##########

            # print("\n[2] Loading Input CHPED file.")

            INPUT_CHPED = pd.read_table(_chped, sep='\t', header=None, index_col=[0, 1, 2, 3, 4, 5], dtype=str)
            INPUT_CHPED.columns = pd.Index([name + '_' + str(j + 1) for name in HLA_names for j in range(0, 2)])

            # print(INPUT_CHPED.head())



        if BRINGING_SEQUENCES:

            ########## <3. Bringing Sequences> ##########

            print("\n[3]: Transforming Allele names to Sequences.")

            l_FOR_OUTPUT = []
            l_FOR_OUTPUT_test = []

            for i in range(0, len(HLA_names)):

                curr_dict = HLA_DICTIONARY_byHLA[HLA_names[i]]

                df_temp = INPUT_CHPED.filter(regex='_'.join([HLA_names[i], '\d{1}']), axis=1)
                # print(df_temp.head())

                df_temp = df_temp.applymap(lambda x : BringSequence(x, curr_dict) if x != "0" else x)
                # print(df_temp.head())

                l_FOR_OUTPUT_test.append(df_temp)

                # print("\n===============\n")

                # Now, we need to iterate on the number of rows of this DataFrame

                COLUMNS_BY_HLA = []
                # COLUMNS_BY_HLA_test = []

                for j in range(0, len(df_temp)):

                    if df_temp.iat[j, 0] != '0' and df_temp.iat[j, 1] != '0':

                        # (Case1) Most normal case - wehn allele_name is found as pair.
                        # ex) allele1 : A*25:01:01  /  allele2 : A*03:01:01:01

                        # seq1 = df_temp.iat[j, 0] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 0][::-1]
                        # seq2 = df_temp.iat[j, 1] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 1][::-1]

                        # (2018. 3. 9) 다시 여기서 reverse안시키는 걸로 바꿈
                        seq1 = df_temp.iat[j, 0]
                        seq2 = df_temp.iat[j, 1]

                        PAIRED = [value for item in zip(seq1, seq2) for value in item]

                    else:

                        # (Case2) when not found as a pair of alleles, but as a only single allele, => 0 is given
                        # (0, 0) will be compensated as length of `HLA_seq`, Due to here, I need to prepared `len(HLA_seq)` beforehand.

                        seq1 = ['0' for z in range(0, HLA_SEQ_LENGTH[HLA_names[i]])]

                        PAIRED = [value for item in zip(seq1, seq1) for value in item]

                    COLUMNS_BY_HLA.append(PAIRED)
                    # COLUMNS_BY_HLA_test.append(''.join(PAIRED))


                l_FOR_OUTPUT.append(pd.DataFrame(COLUMNS_BY_HLA))
                # l_FOR_OUTPUT_test.append(pd.Series(COLUMNS_BY_HLA_test))

                # End of interation. Don't get confused.



        if EXPORTING_OUTPUT_PED:

            ########## <4. Exporting OUTPUT PED file> ##########

            print("\n[4]: Exporting OUTPUT PED file.")

            df_OUTPUT = pd.concat(l_FOR_OUTPUT, axis=1)
            df_OUTPUT.index = INPUT_CHPED.index
            df_OUTPUT.columns = pd.Index(range(0, df_OUTPUT.shape[1]))

            print(df_OUTPUT.head())


            if __return_as_DataFrame:
                return df_OUTPUT
            else:
                ### Final Output ped file.
                if _type == 'AA':
                    df_OUTPUT.to_csv(_out + '.AA.ped', sep='\t', header=False, index=True)
                    return _out + ".AA.ped"
                elif _type == 'SNPS':
                    df_OUTPUT.to_csv(_out + '.SNPS.ped', sep='\t', header=False, index=True)
                    return _out + ".SNPS.ped"


    else:

        ### Generating sequence information based on plain for loops and File writing(File I/O).
        ### When memory limit is tight.


        # process = psutil.Process(os.getpid())
        #
        # mem_start = process.memory_info().rss / 1024**2
        # t1 = time.clock()


        ##### < [1] Loading HLA Dictionary > #####

        HLA_DICTIONARY_byHLA = {HLA_names[i] : {} for i in range(0, len(HLA_names))} # Initialization.
        HLA_SEQ_LENGTH = {HLA_names[i]: -1 for i in range(0, len(HLA_names))}

        with open(_dictionary, "r") as f_dictionary:

            for line in f_dictionary:

                t_line = re.split(r'\s+', line.rstrip('\n'))
                # ex) ['B*58:01:01', 'CTAGTCCTGCTTCAGGGTCCGGGGCCCG...']
                # t_line[0] := allele name, t_line[1] := the sequence

                for i in range(0, len(HLA_names)):

                    if re.match(r'{}\*'.format(HLA_names[i]), t_line[0]):

                        HLA_DICTIONARY_byHLA[HLA_names[i]][t_line[0]] = t_line[1]

                        if HLA_SEQ_LENGTH[HLA_names[i]] == -1:
                            HLA_SEQ_LENGTH[HLA_names[i]] = len(t_line[1])

                        break




        # # Result check
        # for i in range(0, len(HLA_names)):
        #     print("\n{} :\n".format(HLA_names[i]))
        #     for k, v in HLA_DICTIONARY_byHLA[HLA_names[i]].items():
        #         print("{} : {}".format(k, v))
        #
        # for k, v in HLA_SEQ_LENGTH.items():
        #     print("The length of HLA-{} : {}".format(k, v))
        #
        #
        # print("Mem of HLA dictionary : {}(Mb)".format(sys.getsizeof(HLA_DICTIONARY_byHLA)/1024**2))




        ##### < [2] Transforming each HLA alleles to corresponding sequences > #####

        with open(_out + ".{}.ped".format(_type), 'w') as f_output:
            f_output.writelines(GenerateLines(_chped, HLA_DICTIONARY_byHLA, HLA_SEQ_LENGTH))


        # mem_end = process.memory_info().rss / 1024**2
        # t2 = time.clock()
        #
        # print("\n\nMemory before : {}(Mb)\nMemory after : {}(Mb)".format(mem_start, mem_end))
        # print("Difference(Memory usage) : {}(Mb)".format(mem_end-mem_start))
        # print("Time : {}(sec)".format(t2 - t1))


        return _out + ".{}.ped".format(_type)





def BringSequence(_single_allele, _dict):

    try:
        Seq = _dict.loc[_single_allele, "Seqs"]
    except KeyError:
        Seq = "0"

    return Seq


def BringSequence2(_HLA_allele1, _HLA_allele2, _dict, _length):

    # # checking `_dict`
    # count = 0
    #
    # for k,v in _dict.items():
    #     print("{} : {}".format(k, v))
    #
    #     count += 1
    #     if count > 10:
    #         break

    # print("al1 : {}\nal2 : {}".format(_HLA_allele1, _HLA_allele2))

    # Finding the corresponding sequence of `_HLA_allele1`.
    try:
        Seq1 = _dict[_HLA_allele1]
        # Seq1 = _dict.get(_HLA_allele1)
    except KeyError:
        Seq1 = "0"

    # Same job for `_HLA_allele2`.
    try:
        Seq2 = _dict[_HLA_allele2]
        # Seq2 = _dict.get(_HLA_allele2)
    except KeyError:
        Seq2 = "0"


    # print("Corresponding Seqs : \n{}\n{}".format(Seq1, Seq2))


    if Seq1 != "0" and Seq2 != "0":
        return '\t'.join([value for item in zip(Seq1, Seq2) for value in item])
    else:
        return '\t'.join(["0" for z in range(0, 2*_length)])




def GenerateLines(_chped, _dict, _dict_length):

    with open(_chped, "r") as f:

        # mem_p1 = process.memory_info().rss / 1024**2


        for line in f:
            t_line = re.split(r'\s+', line.rstrip('\n'))
            # print(t_line)

            """
            [0,1,2,3,4,5] := ped file information
            [6,7] := HLA-A,
            [8,9] := HLA-B,
            ...,
            [20, 21] := HLA-DRB1
            """

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_part__ = '\t'.join([BringSequence2(t_line[(2 * i + 6)], t_line[(2 * i + 6 + 1)], _dict[HLA_names[i]], _dict_length[HLA_names[i]])
                                          for i in range(0, len(HLA_names))])

            # mem_p2 = process.memory_info().rss / 1024 ** 2
            # print("{}(Mb)".format(mem_p2 - mem_p1))

            yield '\t'.join([__ped_info__, __genomic_part__]) + "\n"




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

    parser.add_argument("-chped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-dict", help="\nHLA dictonary file name(ex. 'HLA_DICTIONARY_AA.txt')\n\n", required=True)
    parser.add_argument("-type", help="\nAA(for Amino Acid) or SNP(for SNPs)\n\n", choices=["AA", "SNPS"], required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)




    ##### <for Test> #####

    # # (2018. 8. 27.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Data/HATK/data/HLA_Analysis/AssociationTest/CancerResearch_Example/data_Rev_merged.4field.ped",
    #                           "-dict", "/Users/wansun/Data/HATK/data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/CHECK_HLAtoSeuqences.txt"])

    # # (2018. 8. 28.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Data/HATK/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/IMGT370_fixed/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Data/HATK/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/IMGT370_fixed/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK"])

    # args = parser.parse_args(["-chped", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/wtccc_filtered_58C_NBS_RA_T1D.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/20181221/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/20181221/wtccc_filtered_58C_NBS_RA.NewHLAtoSequences"])


    ##### <for Publication> #####

    args = parser.parse_args()


    print(args)


    # Implementing Main Function
    HLAtoSequences(args.chped, args.dict, args.type, args.o)


