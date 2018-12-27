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



def encodeVariants(_p_ped, _p_map, _out, __use_pandas=False, __return_as_DataFrame=False):


    if __use_pandas:

        ########## < Control Flags > ##########

        _1_LOADING_MAPFILE = 1
        _2_ALLELE_OVERLAPPING = 1
        _3_MAKING_NEW_PEDFILE = 1
        _4_MAKING_NEW_MAPFILE = 1
        _5_ADDING_DUMMY_MARKER = 1


        if _1_LOADING_MAPFILE:

            ########## < 1. Loading MAP file > ##########

            print(std_MAIN_PROCESS_NAME + "[1] Loading MAP file.\n")

            # Processing map file
            df_mapfile = pd.read_table(_p_map, sep='\t', header=None, dtype=str)
            print(df_mapfile.head())



        if _2_ALLELE_OVERLAPPING:

            ########## < 2. Allele overlapping > ##########

            print(std_MAIN_PROCESS_NAME + "[2] Allele Overlapping.\n")

            df_pedfile = pd.read_table(_p_ped, sep='\t', header=None, dtype=str, index_col=[0, 1, 2, 3, 4, 5])
            print(df_pedfile.head())

            flattened = df_pedfile.apply(set, axis=0)
            print(flattened)

            # alleles = [tuple(flattened.iat[i].union(flattened.iat[i+1]).difference({0, "0"})) for i in range(0, flattened.shape[0], 2)]
            # print(alleles)

            # (2018. 9. 1.)

            alleles = []

            for i in range(0, flattened.shape[0], 2):

                l_temp = list(flattened.iat[i].union(flattened.iat[i+1]).difference({0, "0"}))
                l_temp.sort()

                alleles.append(l_temp)

            # (2018. 8. 28.)
            sr_alleles = pd.Series(alleles, index=pd.Index(df_mapfile.iloc[:, 1]))
            sr_alleles.to_csv(_out+".alleleset", sep='\t', header=False, index=True)

            # Now, `alleles`(list of list) will be used to generate new .ped file.
            print("Making Alleles list is done!")



        if _3_MAKING_NEW_PEDFILE:

            ########## < 3. Making new .ped file > ##########

            print(std_MAIN_PROCESS_NAME + "[3] Making new .ped file.\n")

            # each lines will be converted to lists,
            # and these lists will be gathered to make list of lists. Finally this list of lists will be DataFrame.

            to_DataFrame = []

            # Iterating over each lines of .ped files
            for N_pedline in range(0, df_pedfile.shape[0]):

                curr_line = list(df_pedfile.iloc[N_pedline, :])

                # ID = curr_line[0:6]
                Seq = []

                # for i in range(6, len(curr_line), 2):
                for i in range(0, df_pedfile.shape[1], 2):

                    # (SNP1, SNP2) : (0,1) -> (2,3) -> (4, 5) -> ...

                    SNP1 = curr_line[i]
                    SNP2 = curr_line[i+1]

                    idx_alleles = int(i/2)
                    # idx_Seq = i-6

                    # alleles[idx_alleles] := the number of sort of alleles on each positions.
                    if len(alleles[idx_alleles]) > 2:

                        for j in range(0, len(alleles[idx_alleles])):

                            if SNP1 == "0" or SNP2 == "0":
                                Seq.append("0"); Seq.append("0")

                            else:

                                if alleles[idx_alleles][j] == SNP1:
                                    Seq.append("p")
                                else:
                                    Seq.append("a")

                                if alleles[idx_alleles][j] == SNP2:
                                    Seq.append("p")
                                else:
                                    Seq.append("a")


                        if len(alleles[idx_alleles]) > 3:

                            j_end = 1 if len(alleles[idx_alleles]) == 4 else len(alleles[idx_alleles])

                            for j in range(0, j_end):

                                for k in range(j+1, len(alleles[idx_alleles])):

                                    if SNP1 == "0" or SNP2 == "0":
                                        Seq.append("0"); Seq.append("0")

                                    else:
                                        if alleles[idx_alleles][j] == SNP1 or alleles[idx_alleles][k] == SNP1:
                                            Seq.append("p")
                                        else:
                                            Seq.append("a")

                                        if alleles[idx_alleles][j] == SNP2 or alleles[idx_alleles][k] == SNP2:
                                            Seq.append("p")
                                        else:
                                            Seq.append("a")


                            if len(alleles[idx_alleles]) > 5:

                                j_end = 1 if len(alleles[idx_alleles]) == 6 else len(alleles[idx_alleles])

                                for j in range(0, j_end):
                                    for k in range(j+1, len(alleles[idx_alleles])):
                                        for l in range(k+1, len(alleles[idx_alleles])):

                                            if SNP1 == "0" or SNP2 == "0":
                                                Seq.append("0"); Seq.append("0")

                                            else:
                                                if alleles[idx_alleles][j] == SNP1 or alleles[idx_alleles][k] == SNP1 or alleles[idx_alleles][l] == SNP1:
                                                    Seq.append("p")
                                                else:
                                                    Seq.append("a")

                                                if alleles[idx_alleles][j] == SNP2 or alleles[idx_alleles][k] == SNP2 or alleles[idx_alleles][l] == SNP2:
                                                    Seq.append("p")
                                                else:
                                                    Seq.append("a")



                    else:
                        # Most of Cases have length less than equal 2, they will fall into this if-else block.
                        if SNP1 == "0" or SNP2 == "0":
                            Seq.append("0"); Seq.append("0")
                        else:
                            Seq.append(SNP1); Seq.append(SNP2)


                to_DataFrame.append(Seq)

                # End-point of one line of .ped file.

            df_output_pedfile = pd.DataFrame(to_DataFrame, index=df_pedfile.index)
            # df_output_pedfile.to_csv(_out + '.ped', sep='\t', header=False, index=True)

            print("Making .ped file is done")



        if _4_MAKING_NEW_MAPFILE:

            ########## < 4. Making new .map file > ##########

            print(std_MAIN_PROCESS_NAME + "[4] Making new .map file.\n")

            to_DataFrame = []

            # print(len(alleles))
            # print(len(df_mapfile.index))
            # len(alleles) == len(df_mapfile.index)

            for i in range(0, len(alleles)):

                # print(alleles[i])
                # print(list(df_mapfile.iloc[i, :]))

                curr_line = tuple(df_mapfile.iloc[i, :])
                # [0]: 6, [1]: 'AA_A_9_30126516', [2]: 0, [3]: 30126516

                if len(alleles[i]) > 2:

                    # multi_alleles_1,2,3 => These are all list of lists

                    multi_alleles_1 = [(curr_line[0], curr_line[1]+'_'+alleles[i][j], curr_line[2], curr_line[3]) for j in range(0, len(alleles[i]))]
                    # to_DataFrame += multi_alleles_1
                    to_DataFrame.extend(multi_alleles_1)

                    if len(alleles[i]) > 3:

                        j_end = 1 if len(alleles[i]) == 4 else len(alleles[i])

                        # for j in range(0, j_end):
                        #     for k in range(j+1, len(alleles[i])):
                        #         print([curr_line[0], curr_line[1]+ '_' + alleles[i][j]+alleles[i][k], curr_line[2], curr_line[3]])

                        multi_alleles_2 = [(curr_line[0], curr_line[1]+ '_' + alleles[i][j]+alleles[i][k], curr_line[2], curr_line[3]) for j in range(0, j_end) for k in range(j+1, len(alleles[i]))]
                        # to_DataFrame += multi_alleles_2
                        to_DataFrame.extend(multi_alleles_2)

                        if len(alleles[i]) > 5:

                            j_end = 1 if len(alleles[i]) == 6 else len(alleles[i])

                            # for j in range(0, j_end):
                            #     for k in range(j+1, len(alleles[i])):
                            #         for l in range(k+1, len(alleles[i])):

                            multi_alleles_3 = [(curr_line[0], curr_line[1]+ '_' + alleles[i][j]+alleles[i][k]+alleles[i][l], curr_line[2], curr_line[3]) for j in range(0, j_end) for k in range(j+1, len(alleles[i])) for l in range(k+1, len(alleles[i]))]
                            # to_DataFrame += multi_alleles_3
                            to_DataFrame.extend(multi_alleles_3)

                # Job to divide multi-allele is done.


                else:
                    to_DataFrame.extend([curr_line])

            df_output_mapfile = pd.DataFrame(to_DataFrame)
            # df_output_mapfile.to_csv(_out + '.map', sep='\t', header=False, index=False)

            # print(df_output_mapfile.head())

            print("Making new .map file is done!")




        if _5_ADDING_DUMMY_MARKER:

            ########## < 5. Adding dummy_marker to ped and map files > ##########

            """
            for Compatitability with Plink1.07 version.
            """

            print(std_MAIN_PROCESS_NAME + "[5] Adding dummy_marker to final outputs.\n")


            df_OUTPUT_ped, df_OUTPUT_map = addDummyMarker(df_output_pedfile, df_output_mapfile)


            df_OUTPUT_ped.to_csv(_out + '.ped', sep='\t', header=False, index=True)
            df_OUTPUT_map.to_csv(_out + '.map', sep='\t', header=False, index=False)



    else:

        ### Processing line by line to save memory

        process = psutil.Process(os.getpid())

        mem_start = process.memory_info().rss / 1024**2
        t1 = time.clock()


        HARDCODING = 1

        if not HARDCODING:

            # Acquiring column number
            n_loci = 0
            n_col = 0
            n_row = 0

            with open(_p_ped, 'r') as f:

                for line in f:

                    t_line = re.split(r'\s+', line.rstrip('\n'))
                    # print(t_line[6:26])

                    if n_row == 0:

                        ped_info = t_line[:6]
                        genomic_info = t_line[6:]

                        n_loci = len(genomic_info)
                        if n_loci % 2 == 1: print(std_ERROR_MAIN_PROCESS_NAME + "The number of alleles is not even.\n"); sys.exit()
                        n_col = len(t_line)

                        # list containing factors which appear in each locus.
                        l_factors = [[genomic_info[i]] for i in range(0, int(n_loci/2))] # Initialization
                        # print(l_factors)

                    elif n_row > 0:

                        for i in range(0, int(n_loci/2)):

                            idx1 = 2*i + 6 # index for `t_line`
                            idx2 = idx1 + 1

                            ##### Allele overlapping
                            if t_line[idx1] != "0" and (t_line[idx1] not in l_factors[i]):
                                l_factors[i].append(t_line[idx1])
                            if t_line[idx2] != "0" and (t_line[idx2] not in l_factors[i]):
                                l_factors[i].append(t_line[idx2])


                    n_row += 1

                    # if n_row > 5:
                    #     break


            # print(l_factors[:10])
            print("Memory size of overlapped factor info : {}(Mb)".format(sys.getsizeof(l_factors)/1024**2))

            # Sorting elements of each lists.
            for i in range(0, len(l_factors)):
                l_factors[i].sort()

            # Writing file.
            g1 = ('\t'.join(l_factors[i])+"\n" for i in range(0, len(l_factors)))

            with open(_out+".temp.alleleset", 'w') as f_alleleset:
                f_alleleset.writelines(g1)

        else:

            ### (Warning) Temporary Hardcoding

            l_factors = []

            with open(_out+".temp.alleleset", 'r') as f_alleleset:
                for line in f_alleleset:
                    t_list = re.split(r'\s+', line.rstrip('\n'))
                    t_list.sort()
                    l_factors.append(t_list)


            print(l_factors[:20])

            ### --- `l_factors` done.




        ##### Making New ped file #####

        with open(_out + ".ped", 'w') as f_output:
            f_output.writelines(stackLines(_p_ped, l_factors))




        mem_end = process.memory_info().rss / 1024**2
        t2 = time.clock()

        print("\n\nMemory before : {}(Mb)\nMemory after : {}(Mb)".format(mem_start, mem_end))
        print("Difference(Memory usage) : {}(Mb)".format(mem_end-mem_start))
        print("Time : {}(sec)".format(t2 - t1))







def addDummyMarker(_df_ped, _df_map):

    # adding dummy marker to map file.
    sr_dummy_marker = pd.DataFrame([["6", "dummy_marker", "0", "33999999"]], index=pd.Index([_df_map.shape[0]]), columns=_df_map.columns)
    df_RETURN_map = pd.concat([_df_map, sr_dummy_marker], axis=0)


    print("\nMap file with dummy_marker : \n")
    print(df_RETURN_map.tail())


    # adding two dummy columns to ped file.
    sr_dummy_columns = pd.Series([("d" if bool(i % 2) else "D") for i in range(0, _df_ped.shape[0])], index=_df_ped.index)
    df_RETURN_ped = pd.concat([_df_ped, sr_dummy_columns, sr_dummy_columns], axis=1)

    print("\nPed file with dummy_marker : \n")
    print(df_RETURN_ped.head())


    return [df_RETURN_ped, df_RETURN_map]




def divideToBinaryMarkers(_SNP1, _SNP2, _factors):

    Seq = []

    if len(_factors) > 2:

        for j in range(0, len(_factors)):

            if _SNP1 == "0" or _SNP2 == "0":
                Seq.append("0"); Seq.append("0")

            else:

                if _factors[j] == _SNP1:
                    Seq.append("p")
                else:
                    Seq.append("a")

                if _factors[j] == _SNP2:
                    Seq.append("p")
                else:
                    Seq.append("a")

        if len(_factors) > 3:

            j_end = 1 if len(_factors) == 4 else len(_factors)

            for j in range(0, j_end):

                for k in range(j + 1, len(_factors)):

                    if _SNP1 == "0" or _SNP2 == "0":
                        Seq.append("0"); Seq.append("0")

                    else:
                        if _factors[j] == _SNP1 or _factors[k] == _SNP1:
                            Seq.append("p")
                        else:
                            Seq.append("a")

                        if _factors[j] == _SNP2 or _factors[k] == _SNP2:
                            Seq.append("p")
                        else:
                            Seq.append("a")

            if len(_factors) > 5:

                j_end = 1 if len(_factors) == 6 else len(_factors)

                for j in range(0, j_end):
                    for k in range(j + 1, len(_factors)):
                        for l in range(k + 1, len(_factors)):

                            if _SNP1 == "0" or _SNP2 == "0":
                                Seq.append("0"); Seq.append("0")

                            else:
                                if _factors[j] == _SNP1 or _factors[k] == _SNP1 or _factors[l] == _SNP1:
                                    Seq.append("p")
                                else:
                                    Seq.append("a")

                                if _factors[j] == _SNP2 or _factors[k] == _SNP2 or _factors[l] == _SNP2:
                                    Seq.append("p")
                                else:
                                    Seq.append("a")



    else:
        # Most of Cases have length less than equal 2, they will fall into this if-else block.
        if _SNP1 == "0" or _SNP2 == "0":
            Seq.append("0"); Seq.append("0")
        else:
            Seq.append(_SNP1); Seq.append(_SNP2)

    return '\t'.join(Seq)



def stackLines(_p_ped, _l_factors):

    count = 0

    with open(_p_ped, 'r') as f:
        for line in f:
            t_line = re.split(r'\s+', line.rstrip('\n'))

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                divideToBinaryMarkers(t_line[2 * i + 6], t_line[2 * i + 7], _l_factors[i]) for i in range(0, len(_l_factors))
            ])

            # add Dummy Markers.
            dummy_markers = '\t'.join(['d', 'D'] if bool(count % 2) else ['D', 'd'])


            yield '\t'.join([__ped_info__, __genomic_info__, dummy_markers]) + "\n"





if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        Original Author : Sherman Jia, 2012
     
        encodeVariants.py
     
        - This script encodes multi-allelic PLINK markers (amino acids and SNPs) into bi-allelic 
            markers 
        - Input files include a normal PLINK .ped and .map file (where the .ped file contains 
            multi-allelic positions).

    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-map", help="\nMap file for given \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)


    ##### <for Test> #####

    # (2018. 2. 26)
    # args = parser.parse_args(["TEST_v2.AA.ped", "TEST_v2.AA.map", "TEST_v2"])

    # (2018. 7. 13.)

    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.enCODED"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/CHECK_HLAtoSeuqences.txt.AA.ped",
    #                           "-map", "/Users/wansun/Data/HATK/data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.map",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/CHECK_HLAtoSeuqences.txt"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/PED_onlyHLA_C.txt.AA.ped",
    #                           "-map", "/Users/wansun/Data/HATK/data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.map",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/PED_onlyHLA_C.txt"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK.AA.map",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK.AA.CODED"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.enCODED"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK.SNPS.map",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK.SNPS.CODED"])

    # # (2018. 12. 27.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/20181221/wtccc_filtered_58C_NBS_RA.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/20181221/wtccc_filtered_58C_NBS_RA.SNPS.map",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/20181221/wtccc_filtered_58C_NBS_RA.SNPS.CODED"])



    ##### <for Publication> #####

    args = parser.parse_args()



    print(args)


    # Implementing Main Function.
    encodeVariants(args.ped, args.map, args.o)