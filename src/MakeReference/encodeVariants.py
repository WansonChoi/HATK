# -*- coding: utf-8 -*-

# (2017/11/25) recoded by Wanson Choi

###########################################################################################
#
# Created by Sherman Jia, 2012
# encodeVariants.pl
#
# This script encodes multi-allelic PLINK markers (amino acids and SNPs) into bi-allelic markers
# Input files include a normal PLINK .ped and .map file (where the .ped file contains
# multi-allelic positions).
#
###########################################################################################

import os
import pandas as pd
import argparse, textwrap

def encodeVariants(_p_ped, _p_map, _out):


    # Module Name for stdout
    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.")


    ########## < Control Flags > ##########

    _1_LOADING_MAPFILE = 1
    _2_ALLELE_OVERLAPPING = 1
    _3_MAKING_NEW_PEDFILE = 1
    _4_MAKING_NEW_MAPFILE = 1


    if _1_LOADING_MAPFILE:

        ########## < 1. Loading MAP file > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Loading MAP file.\n")

        # Processing map file
        df_mapfile = pd.read_table(_p_map, sep='\t', header=None, dtype=str)
        print(df_mapfile.head())



    if _2_ALLELE_OVERLAPPING:

        ########## < 2. Allele overlapping > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Allele Overlapping.\n")

        df_pedfile = pd.read_table(_p_ped, sep='\t', header=None, dtype=str, index_col=[0,1,2,3,4,5])
        print(df_pedfile.head())

        flattened = df_pedfile.apply(set, axis=0)
        print(flattened)

        alleles = [tuple(flattened.iat[i].union(flattened.iat[i+1]).difference({0, "0"})) for i in range(0, flattened.shape[0], 2)]
        print(alleles)

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
                                Seq.append("P")
                            else:
                                Seq.append("A")

                            if alleles[idx_alleles][j] == SNP2:
                                Seq.append("P")
                            else:
                                Seq.append("A")


                    if len(alleles[idx_alleles]) > 3:

                        j_end = 1 if len(alleles[idx_alleles]) == 4 else len(alleles[idx_alleles])

                        for j in range(0, j_end):

                            for k in range(j+1, len(alleles[idx_alleles])):

                                if SNP1 == "0" or SNP2 == "0":
                                    Seq.append("0"); Seq.append("0")

                                else:
                                    if alleles[idx_alleles][j] == SNP1 or alleles[idx_alleles][k] == SNP1:
                                        Seq.append("P")
                                    else:
                                        Seq.append("A")

                                    if alleles[idx_alleles][j] == SNP2 or alleles[idx_alleles][k] == SNP2:
                                        Seq.append("P")
                                    else:
                                        Seq.append("A")


                        if len(alleles[idx_alleles]) > 5:

                            j_end = 1 if len(alleles[idx_alleles]) == 6 else len(alleles[idx_alleles])

                            for j in range(0, j_end):
                                for k in range(j+1, len(alleles[idx_alleles])):
                                    for l in range(k+1, len(alleles[idx_alleles])):

                                        if SNP1 == "0" or SNP2 == "0":
                                            Seq.append("0"); Seq.append("0")

                                        else:
                                            if alleles[idx_alleles][j] == SNP1 or alleles[idx_alleles][k] == SNP1 or alleles[idx_alleles][l] == SNP1:
                                                Seq.append("P")
                                            else:
                                                Seq.append("A")

                                            if alleles[idx_alleles][j] == SNP2 or alleles[idx_alleles][k] == SNP2 or alleles[idx_alleles][l] == SNP2:
                                                Seq.append("P")
                                            else:
                                                Seq.append("A")



                else:
                    # Most of Cases have length less than equal 2, they will fall into this if-else block.
                    if SNP1 == "0" or SNP2 == "0":
                        Seq.append("0"); Seq.append("0")
                    else:
                        Seq.append(SNP1); Seq.append(SNP2)


            to_DataFrame.append(Seq)

            # End-point of one line of .ped file.

        df_output_pedfile = pd.DataFrame(to_DataFrame, index=df_pedfile.index)
        df_output_pedfile.to_csv(_out + '.ped', sep='\t', header=False, index=True)

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
        df_output_mapfile.to_csv(_out + '.map', sep='\t', header=False, index=False)

        # pd.Series(alleles).to_csv('alleles.txt', sep='\t', header=False)

        # print(df_output_mapfile.head())

        print("Making new .map file is done!")




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


    ##### <for Publication> #####

    args = parser.parse_args()



    print(args)


    # Implementing Main Function.
    encodeVariants(args.ped, args.map, args.o)