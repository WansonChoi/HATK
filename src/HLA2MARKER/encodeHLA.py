import os, re
import pandas as pd
import argparse, textwrap
from collections import OrderedDict


def encodeHLA(_INPUT_PED, _OUTPUT, _hg = "19"):


    ########## <Core Variables> ##########

    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    # isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}

    # (2018. 9. 25.) Replaced by lift-over values.
    genepos_hg ={"18" : {"A": 30018226, "C": 31344505, "B": 31429628, "DRB1": 32654525, "DQA1": 32713161, "DQB1": 32735219, "DPA1": 33140324, "DPB1": 33151681},
                 "19" : {"A": 29910247, "C": 31236526, "B": 31321649, "DRB1": 32546547, "DQA1": 32605183, "DQB1": 32627241, "DPA1": 33032346, "DPB1": 33043703},
                 "38" : {"A": 29942470, "C": 31268749, "B": 31353872, "DRB1": 32578770, "DQA1": 32637406, "DQB1": 32659464, "DPA1": 33064569, "DPB1": 33075926}}

    # Module Name for stdout
    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.")


    # Set of allele names ocurring in each HLA columns of given ped file.
    ALLELE_TABLES = OrderedDict()
    ALLELE_TABLES_1field = OrderedDict()

    ALL_ALLELES = []
    dict_ALL_ALLELES = {}


    df_OUTPUT_map = pd.DataFrame()
    df_OUTPUT_ped = pd.DataFrame()




    ########## < Control Flags > ##########


    LOADING_PEDFILE = 1
    MAKING_ALLELE_TABLE = 1
    MAKING_OUTPUT_MAP = 1
    MAKING_OUTPUT_PED = 1
    ADDING_DUMMY_MARKER = 1





    if LOADING_PEDFILE:

        ########## < 1. Loading Input PED file > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Loading Input PED file.\n")

        INPUT_PED = pd.read_table(_INPUT_PED, sep='\t', header=None, dtype=str,
                                  names = ['Fam_ID', 'Sample_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phe'] + [''.join([HLA_names[i], '_', str(j)]) for i in range(0, len(HLA_names)) for j in range(1,3)],
                                  index_col=[0, 1, 2, 3, 4, 5]
                                  )

        print(INPUT_PED.head())



    if MAKING_ALLELE_TABLE:

        ########## < 2. Making Allele Table > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Making Allele Table.\n")

        """
        for i in range(0, len(INPUT_PED.index)):

            line = tuple(INPUT_PED.iloc[i, :])

            alleles["HLA_A_"+line[6]] = genepos["HLA_A"];       alleles["HLA_A_"+line[7]] = genepos["HLA_A"];
            alleles["HLA_A_"+line[6][0:2]] = genepos["HLA_A"];  alleles["HLA_A_"+line[7][0:2]] = genepos["HLA_A"];

            alleles["HLA_B_"+line[8]] = genepos["HLA_B"];       alleles["HLA_B_"+line[9]] = genepos["HLA_B"];
            alleles["HLA_B_"+line[8][0:2]] = genepos["HLA_B"];  alleles["HLA_B_"+line[9][0:2]] = genepos["HLA_B"];

            alleles["HLA_C_"+line[10]] = genepos["HLA_C"];       alleles["HLA_C_"+line[11]] = genepos["HLA_C"];
            alleles["HLA_C_"+line[10][0:2]] = genepos["HLA_C"];  alleles["HLA_C_"+line[11][0:2]] = genepos["HLA_C"];

            alleles["HLA_DPA1_"+line[12]] = genepos["HLA_DPA1"];       alleles["HLA_DPA1_"+line[13]] = genepos["HLA_DPA1"];
            alleles["HLA_DPA1_"+line[12][0:2]] = genepos["HLA_DPA1"];  alleles["HLA_DPA1_"+line[13][0:2]] = genepos["HLA_DPA1"];

            alleles["HLA_DPB1_"+line[14]] = genepos["HLA_DPB1"];       alleles["HLA_DPB1_"+line[15]] = genepos["HLA_DPB1"];
            alleles["HLA_DPB1_"+line[14][0:2]] = genepos["HLA_DPB1"];  alleles["HLA_DPB1_"+line[15][0:2]] = genepos["HLA_DPB1"];

            alleles["HLA_DQA1_"+line[16]] = genepos["HLA_DQA1"];       alleles["HLA_DQA1_"+line[17]] = genepos["HLA_DQA1"];
            alleles["HLA_DQA1_"+line[16][0:2]] = genepos["HLA_DQA1"];  alleles["HLA_DQA1_"+line[17][0:2]] = genepos["HLA_DQA1"];

            alleles["HLA_DQB1_"+line[18]] = genepos["HLA_DQB1"];       alleles["HLA_DQB1_"+line[19]] = genepos["HLA_DQB1"];
            alleles["HLA_DQB1_"+line[18][0:2]] = genepos["HLA_DQB1"];  alleles["HLA_DQB1_"+line[19][0:2]] = genepos["HLA_DQB1"];

            alleles["HLA_DRB1_"+line[20]] = genepos["HLA_DRB1"];       alleles["HLA_DRB1_"+line[21]] = genepos["HLA_DRB1"];
            alleles["HLA_DRB1_"+line[20][0:2]] = genepos["HLA_DRB1"];  alleles["HLA_DRB1_"+line[21][0:2]] = genepos["HLA_DRB1"];


        for k,v in alleles.items():
            print("key : {0} / value : {1}".format(k, v))
            
        """

        # In the past, `ALLELE_TABLES` was created based on dictionary data structure, but now it is created by "apply()" function with "set()" function.

        for i in range(0, len(HLA_names)):
            # for i in range(0, 1):

            temp = INPUT_PED.filter(regex=HLA_names[i]+'_\d', axis=1).apply(set, axis=0).apply(lambda x : x.difference({0, "0"}))

            # print(temp)

            # Column 1 and 2
            set_al1 = temp.iat[0] # ex)     A_1 := {A*32:01:01, A*23:01:01, A*26:01:01, A*02:01:0...
            set_al2 = temp.iat[1] # ex)     A_2 := {A*24:02:01:01, A*01:01:01:01, A*02:05:01, A*3...

            # sr_Unioned_Set = pd.Series(list(set_al1.union(set_al2))).sort_values() # sorting

            l_Unioned_Set = list(set_al1.union(set_al2))
            l_Unioned_Set.sort() # sorting

            # l_Unioned_Set = pd.Series(l_Unioned_Set)
            print(l_Unioned_Set)

            ALLELE_TABLES[HLA_names[i]] = l_Unioned_Set


            ##### Dealing with 1-field #####

            if len(l_Unioned_Set) > 0:

                ### The case where the union of set of each two column has at least 1 element.

                sr_temp_1field = pd.Series(l_Unioned_Set).apply(lambda x : re.match(pattern='\*'.join([HLA_names[i], '\d{2,3}']), string=x).group()).unique() # Going through "unique()" function.

                # print("\nsr_temp_1field\n")
                # print(sr_temp_1field)

                ALLELE_TABLES_1field[HLA_names[i]] = sr_temp_1field.tolist()

            else:
                ### 집합 원소의 개수가 0 일때,(아예 input_ped에서 allele_name이 '0'으로 주어져서 없는 경우)
                ALLELE_TABLES_1field[HLA_names[i]] = l_Unioned_Set


        # for k, v in ALLELE_TABLES.items():
        #
        #     print("\n===============\n")
        #     print("{0} : \n{1}".format(k, v))



    if MAKING_OUTPUT_MAP:

        ########## < 3. Making OUTPUT .map file > ##########

        print(std_MAIN_PROCESS_NAME + "[3] Making OUTPUT .map file.\n")

        """        
        to_df_OUTPUT_map = []


        for name in HLA_names:
            for k in sorted_keys:
                temp = k.split('_')
                al = temp[2]

                if (temp[1] == name) and (al != "NA") and (al != "") and (al != "0") and (al != "0 0"):
                    # to_df_OUTPUT_map += [ ["6", k, "0", genepos['_'.join([temp[0], temp[1]])]] ]
                    to_df_OUTPUT_map.extend([("6", k, "0", genepos['_'.join([temp[0], temp[1]])])])

        df_OUTPUT_map = pd.DataFrame(to_df_OUTPUT_map)
        df_OUTPUT_map.to_csv(_OUTPUT + '.HLA.map', sep='\t', header=False, index=False)
                
        """


        """ 
        
        (1) map_LABELS
        (2) map_CHR
        (3) map_GENETIC_DISTANCE
        (4) map_POS
        
        """

        ##### Making Label for *.map file. #####

        # ALL_ALLELES = pd.concat([ALLELE_TABLES[HLA_names[i]].append(ALLELE_TABLES_1field[HLA_names[i]]) for i in range(0, len(HLA_names))]).sort_values()

        ALL_ALLELES = ALLELE_TABLES[HLA_names[0]]
        ALL_ALLELES.extend(ALLELE_TABLES_1field[HLA_names[0]])
        # ALL_ALLELES = [ALLELE_TABLES[HLA_names[i]].append(ALLELE_TABLES_1field[HLA_names[i]]) for i in range(0, len(HLA_names))]

        for i in range(1, len(HLA_names)):
            ALL_ALLELES.extend(ALLELE_TABLES[HLA_names[i]])
            ALL_ALLELES.extend(ALLELE_TABLES_1field[HLA_names[i]])

        ALL_ALLELES.sort()

        print("\nALL_ALLELES\n")
        print(ALL_ALLELES)


        ### HLA_index(Mining HLA gene name)
        sr_HLA = pd.Series(ALL_ALLELES).apply(lambda x : re.search(pattern='\w+\*', string=x).group().rstrip('*'))
        # print(sr_HLA)

        ### map_LABELS & map_POS ###
        map_LABELS = pd.Series(['HLA_' + ALL_ALLELES[i] for i in range(0, len(ALL_ALLELES))])
        # print("\n`map_LABELS`\n")
        # print(map_LABELS)

        map_POS = [genepos_hg[_hg][sr_HLA.iat[i]] for i in range(0, len(sr_HLA))]
        # print(map_POS)


        dict_ALL_ALLELES = {HLA_names[i] : [ALL_ALLELES[j] for j in range(0, len(ALL_ALLELES)) if (HLA_names[i]+'*' in ALL_ALLELES[j])] for i in range(0, len(HLA_names))}

        # print("\nsegmented `ALL_ALLELES`\n")
        # for k,v in dict_ALL_ALLELES.items():
        #     print("\n============\n")
        #     print("{0} : \n{1}".format(k, v))


        ##### map_Label을 제외한 나머지 map파일 항목 만들기 #####

        map_CHR = ['6' for i in range(0, len(map_LABELS))]
        map_GENETIC_DISTANCE = ['0' for i in range(0, len(map_LABELS))]


        df_OUTPUT_map = pd.DataFrame.from_dict({"Chr" : map_CHR, "Name" : map_LABELS.tolist(), "GD" : map_GENETIC_DISTANCE, "POS" : map_POS}).loc[:, ["Chr", "Name", "GD", "POS"]]
        print(std_MAIN_PROCESS_NAME + "Output .map file.\n")
        print(df_OUTPUT_map.head(50))

        # df_OUTPUT_map.to_csv('.'.join([_OUTPUT,'HLA.map']), sep='\t', header=False, index=False)



    if MAKING_OUTPUT_PED:

        ########## < 4. Making OUTPUT.ped file > ##########

        print(std_MAIN_PROCESS_NAME + "[4] Making .ped file.\n")


        """
        
                to_df_OUTPUT_ped = []
        
                for i in range(0, len(INPUT_PED.index)):
        
                    line = tuple(INPUT_PED.iloc[i, :])
        
                    to_df_OUTPUT_ped.extend([
                        line[0:6] +
                        PrintGenotypes("A", line[6], line[7], sorted_keys) +
                        PrintGenotypes("C", line[10], line[11], sorted_keys) +
                        PrintGenotypes("B", line[8], line[9], sorted_keys) +
                        PrintGenotypes("DRB1", line[20], line[21], sorted_keys) +
                        PrintGenotypes("DQA1", line[16], line[17], sorted_keys) +
                        PrintGenotypes("DQB1", line[18], line[19], sorted_keys) +
                        PrintGenotypes("DPA1", line[12], line[13], sorted_keys) +
                        PrintGenotypes("DPB1", line[14], line[15], sorted_keys)
                    ])
                        
                df_OUTPUT_ped = pd.DataFrame(to_df_OUTPUT_ped)
                df_OUTPUT_ped.to_csv(_OUTPUT + '.HLA.ped', sep='\t', header=False, index=False)
                    
                """


        to_df_OUTPUT_ped = []

        # for i in range(0, 5):
        for i in range(0, INPUT_PED.shape[0]):

            # print("\n================\n")

            line_INPUT_PED = tuple(INPUT_PED.iloc[i, :])
            # print(line_INPUT_PED)


            t_line_OUTPUT_PED = [PrintGenotypes3(line_INPUT_PED[2*j], line_INPUT_PED[2*j+1], dict_ALL_ALLELES[HLA_names[j]]) for j in range(0, len(HLA_names))]
            # print(t_line_OUTPUT_PED)
            # print(pd.Series(dict_ALL_ALLELES["A"]))

            # Flattening
            line_OUTPUT_PED = [item for eachlist in t_line_OUTPUT_PED for item in eachlist]

            # print("\nFlattened t_line_OUTPUT_PED is \n")
            # print(line_OUTPUT_PED)

            to_df_OUTPUT_ped.append(line_OUTPUT_PED)


        df_OUTPUT_ped = pd.DataFrame(to_df_OUTPUT_ped)
        df_OUTPUT_ped.index = INPUT_PED.index

        print(df_OUTPUT_ped.head())


        # df_OUTPUT_ped.to_csv('.'.join([_OUTPUT, 'HLA.ped']), sep='\t', header=False, index=True)




    if ADDING_DUMMY_MARKER:

        ########## < 5. Adding dummy_marker to ped and map files > ##########

        """
        for Compatitability with Plink1.07 version.
        """

        print(std_MAIN_PROCESS_NAME + "[5] Adding dummy_marker to final outputs.\n")


        df_OUTPUT_ped, df_OUTPUT_map = addDummyMarker(df_OUTPUT_ped, df_OUTPUT_map)

        df_OUTPUT_map.to_csv('.'.join([_OUTPUT,'HLA.map']), sep='\t', header=False, index=False)
        df_OUTPUT_ped.to_csv('.'.join([_OUTPUT, 'HLA.ped']), sep='\t', header=False, index=True)




    return 0






def PrintGenotypes3(_allele1, _allele2, _seg_ALL_ALLELES):

    l_output = []

    # print("\nAlleles : {0} and {1}".format(_allele1, _allele2))
    # print("\ndict_ALL_ALLELES: \n{0}".format(_seg_ALL_ALLELES))

    if len(_seg_ALL_ALLELES) > 0:

        if ((_allele1 != "0") and (_allele2 != "0")):  # The condition for checking integer value 0 won't be included here because .ped file was read with "dtype=str" option.

            t_sr1 = pd.Series(_seg_ALL_ALLELES, index=pd.Index(_seg_ALL_ALLELES)).apply(lambda x: "p" if (x in _allele1) else "a")
            # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr1))
            t_sr2 = pd.Series(_seg_ALL_ALLELES, index=pd.Index(_seg_ALL_ALLELES)).apply(lambda x: "p" if (x in _allele2) else "a")
            # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr2))

            for i in range(0, len(t_sr1)):
                l_output.append(t_sr1.iat[i])
                l_output.append(t_sr2.iat[i])

            # print(l_output)
            return l_output


        else:

            # Comments... at most below part.
            # print("At least one of allele is 0")

            for i in range(0, len(_seg_ALL_ALLELES)):
                l_output.append("0")
                l_output.append("0")

            return l_output

    else:
        # In cases such as "DPA1" or "DPB1 where any alleles don't appear, Just skip.
        # Then just return NULL
        return l_output




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






if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        Sherman Jia, 2012
        encodeHLA.py

        This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
        The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
                                            2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

    ###########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")



    ##### <for Test> #####

    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./TEST_0228", "--hg", "19"])

    # args = parser.parse_args(["./HAPMAP_CEU_HLA.ped", "./TEST_0228", "--hg", "19"])
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./TEST_0305", "-hg", "19"])

    # (2018. 7. 16.)

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "19",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])



    ##### <for Publication> #####

    args = parser.parse_args()


    print(args)

    # Implementing Main Function
    encodeHLA(args.ped, args.o, args.hg)
