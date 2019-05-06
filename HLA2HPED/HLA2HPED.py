# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
import json

########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

### General
HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}




def HATK_HLA2HPED(_rhped, _out, _platform):

    if not bool(_platform):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--platform\" wans't given. Please check it again.\n")
        sys.exit()

    if not bool(_rhped):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"-rhped\" wasn't given. Please check it again.\n")
        sys.exit()
    else:
        if not isinstance(_rhped, list):
            print(std_ERROR_MAIN_PROCESS_NAME + "The parameter \"_rhped\" to make \"hped\" wans't given as list. Please check the argument \"-rhped\" again.\n")
            sys.exit()

    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + "Output prefix wasn't given properly. Please check the argument \"--out\" again.\n")



    return HLA2HPED(_rhped, _out, _platform)



def HLA2HPED(_rhped, _out, _platform):

    ########## < Core Variables > ##########

    ### General

    # RETURN value
    OUTPUT_RETURN = None



    ########## < Argument Checking. > ##########

    # Preparing intermediate paths.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)



    ########## < Main conversion. > ##########

    if _platform == "AXIOM":

        print(std_MAIN_PROCESS_NAME + "Converting output files from \"AXIOM\" to \".hped\".\n")

        OUTPUT_RETURN = _convert_AXIOM(_rhped, _out)


    elif _platform == "HIBAG":

        print(std_MAIN_PROCESS_NAME + "Converting output files of \"HIBAG\" to \".hped\".\n")

        OUTPUT_RETURN = _convert_HIBAG(_rhped, _out)


    elif _platform == "xHLA":

        print(std_MAIN_PROCESS_NAME + "Converting output files of \"xHLA\" to \".hped\".\n")

        OUTPUT_RETURN = _convert_xHLA(_rhped, _out)

    else:

        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong platform specification.({0}) Please check it again.\n".format(_platform))
        sys.exit()



    # print("\n`OUTPUT_RETURN` is {0}\n".format(OUTPUT_RETURN))



    return OUTPUT_RETURN





def _convert_AXIOM(_i_HLA, _out):

    if len(_i_HLA) != 8:
        print(std_ERROR_MAIN_PROCESS_NAME + "Platform \"{0}\" needs 8 input files. Please check it again.\n".format("AXIOM"))
        sys.exit()


    DICT_HLA_INPUTS = \
        {HLA_names[i]: pd.read_table(_i_HLA[i], header=None, sep='\t|\s+', engine='python', dtype=str,
                                     names=["Label", "No", "Allele", "v1", "v2"]) for i in range(0, len(HLA_names))}

    N_rows = []


    for i in range(0, len(HLA_names)):

        print("\nGiven input(HLA : {0})".format(HLA_names[i]))
        print(DICT_HLA_INPUTS[HLA_names[i]].head())

        if DICT_HLA_INPUTS[HLA_names[i]].shape[0] % 2 == 1:
            # When given input has odd number of entities.
            print(std_WARNING_MAIN_PROCESS_NAME + "Input file of HLA-{0} has odd number of rows."
                                                  "Appendig last entity one more time by force.\n")

        N_rows.append(DICT_HLA_INPUTS[HLA_names[i]].shape[0])


    print("\nThe numbers of each HLA input files are : {0}".format(N_rows))



    l_df_RETURN = []

    for i in range(0, len(HLA_names)):

        df_temp = DICT_HLA_INPUTS[HLA_names[i]]

        l_temp = []

        for j in range(0, df_temp.shape[0], 2):

            l_temp.append([df_temp.iat[j, 2], df_temp.iat[j+1, 2]])


        df_temp2 = pd.DataFrame(l_temp)
        l_df_RETURN.append(df_temp2)

    df_RETURN = pd.concat(l_df_RETURN, axis=1)

    ### Index
    l_idx_RETURN = []

    # (2) SampleID
    idx_temp = df_temp.iloc[:, 0].tolist()
    idx_SampleID = [idx_temp[i] for i in range(0, len(idx_temp), 2)]

    # (1) FamID
    l_idx_RETURN.append(["FamID_"+str(i) for i in range(0, len(idx_SampleID))])
    l_idx_RETURN.append(idx_SampleID)

    # (3), (4), (5) : P_ID, M_ID, Sex(Unknown default)
    idx_temp = ["0" for i in range(0, len(idx_SampleID))]
    l_idx_RETURN.append(idx_temp)
    l_idx_RETURN.append(idx_temp)
    l_idx_RETURN.append(idx_temp)

    # (6) Phe (Unknown default.)
    idx_temp = ["-9" for i in range(0, len(idx_SampleID))]
    l_idx_RETURN.append(idx_temp)


    idx_RETURN = pd.MultiIndex.from_arrays(l_idx_RETURN)
    df_RETURN.index = idx_RETURN



    print("\ndf_AXIOM is \n")
    print(df_RETURN.head())

    # Exporting(File Writing)
    df_RETURN.to_csv(_out+".hped", sep='\t', header=False, index=True)


    return _out + ".hped"




def _convert_HIBAG(_i_HLA, _out):


    if len(_i_HLA) != 8:
        print(std_ERROR_MAIN_PROCESS_NAME + "Platform \"{0}\" needs 8 input files. Please check it again.\n".format("HIBAG"))
        sys.exit()


    DICT_HLA_INPUTS = \
        {HLA_names[i]: pd.read_table(_i_HLA[i], header=None, sep='\t|\s+', engine='python', dtype=str,
                                     names=["Label", "Allele1", "Allele2", "v1"]) for i in range(0, len(HLA_names))}


    N_rows = []
    l_df_RETURN = []

    for i in range(0, len(HLA_names)):

        print("\nGiven input(HLA : {0})".format(HLA_names[i]))
        print(DICT_HLA_INPUTS[HLA_names[i]].head())

        N_rows.append(DICT_HLA_INPUTS[HLA_names[i]].shape[0])

        l_df_RETURN.append(DICT_HLA_INPUTS[HLA_names[i]].iloc[:, [1, 2]])


    print("\nThe numbers of each HLA input files are : {0}\n".format(N_rows))


    df_RETURN = pd.concat(l_df_RETURN, axis=1)


    ### Index
    l_idx_RETURN = []

    # (1) FamID
    l_idx_RETURN.append(["FamID_"+str(i) for i in range(0, N_rows[-1])])

    # (2) SampleID
    idx_SampleID = DICT_HLA_INPUTS["A"].iloc[:, 0].tolist()
    l_idx_RETURN.append(idx_SampleID)

    # (3), (4), (5) : P_ID, M_ID, Sex(Unknown default)
    idx_temp = ["0" for i in range(0, len(idx_SampleID))]
    l_idx_RETURN.append(idx_temp)
    l_idx_RETURN.append(idx_temp)
    l_idx_RETURN.append(idx_temp)

    # (6) Phe (Unknown default.)
    idx_temp = ["-9" for i in range(0, len(idx_SampleID))]
    l_idx_RETURN.append(idx_temp)

    idx_REUTURN = pd.MultiIndex.from_arrays(l_idx_RETURN)
    df_RETURN.index = idx_REUTURN



    print("\ndf_HIBAG is \n")
    print(df_RETURN.head())

    # Exporting(File Writing)
    df_RETURN.to_csv(_out+".hped", sep='\t', header=False, index=True)

    return _out+".hped"



def _convert_xHLA(_i_HLA, _out):

    # It Assumes that json files are more than 1.

    DICT_json = []

    for i in range(0, len(_i_HLA)):

        with open(_i_HLA[i]) as f:

            DICT_json.append(json.load(f))


    """
    The output file from xHLA has 6 keys.
    
    (1) "subject_id"
    (2) "creation_time"
    (3) "report_version"
    (4) "report_type"
    (5) "sample_id"
    (6) "hla"
    
    Here, It will be assuemd that "subject_id" is same to "sample_id". 
    """

    l_SampleID = []
    l_HLA = []

    for i in range(0, len(_i_HLA)):

        l_SampleID.append(DICT_json[i]["subject_id"])

        l_temp = DICT_json[i]["hla"]["alleles"]
        l_eachRows = []

        for j in range(0, len(HLA_names)):

            l_temp2 = list(filter(lambda x : re.match(''.join(["^", HLA_names[j], "\*"]), x), l_temp))
            l_temp2 = ["0", "0"] if not bool(l_temp2) else l_temp2
            l_eachRows.extend(l_temp2)

        l_HLA.append(l_eachRows)

    print(l_SampleID)
    print(l_HLA)


    ### DataFrame
    df_RETURN = pd.DataFrame(l_HLA)

    ### Index
    l_idx_RETURN = []

    # (1) FamID
    l_idx_RETURN.append(["FamID_"+str(i) for i in range(0, len(l_SampleID))])

    # (2) SampleID
    l_idx_RETURN.append(l_SampleID)

    # (3), (4), (5) : P_ID, M_ID, Sex(Unknown default)
    idx_temp = ["0" for i in range(0, len(l_SampleID))]
    l_idx_RETURN.append(idx_temp)
    l_idx_RETURN.append(idx_temp)
    l_idx_RETURN.append(idx_temp)

    # (6) Phe (Unknown default.)
    idx_temp = ["-9" for i in range(0, len(l_SampleID))]
    l_idx_RETURN.append(idx_temp)

    idx_REUTURN = pd.MultiIndex.from_arrays(l_idx_RETURN)
    df_RETURN.index = idx_REUTURN


    print("\ndf_xHLA is \n")
    print(df_RETURN.head())



    ### Exporting(File Writing)
    df_RETURN.to_csv(_out+".hped", sep='\t', header=False, index=True)

    return _out+".hped"






def _convert_HISAT(_i_HLA, _out):

    # It is going to be introduced soon.

    return _out+".hped"






if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        < HLA2HPED.py >

         This Converts the output results from other HLA imputation software or framework to *.hped 
        format.
        
         The number of raw input file(s) to make a hped file could be 1, 8 or More.
        
         In case of 8 raw input files are given by "-rhped" argument(ex. AXIOM, HIBAG), those file 
        must be given with the below order.
        (8 HLA genes - "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1").
        
        
         List of available HLA sofware is:
         
            (1) HLA*IMP(Axiom)
            (2) HIBAG
            (3) xHLA
            (4) HISAT
        
        Other HLA software can be added later.
        
        

    #################################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-rhped", help="\nInput Data file(Output result(s) from other HLA related software)\n\n", nargs='*')
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--platform", "-p", help="\nSoftware platform.\n\n", required=True,
                        choices=["AXIOM", "HIBAG", "xHLA", "HISAT"])





    ##### < for Test > #####

    # args = parser.parse_args(["-p", "AXIOM", "-rhped",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_A_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_B_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_C_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_DPA1_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_DPB1_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_DQA1_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_DQB1_Results.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/AxiomHLA_TestResult/AxiomHLA_4dig_DRB1_Results.txt",
    #                           "-o", "HLA2HPED_removethis"
    #                           ])

    # args = parser.parse_args(["-p", "HIBAG", "-rhped",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/HIBAG_TestResult_HLA-A.txt",
    #                           "-o", "HLA2HPED_removethis_HIBAG"
    #                           ])

    # args = parser.parse_args(["-p", "xHLA", "-rhped",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/test.json",
    #                           "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/HLA2HPED/test2.json",
    #                           "-o", "HLA2HPED_removethis_xHLA"
    #                           ])


    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    # Main function execution
    HLA2HPED(args.rhped, args.out, args.platform)