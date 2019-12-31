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



class HATK_HLA2HPED(object):

    def __init__(self, _rhped, _out, _platform):


        if not _rhped:
            print(std_ERROR_MAIN_PROCESS_NAME + "Raw hped file(s) wasn't(weren't) given. Please check '--rhped' argument again.")
            sys.exit()

        if not _platform:
            print(std_ERROR_MAIN_PROCESS_NAME + "Platform(ex. AXIOM, HIBAG, ...) infor wasn't given. Please check '--platform' argument again.")
            sys.exit()


        self.result = HLA2HPED(_rhped, _out, _platform)


    def getResult(self):
        return self.result



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

        # print(std_MAIN_PROCESS_NAME + "Converting output files from \"AXIOM\" to \".hped\".\n")

        OUTPUT_RETURN = _convert_AXIOM(_rhped, _out)


    elif _platform == "HIBAG":

        # print(std_MAIN_PROCESS_NAME + "Converting output files of \"HIBAG\" to \".hped\".\n")

        OUTPUT_RETURN = _convert_HIBAG(_rhped, _out)


    elif _platform == "xHLA":

        # print(std_MAIN_PROCESS_NAME + "Converting output files of \"xHLA\" to \".hped\".\n")

        OUTPUT_RETURN = _convert_xHLA(_rhped, _out)

    else:

        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong platform specification.({0}) Please check it again.\n".format(_platform))
        sys.exit()



    # print("\n`OUTPUT_RETURN` is {0}\n".format(OUTPUT_RETURN))



    return OUTPUT_RETURN





def _convert_AXIOM(_rhped, _out):

    if len(_rhped) != 8:
        print(std_ERROR_MAIN_PROCESS_NAME + "Platform \"{0}\" needs 8 input files. Please check it again.\n".format("AXIOM"))
        sys.exit()


    ### Exception handling 1 - When every rhped is 'NA'.

    f_isNA = list(map(lambda x : x == 'NA', _rhped))

    if all(f_isNA):
        print(std_ERROR_MAIN_PROCESS_NAME + "No any AXIOM output file has been given. Please check the '--rhped' argument again.")
        sys.exit()

    first_appear = -1

    for i in range(len(HLA_names)):
        if not f_isNA[i]:
            first_appear = i
            break


    DICT_rhped = {HLA_names[i]: pd.read_csv(_rhped[i], header=None, sep='\s+', dtype=str, names=["IID", "idx", "4digit", "p1", "p2"], usecols=['IID', 'idx', '4digit']) if not f_isNA[i] else pd.DataFrame([]) for i in range(0, len(HLA_names))}
    # print(DICT_rhped['B'])


    ### Restructuring rhped
    for i in range(len(HLA_names)):
        if not f_isNA[i]:
            DICT_rhped[HLA_names[i]] = DICT_rhped[HLA_names[i]].set_index(['IID', 'idx']).unstack(['idx'])
            # print(DICT_rhped[HLA_names[i]])


    ### Exception handling 2 - When there is any rhped which has different number of rows.
    L = [DICT_rhped[HLA_names[i]].shape[0] for i in range(len(HLA_names)) if not f_isNA[i]]
    if len(set(L)) > 1:
        print(std_ERROR_MAIN_PROCESS_NAME + "There is an HIBAG output file which has different number of rows."
                                            "Please check the HIBAG output files given to the '--rhped' argument.")
        sys.exit()
    else:
        L = L.pop()



    ### The Main step to generate hped file.

    l_HLA = []

    for i in range(len(HLA_names)):

        if not f_isNA[i]:
            df_temp = DICT_rhped[HLA_names[i]].reset_index(drop=True)
        else:
            df_temp = pd.DataFrame([['0', '0'] for z in range(L)])

        df_temp.columns = [HLA_names[i]+'_1', HLA_names[i]+'_2']
        # print(df_temp)
        l_HLA.append(df_temp)


    df_Right = pd.concat(l_HLA, axis=1)
    # print(df_Right)



    ### Meta information(Left 6 columns)

    #
    sr_PID = pd.Series(['0' for z in range(L)])
    sr_MID = pd.Series(['0' for z in range(L)])
    sr_Sex = pd.Series(['0' for z in range(L)])
    sr_Phe = pd.Series(['-9' for z in range(L)])

    df_Left = pd.concat([DICT_rhped[HLA_names[first_appear]].index.to_frame(index=False),
                         DICT_rhped[HLA_names[first_appear]].index.to_frame(index=False),
                         sr_PID, sr_MID, sr_Sex, sr_Phe], axis=1)
    # print(df_Left)


    df_HPED_axiom = pd.concat([df_Left, df_Right], axis=1)
    # print(df_HPED_axiom)
    df_HPED_axiom.to_csv(_out+'.hped', sep='\t', header=False, index=False)

    return _out+'.hped'



def _convert_HIBAG(_rhped, _out):


    if len(_rhped) != 8:
        print(std_ERROR_MAIN_PROCESS_NAME + "Platform \"{0}\" needs 8 input files. Please check it again.\n".format("HIBAG"))
        sys.exit()

    ### Exception handling 1 - When every rhped is 'NA'.

    f_isNA = list(map(lambda x : x == 'NA', _rhped))

    if all(f_isNA):
        print(std_ERROR_MAIN_PROCESS_NAME + "No any HIBAG output file has been given. Please check the '--rhped' argument again.")
        sys.exit()

    first_appear = -1

    for i in range(len(HLA_names)):
        if not f_isNA[i]:
            first_appear = i
            break

    DICT_rhped = {HLA_names[i]: pd.read_csv(_rhped[i], header=None, sep='\s+', dtype=str, names=["FID", "IID", "HLA", "2digit", "4digit", "p1", "p2"]) if not f_isNA[i] else pd.DataFrame([]) for i in range(0, len(HLA_names))}


    ### Exception handling 2 - When there is any rhped which has different number of rows.
    L = [DICT_rhped[HLA_names[i]].shape[0] for i in range(len(HLA_names)) if not f_isNA[i]]
    if len(set(L)) > 1:
        print(std_ERROR_MAIN_PROCESS_NAME + "There is an HIBAG output file which has different number of rows."
                                            "Please check the HIBAG output files given to the '--rhped' argument.")
        sys.exit()
    else:
        L = L.pop()


    ### The Main step to generate hped file.

    l_HLA = []  # the list to be pd.concat(..., axis=1)

    for i in range(len(HLA_names)):

        if not f_isNA[i]:
            df_temp = DICT_rhped[HLA_names[i]]['4digit'].str.extract(r'(\d{4,5}),(\d{4,5})')
        else:
            df_temp = pd.DataFrame([['0', '0'] for z in range(L)])

        df_temp.columns = [HLA_names[i]+'_1', HLA_names[i]+'_2']
        # print(df_temp)
        l_HLA.append(df_temp)


    df_RIGHT = pd.concat(l_HLA, axis=1)
    # print(df_RIGHT)



    ### Meta information(Left 6 columns)

    sr_PID = pd.Series(['0' for z in range(L)])
    sr_MID = pd.Series(['0' for z in range(L)])
    sr_Sex = pd.Series(['0' for z in range(L)])
    sr_Phe = pd.Series(['-9' for z in range(L)])

    df_Left = pd.concat([DICT_rhped[HLA_names[first_appear]][['FID', 'IID']], sr_PID, sr_MID, sr_Sex, sr_Phe], axis=1)
    # print(df_Left)

    df_HPED_HIBAG = pd.concat([df_Left, df_RIGHT], axis=1)
    df_HPED_HIBAG.to_csv(_out+'.hped', sep='\t', header=False, index=False)

    return _out + '.hped'



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

    # print(l_SampleID)
    # print(l_HLA)


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


    # print("\ndf_xHLA is \n")
    # print(df_RETURN.head())



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
    #                           "NA",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/Axiom/AxiomHLA_4dig_B_Results.txt",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/Axiom/AxiomHLA_4dig_C_Results.txt",
    #                           "NA",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/Axiom/AxiomHLA_4dig_DPB1_Results.txt",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/Axiom/AxiomHLA_4dig_DQA1_Results.txt",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/Axiom/AxiomHLA_4dig_DQB1_Results.txt",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/Axiom/AxiomHLA_4dig_DRB1_Results.txt",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/MyHLA2HPED_AXIOM/AXIOM"
    #                           ])

    # args = parser.parse_args(["-p", "HIBAG", "-rhped",
    #                           "NA",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-B.out",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-C.out",
    #                           "NA",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DPB1.out",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DQA1.out",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DQB1.out",
    #                           "/Users/wansun/Git_Projects/HATK/example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DRB1.out",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/MyHLA2HPED_HIBAG/HIBAG"
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