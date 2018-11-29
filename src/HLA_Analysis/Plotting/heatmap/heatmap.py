#-*- coding: utf-8 -*-


import os, sys, re
import argparse, textwrap
import pandas as pd
import numpy as np
from mpmath import log10
from shutil import which


# Paths
p_Rscript = which("Rscript")



########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def HATK_HEATMAP(_hla_name, _out, _p_maptable, _p_assoc_logistic_AA, _p_assoc_logistic_HLA,
                 _4field=False, _oldv=False, _p_Rscript=p_Rscript,
                 _p_heatmapR="./src/HLA_Analysis/plot/heatmap/8b_plot_WS.R"):


    # _hla_name ("--HLA")
    if not bool(_hla_name):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"{0}\" wasn't given. Please check it again.\n".format("--HLA"))
        sys.exit()

        # Checking whether `_hla_name` belongs to ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'] is done in argparse codes.


    # _out ("--out")
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"{0}\" wasn't given. Please check it again.\n".format("--out"))
        sys.exit()


    # _p_maptable ("--maptable")
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"{0}\" wasn't given. Please check it again.\n".format("--maptable"))
        sys.exit()
    else:
        # HLA name and Maptable.
        p = re.compile(r"HLA_MAPTABLE_{0}".format(_hla_name))
        m = p.search(_p_maptable)
        print(_p_maptable)

        if not m:
            print(std_ERROR_MAIN_PROCESS_NAME + "Given maptable file(\"{0}\") doesn't match to given HLA gene(\"{1}\"). Please check them again.\n".format(_p_maptable, _hla_name))
            sys.exit()


    # Logistic regression results
    if not (bool(_p_assoc_logistic_AA) and bool(_p_assoc_logistic_HLA)):
        print(std_ERROR_MAIN_PROCESS_NAME + "Necessary logistic regression results weren't given properly. Please check them again.\n(\"{0}\",\"{1}\")".format(_p_assoc_logistic_AA, _p_assoc_logistic_HLA))
        sys.exit()






    return HEATMAP(_hla_name, _out, _p_maptable, _p_assoc_logistic_AA, _p_assoc_logistic_HLA,
                 _p_assoc_logistic_MERGED=None, _4field=False, _oldv=False, _p_Rscript=p_Rscript,
                 _p_heatmapR="src/HLA_Analysis/Plotting/heatmap/8b_plot_WS.R")






def HEATMAP(_hla_name, _out, _p_maptable, _p_assoc_logistic_AA, _p_assoc_logistic_HLA, _p_assoc_logistic_MERGED=None,
            _4field=False, _oldv=False, _p_Rscript=p_Rscript, _No_Intermediates=True,
            _p_heatmapR="src/HLA_Analysis/Plotting/heatmap/8b_plot_WS.R"):

    """

    This script is replacement of '7_prepare.R'.
    The main job of that script is to subset and filter HLA marker data given as file '6_table.txt' which is equivalent to HLA marker dicitionary file.

    Finally, after this script implemented, 3 data file will be made.

        (1) *_map.txt
        (2) *_assoc.txt
        (3) *_alleleP.txt


    (2018. 8. 20.)
    The argument "_4field" and "_oldv" are deprecated.

    """

    ########## < Core Variables > ##########

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting HeatMap Plotting.\n\n")

    H_MARKERS = pd.DataFrame()
    BIM_HLA = pd.DataFrame()
    BIM_AA = pd.DataFrame()
    ASSOC_LOGISTIC = pd.DataFrame()

    # MARKERS = pd.DataFrame()


    # Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)






    ##### < Control Flags > #####

    LOADING_HLA_MARKERS_DICTIONARY = 1
    LOADING_ASSOC = 1
    SUBSETTING_MAPTABLE1 = 1
    SUBSETTING_MAPTABLE2 = 1

    MAKING_NEW_ASSOC = 1
    MAKING_ASSOC_P = 1
    EXPORTING_OUTPUT = 1
    PLOT_HEATMAP = 1




    if LOADING_HLA_MARKERS_DICTIONARY:

        ########## < [1] Loading HLA_MARKERS_DICTIONARY("maptable") file > ##########
        print(std_MAIN_PROCESS_NAME + "[1] Loading HLA_MARKERS_DICTIONARY(\"maptable\") file\n")

        H_MARKERS = pd.read_table(_p_maptable, sep=' |\t', engine='python', header=[0, 1, 2], index_col=0).filter(regex=_hla_name + "\*", axis=0)
        # filter() due to the case of "DRB2", "DRB3", etc.

        print(H_MARKERS.head())

        # So far, Only by loading HLA marker dictionary file, we finished preparing `maptable`.
        # (2018. 6. 12) Anyway, `H_MARKERS` corresponds to "maptable" in Professor Han's pipeline.



    if LOADING_ASSOC:

        ########## < [2] Loading *.assoc.logistic files of AA and HLA markers. > ##########

        print(std_MAIN_PROCESS_NAME + "Loading '*.assoc.logistc' files of AA and HLA markers.\n")

        # ASSOC_LOGISTIC = pd.read_table(_p_assoc_logistic, header=0, sep='\s+', usecols=["SNP", "P", "OR", "A1"]
        #                                ).set_index("SNP", drop=False).filter(regex="_".join(["AA", _hla_name]), axis=0).reset_index(drop=True)

        p_AA_marker = re.compile(r"AA_{0}_".format(_hla_name))
        p_HLA_marker = re.compile(r"HLA_{0}".format(_hla_name))

        # Just classify the cases where *.assoc.logistic files are given separately or not.
        if (bool(_p_assoc_logistic_AA) and bool(_p_assoc_logistic_HLA)) and not bool(_p_assoc_logistic_MERGED):

            # Respective
            print("\nLoading \"*.assoc.logistic\" of AA and HLA markers.\n")

            ASSOC_LOGISTIC_AA = pd.read_table(_p_assoc_logistic_AA, header=0, sep='\s+',
                                              usecols=["SNP", "BP", "A1", "OR", "P"]
                                              ).set_index("SNP").filter(regex=p_AA_marker, axis=0).reset_index(drop=False)

            ASSOC_LOGISTIC_HLA = pd.read_table(_p_assoc_logistic_HLA, header=0, sep='\s+',
                                               usecols=["SNP", "BP", "A1", "OR", "P"]
                                               ).set_index("SNP").filter(regex=p_HLA_marker, axis=0).reset_index(drop=False)


        elif not (bool(_p_assoc_logistic_AA) or bool(_p_assoc_logistic_HLA)) and bool(_p_assoc_logistic_MERGED):

            # When the results of logistic regression for AA and HLA are given as merged single file.

            print("\nLoading \"*.assoc.logistic\" of AA and HLA markers as single file.\n")

            df_temp = pd.read_table(_p_assoc_logistic_MERGED, header=0, sep='\s+',
                                    usecols=["SNP", "BP", "A1", "OR", "P"]).set_index("SNP")

            ASSOC_LOGISTIC_AA = df_temp.filter(regex="AA_{0}_".format(_hla_name)).reset_index(drop=False)
            ASSOC_LOGISTIC_HLA = df_temp.filter(regex="HLA_{0}".format(_hla_name)).reset_index(drop=False)


        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Check the file names of \"*.assoc.logistic.\"")


        # Checking subsetted logistic regression result.
        if ASSOC_LOGISTIC_AA.shape[0] == 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "The subsetted logistic regression result of AA markers has 0 rows(0 markers). Please check next 2 possible cases.\n"
                                                "(1) Given logistic regression file for AA markers(\"{0}\") has 0 markers at first, or\n"
                                                "(2) Given logistic regression file for AA markers doesn't have markers which belongs to target HLA gene(\"{1}\"). "
                                                .format(_p_assoc_logistic_HLA, _hla_name))
            return -1

        if ASSOC_LOGISTIC_HLA.shape[0] == 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "The subsetted logistic regression result of HLA markers has 0 rows(0 markers). Please check next 2 possible cases.\n"
                                                "(1) Given logistic regression file for HLA markers(\"{0}\") has 0 markers at first, or\n"
                                                "(2) Given logistic regression file for HLA markers doesn't have markers which belongs to target HLA gene(\"{1}\"). "
                                                .format(_p_assoc_logistic_HLA, _hla_name))

            return -1



        HLA_MARKERS = ASSOC_LOGISTIC_HLA.loc[:, "SNP"]
        AA_MARKERS = ASSOC_LOGISTIC_AA.loc[:, "SNP"]


        print("\nLogistic Regression of AA :")
        print(ASSOC_LOGISTIC_AA.head())
        # print(ASSOC_LOGISTIC_AA.tail())

        print("\nLogistic Regression of HLA :")
        print(ASSOC_LOGISTIC_HLA.head())
        # print(ASSOC_LOGISTIC_HLA.tail())



    if SUBSETTING_MAPTABLE1:

        ########## < [2] Subsetting maptable ((1) HLA alleles(index), (2) Relative p(column)) > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Subsetting maptable ((1) HLA alleles(index), (2) Relative p(column))")


        ### Subsetting HLA dictionary(`H_MARKERS`) DataFrame

        # (1st Filtering condition - Row) Find the markers of *.bim file in HLA dict(`H_MARKERS`) DataFrame.

        p = re.compile(r'HLA_')

        HLA_alleles = HLA_MARKERS.apply(lambda x : p.sub(repl='', string=x)).tolist()
        print("\nHLA_alleles in \"*.assoc.logistic\" of HLA markers : \n{0}\n".format(HLA_alleles))

        flag_LABELSinDICT = H_MARKERS.index.to_series().isin(HLA_alleles)

        ### Filtering `H_MARKERS` with 1st condition.
        sub_H_MARKERS = H_MARKERS.loc[flag_LABELSinDICT]

        # sub_H_MARKERS.to_csv(_out+".1st.txt", sep='\t', header=True, index=True)


        # (2nd Filtering Condition 2 - Column) Find the spots which have more than equal 2 kind of AA markers.

        flag_isPolymorphic = sub_H_MARKERS.apply(lambda x : len(set(x)), axis=0) > 1

        ### Filtering `H_MARKERS` with 2nd condition.
        sub_H_MARKERS = sub_H_MARKERS.loc[:, flag_isPolymorphic]


        print(sub_H_MARKERS.head())
        print(sub_H_MARKERS.tail())

        # sub_H_MARKERS.to_csv(_out+".2nd.txt", sep='\t', header=True, index=True)




    if SUBSETTING_MAPTABLE2:

        ########## < [3] Subsetting maptable ((3) AA marker of "*.assoc.logistic" of AA) > ##########

        print(std_MAIN_PROCESS_NAME + "[3] Subsetting maptable ((3) AA marker of \"*.assoc.logistic\" of AA)")

        sub_H_MARKERS_columns = pd.Series(["_".join(["AA", _hla_name, item[0]]) for item in sub_H_MARKERS.columns.tolist()]) # (*****) This part is core to subset
        # Using "relative_position" information to extract markers in both "sub_H_MARKERS" and "ASSOC_LOGISTIC".


        # (3rd Filtering condition - Overlapping relative position)

        p_relPOS = re.compile("(AA_{0}_-?\d+)".format(_hla_name))

        flag_valid_relPOS = sub_H_MARKERS_columns.isin(AA_MARKERS.str.extract(p_relPOS, expand=False).tolist()).tolist()


        ### Filtering `H_MARKERS` with 3rd condition.
        sub_H_MARKERS = sub_H_MARKERS.loc[:, flag_valid_relPOS]

        print("\nFinally subsetted HLA dictionary file(\"maptable\").\n")
        print(sub_H_MARKERS.head())

        # sub_H_MARKERS.to_csv(_out+".{0}.3rd.txt".format(_hla_name), sep='\t', header=True, index=True)


        ### Reindexing `ASSOC_LOGISTIC_AA` DataFrame with "relative_position".

        p_relPOS = re.compile("AA_{0}_(-?\d+)".format(_hla_name))

        ASSOC_LOGISTIC_AA.index = pd.Index(AA_MARKERS.str.extract(p_relPOS, expand=False).tolist())

        # This re-indexing will be used in next code block.

        print(ASSOC_LOGISTIC_AA.head())




    ########## < Main job to process maptable and Making new association > ##########




    if MAKING_NEW_ASSOC:

        ########## < [4] Making new *.assoc file for AA markers. > ##########


        print(std_MAIN_PROCESS_NAME + "Making new *.assoc file for AA markers.")

        # The main process will be done iterating the column of `sub_H_MARKERS` DataFrame.
        COLNAMES = sub_H_MARKERS.columns.tolist() # cf) len(COLNAMES) == len(col of `sub_H_MARKERS`)


        for_new_assoc = []

        for i in range(0, sub_H_MARKERS.shape[1]):
        # for i in range(0, 5):

            print("=================================%d iteration=================================" % i, end='\r')

            # (variable in R) - (2) : `AAname`
            AAname = COLNAMES[i] # ex) ('-25', '32557506', 'AAG')
            # print("\nAAname")
            # print(AAname)


            """
            (2018. 8. 20.)
            Remember that there could be more than 1 marker in one relative position.
            That's why one more searching job should be done to `ASSOC_LOGISTIC_AA` DataFrame.
            
            ex) HLA-A, rel:-15
            1. AA_A_-15_29910359_L
            2. AA_A_-15_29910359_V
            3. AA_A_-15_29910359_x

            => Three markers in -15 relative position of Amino Acids.
            
            """

            ### Bringing markers in each relative position to `df_temp`.
            # (2018. 8. 20.) This temporary DataFrame `df_temp` will do a central role in each iteration.
            df_temp = ASSOC_LOGISTIC_AA.loc[[AAname[0]]]
            # print("\ndf_temp : \n")
            # print(df_temp)



            if df_temp.shape[0] > 0:

                """
                Introducing this condition (len(AAvariants) > 0) might seem trivial, but it filters unexpected markers
                which don't appear in "*.assoc.logistic(`ASSOC_LOGISTIC_AA`)" but in "HLA dictionary(`sub_H_MARKERS`).

                cf) len(AAvariants) == df_temp.shape[0] == The number of rows of `df_temp`


                As iterating over AA poisition of "HLA dictionary(`sub_H_MARKERS`)" and checking how many markers of
                "*.assoc.logistc(`ASSOC_LOGISTIC`)" there are, the process should do it considering next two major cases.

                (1) Bi-allelic,
                (2) more than Tri-allelic.

                If the markers of `ASSOC_LOGISTIC_AA` at each relative position is given as `df_temp`, then it would looks
                like this.


                ### Tri-allelic ###

                AAvar_assoc
                SNP
                AA_DRB1_-25_32665484_K    0.3133
                AA_DRB1_-25_32665484_R    0.3173
                AA_DRB1_-25_32665484_x    0.9150
                Name: P, dtype: float64


                ### Bi-allelic ###

                AAvar_assoc
                SNP
                AA_DRB1_-25_32665484    0.3133    # The bi-allelic marker doesn't have AA character('K', 'R', ...) in the label.
                Name: P, dtype: float64

                """

                ### 1. P-value
                # AAvar_assoc (P-value column)
                AAvar_assoc = df_temp.loc[:, "P"]
                # print("\nAAvar_assoc(P-value) : \n")
                # print(AAvar_assoc)

                ### 2. AA character in "HLA dictionary(sub_H_MARKERS)" - (variable in R : `AAs`)
                AAs = sub_H_MARKERS.iloc[:, i]
                # print("\nAAs")
                # print(AAs)

                ### 3. OR
                t_OR = df_temp.loc[:, "OR"]
                # print("\nOR")
                # print(t_OR)


                ### Preparing index for `AAvar_assoc`. ###

                if df_temp.shape[0] > 1:

                    ##### Tri-allelic or more #####

                    ### Newly index `AAvar_assoc`.

                    # (variable in R) - (3) : `AAvariants`
                    AAvariants = df_temp.loc[:, "SNP"]

                    AAvar_assoc.index = pd.Index([item.split('_')[-1] for item in AAvariants.tolist()])
                    # ex) ['AA_DRB1_-25_32665484_K', 'AA_DRB1_-25_32665484_R', 'AA_DRB1_-25_32665484_x']
                    # => ['K', 'R', 'x']

                    # print("\nNew AAvar_assoc : \n")
                    # print(AAvar_assoc)

                    ### Transforming AA character seq. of `sub_H_MARKERS` to '-log(x)' value.
                    AAs2 = AAs.apply(lambda x : -log10(AAvar_assoc.loc[x]))
                    # This part actually needs exception handling for KeyError...

                    ### Flippng based on OR
                    t_OR.index = AAvar_assoc.index
                    AAs3 = AAs.apply(lambda x : (2*int(t_OR.loc[x] > 1) - 1))


                    """
                    So, when Tri-allelic or more is given like this,

                    SNP
                    AA_DRB1_-25_32665484_K    0.3133
                    AA_DRB1_-25_32665484_R    0.3173
                    AA_DRB1_-25_32665484_x    0.9150
                    Name: P, dtype: float64


                    I will make this to the next form.

                    K    0.3133
                    R    0.3173
                    x    0.9150
                    Name: P, dtype: float64


                    I want to use the `df_temp` with that index(Only AA character).
                    Because, if that DataFrame is prepared, when `sub_H_MARKERS` is given like this,

                    genomic_position  32557506 32557503 32557482 32557479 32557437 32557434  \
                    relative_position      -25      -24      -17      -16       -2       -1
                    codon                  AAG      CTC      ACA      GCG      TTG      GCT
                    DRB1*01:01:01            K        L        T        A        L        A
                    DRB1*03:01:01:01         R        L        A        V        L        A
                    DRB1*04:01:01            K        F        A        A        L        A
                    DRB1*04:03:01            K        F        A        A        L        A
                    DRB1*04:04:01            K        F        A        A        L        A
                    DRB1*04:05:01            K        F        A        A        L        A

                    We can make `NEW_ASSOC` just by this single command.

                    > AAvar_assoc.loc[x])

                    It will be much easier with .loc[] operator.

                    """

                elif df_temp.shape[0] == 1:

                    ##### Bi-allelic #####

                    ### Newly indexing `AAvar_assoc`.

                    refA = df_temp.loc[:, "A1"].iat[0]
                    AAvar_assoc.index = pd.Index([refA])

                    # print("\nNew AAvar_assoc : \n")
                    # print(AAvar_assoc)

                    ### Transforming AA character seq. of `sub_H_MARKERS` to '-log(x)' value.
                    AAs2 = AAs.apply(lambda x: -log10(AAvar_assoc.loc[refA]))

                    ### Flippng based on OR
                    t_OR.index = AAvar_assoc.index
                    AAs3 = AAs.apply(lambda x : (2*int(t_OR.loc[refA] > 1) - 1)*(2*(x == refA)) - 1)



                ### Finally processed AA character seq. of `sub_H_MARKERS`.
                AAs4 = AAs2*AAs3
                # print("\nAAs4 : \n")
                # print(AAs4)

                for_new_assoc.append(AAs4)


        NEW_ASSOC = pd.DataFrame(for_new_assoc).transpose()

        # print("\nNEW_ASSOC : \n{0}\n\n".format(NEW_ASSOC.head()))



    if MAKING_ASSOC_P:

        ########## < [5] Making Assoc_P file. > ##########

        print(std_MAIN_PROCESS_NAME + "Making Assoc_P file.\n\n")

        print("\nHLA markers :\n")
        print(HLA_MARKERS.head())

        flag_subHLAs = sub_H_MARKERS.index.to_series().apply(lambda x : "HLA_"+x).isin(HLA_MARKERS).tolist()

        sub_ASSOC_LOGISTIC_HLA = ASSOC_LOGISTIC_HLA.loc[flag_subHLAs, ["P", "OR"]]

        HLA_P = sub_ASSOC_LOGISTIC_HLA.loc[:, "P"].values # as numpy array
        print("\n\"P-values\" of HLA markers : \n")
        print(HLA_P)
        print(type(HLA_P))

        HLA_OR = sub_ASSOC_LOGISTIC_HLA.loc[:, "OR"].values # as numpy array
        print("\n\"OR\" of HLA markers : \n")
        print(HLA_OR)


        t_alleleP = np.apply_along_axis(lambda x : -np.log10(x), 0, HLA_P)
        t_OR = np.apply_along_axis(lambda x : (2*(x > 1) - 1), 0, HLA_OR)

        ### Processed "alleleP" for HLA markers.
        alleleP = t_alleleP*t_OR

        print("\nprocessed `alleleP` : \n")
        print(alleleP)




    if EXPORTING_OUTPUT:

        ########## < [6] Exporting output files. > ##########

        print(std_MAIN_PROCESS_NAME + "[6] Exporting output files.")

        """
        Files to export.
        
        (1) *.map.txt (differnt to the ones used by plink.)
        (2) *.allleP.txt
        (3) *.assoc.txt
        

        (2018. 8. 23.)
        In this code block, the proper index should be made and assigned to `sub_H_MARKERS`, `alleleP` and `NEW_ASSOC`.       
        Also, as professor Han requested, the HLA allele names in final index should be in the form of 2-field.
        So, i need the marker labels in the next forms
        
        (1) AA_A_-?\d+_
        (2) HLA_A_\d{2,3}\:\d{2,3}
        
        """

        ##### Processing HLA markers.
        p_HLA_marker = re.compile("(%s\*\d{2,3}\:\d{2,3}\w?)" % (_hla_name))
        t_hla_markers = sub_H_MARKERS.index.to_series().str.extract(p_HLA_marker, expand=False)


        ##### Processing AA markers.
        # t_aa_markers = sub_H_MARKERS.columns.to_frame().loc[:, "relative_position"]
        t_aa_markers = sub_H_MARKERS.columns.to_frame().loc[:, "AA_rel_pos"]

        print("\nt_hla_markers :\n{0}\nt_aa_markers :\n{1}\n".format(t_hla_markers, t_aa_markers))

        ##### Assigning processed indexes.

        # `sub_H_MARKERS`
        sub_H_MARKERS.index = t_hla_markers
        sub_H_MARKERS.columns = t_aa_markers
        sub_H_MARKERS.columns.name = None

        # `alleleP`
        alleleP = pd.DataFrame(alleleP, index=t_hla_markers)
        alleleP.columns.name = "0" # dummy index to make it loaded in R more efficiently.

        # `NEW_ASSOC`
        NEW_ASSOC.index = t_hla_markers
        NEW_ASSOC.columns = t_aa_markers
        NEW_ASSOC.columns.name = None



        print("\n(1) sub_H_MARKERS\n")
        print(sub_H_MARKERS.head())

        print("\n(2) alleleP\n")
        print(alleleP.head())

        print("\n(3) NEW_ASSOC\n")
        print(NEW_ASSOC.head())



        sub_H_MARKERS.to_csv(_out+'.map.txt', sep='\t', header=True, index=True)
        alleleP.to_csv(_out+".alleleP.txt", sep='\t', header=True, index=True)
        NEW_ASSOC.to_csv(_out+".assoc.txt", sep='\t', header=True, index=True)





    if PLOT_HEATMAP:

        ########## < [7] Plotting heatmap. > ##########

        print(std_MAIN_PROCESS_NAME + "[7] Plotting heatmap.")

        """
        (Argument Example for "8b_plot_WS.R")

        # static argument preparation for testing(set1 : "UC(Ulcerative Colitis)")
        # args1.disease.map_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_testws_UC_map.txt"
        # args2.disease.assoc_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_testws_UC_assoc.txt"
        # args3.disease.alleleP_ = "/Uers/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_testws_UC_alleleP.txt"
        # args4.HLA_name_ = "DRB1"
        # args5.plot.outf_ = "8_testws_UC"

        """

        command = ' '.join([p_Rscript,
                            _p_heatmapR,
                            _out+".map.txt",
                            _out+".assoc.txt",
                            _out+".alleleP.txt",
                            _hla_name,
                            _out])

        print(command)

        os.system(command)


    if _No_Intermediates:

        l_remove = [".map.txt", ".alleleP.txt", ".assoc.txt"]

        for item in l_remove:

            command = ' '.join(["rm", _out+item])
            print(command)
            os.system(command)



    return _out + ".pdf"


if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################
    
        heatmap.py
    
        - Originally, Plotting heatmap entailed too many steps based on bash, python and R. Now it is
        integrated to this python script. It will be used in "HLA_Analysis.py" script.
    
    #################################################################################################
                                     '''),
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-o", help="\nOutput file name prefix.\n\n", required=True)

    parser.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n", choices = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'])
    parser.add_argument("--maptable", "-mt", help="\nMarker Dictionary file(Maptable) generated by 'IMGTt2Sequence'.\n\n")
    parser.add_argument("--result-assoc-AA", "-lrA", help="\nAA logistic regression result file.\n\n")
    parser.add_argument("--result-assoc-HLA", "-lrH", help="\nHLA logistic regression result file.\n\n")
    parser.add_argument("--p_heatmapR", help="\nheatmap R source.\n\n")




    ##### < for Test > #####

    # (2018. 8. 20.)
    # args = parser.parse_args(["-hla", "DRB1",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_in_HATK2",
    #                           "-hd", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/WRAPER_TEST_DRB1.AA.markers.trim.labeled.txt",
    #                           "-lr", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/forHEATMAP_test.assoc.logistic",
    #                           "-bm", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/merged"
    #                           ])

    # # No more *.bim file.
    # args = parser.parse_args(["-hla", "A",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer",
    #                           "-hd", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/forHATK_A.AA.markers.trim.labeled.txt",
    #                           "-lrA", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.AA.CODED.assoc.logistic",
    #                           "-lrH", "/Users/wansun/Data/HATK/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.HLA.assoc.logistic",
    #                           ])


    # args = parser.parse_args(["-hla", "C",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer_C",
    #                           "-hd", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/forHATK_C.AA.markers.trim.labeled.txt",
    #                           "-lrA", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.AA.CODED.assoc.logistic",
    #                           "-lrH", "/Users/wansun/Data/HATK/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.HLA.assoc.logistic",
    #                           ])

    # # (2018. 10. 29.) HATK Integration test
    # args = parser.parse_args(["-hla", "A",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.A",
    #                           "-hd", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/HLA_MAPTABLE.A.hg18.imgt370.txt",
    #                           "-lrA", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.AA.CODED.assoc.logistic",
    #                           "-lrH", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.HLA.assoc.logistic",
    #                           ])


    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)
    # print(vars(args))


    ##### < Additional Argument Processing. > #####



    # # main function execution. (old)
    # HEATMAP(_hla_name=args.hla, _out=args.o, _p_maptable=args.hla_dict,
    #         _p_assoc_logistic_AA=args.logistic_result_AA, _p_assoc_logistic_HLA=args.logistic_result_HLA,
    #         _p_assoc_logistic_MERGED=args.logistic_result_MERGED,
    #         _p_heatmapR=args.p_heatmapR,
    #         _No_Intermediates=False)

    # main function execution.
    HEATMAP(_hla_name=args.HLA, _out=args.o, _p_maptable=args.maptable,
            _p_assoc_logistic_AA=args.result_assoc_AA, _p_assoc_logistic_HLA=args.result_assoc_HLA,
            _p_heatmapR=args.p_heatmapR,
            _No_Intermediates=False)
