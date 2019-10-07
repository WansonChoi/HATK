#-*- coding: utf-8 -*-


import os, sys, re
import argparse, textwrap
import pandas as pd
from math import log10
from shutil import which


# Paths
p_Rscript = which("Rscript")



########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


class HATK_Heatmap(object):

    def __init__(self, _hla_name, _out, _p_maptable, _p_assoc_result, *args, **kwargs):


        if not _hla_name:
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--HLA' argument again.")
            sys.exit()

        if not _out:
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--out' argument again.")
            sys.exit()

        if not _p_maptable:
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--maptable' argument again.")
            sys.exit()


        t_single_assoc_result = None

        if not _p_assoc_result:
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--assoc-result/-ar' argument again.")
            sys.exit()
        else:
            # Supposed to be a list with only one element
            if isinstance(_p_assoc_result, list):

                if len(_p_assoc_result) > 1:
                    print(std_WARNING_MAIN_PROCESS_NAME + "More than 1 association test result was given.\n"
                                                          "Only 1st item will be used to plot HLA Heatmap.")

                t_single_assoc_result = _p_assoc_result[0]

            elif isinstance(_p_assoc_result, str):

                t_single_assoc_result = _p_assoc_result


        self.results = HEATMAP(_hla_name, _out, _p_maptable, t_single_assoc_result,
                               __save_intermediates=kwargs["__save_intermediates"], _p_src=kwargs["_p_src"],
                               _p_data=kwargs["_p_data"])


        self.removeIntermediates(_out)




    def getResults(self):
        return self.results


    def removeIntermediates(self, _out):

        # *.log
        if os.path.exists(_out+'.log'):
            os.system("rm {}".format(_out+'.log'))




def HEATMAP(_hla_name, _out, _p_maptable, _p_assoc_result, __save_intermediates=False, _p_Rscript=p_Rscript,
            _p_src="./src", _p_data="./data"):

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

    H_MARKERS = pd.DataFrame()
    __MAPTABLE__ = None
    __ASSOC_RESULT__ = pd.DataFrame()
    __ASSOC_RESULT_AA__ = pd.DataFrame()
    __ASSOC_RESULT_HLA__ = pd.DataFrame()


    # Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)


    _p_heatmapR = os.path.join(_p_src, "8b_plot_WS.R")



    ##### < Control Flags > #####

    LOADING_MAPTABLE = 1
    LOADING_ASSOC_RESULT = 1
    PREPROCESSING_MAPTABLE = 1

    MAKING_NEW_ASSOC = 1
    MAKING_ASSOC_P = 1
    EXPORTING_OUTPUT = 1
    PLOT_HEATMAP = 1




    if LOADING_MAPTABLE:

        ########## < [1] Loading HLA_MAPTABLE file > ##########

        # print(std_MAIN_PROCESS_NAME + "[0-1] Loading 'Maptable' file. ({})\n".format(_p_maptable))

        H_MARKERS = pd.read_csv(_p_maptable, sep='\s+', header=[0, 1, 2], index_col=0).filter(regex=_hla_name + "\*", axis=0)
        # filter() due to the case of "DRB2", "DRB3", etc.

        # print(H_MARKERS.head())
        # (2018. 6. 12) Anyway, `H_MARKERS` corresponds to "maptable" in Professor Han's pipeline.



    if LOADING_ASSOC_RESULT:

        ########## < [2] Loading association test result file(ex. *.assoc.logistic). > ##########

        # print(std_MAIN_PROCESS_NAME + "[0-2] Loading 'Association Test Result' file of AA and HLA markers. (Dropping 'NA')\n")

        __ASSOC_RESULT__ = pd.read_csv(_p_assoc_result, header=0, sep='\s+', usecols=["SNP", "A1", "OR", "P"]).dropna() # dropna() introduced (2019. 07. 02.)
        # print(__ASSOC_RESULT__.head())

        p_AA = re.compile(r"^AA_{0}_".format(_hla_name))
        p_HLA = re.compile(r"^HLA_{0}".format(_hla_name))

        f_AA = __ASSOC_RESULT__.loc[:, "SNP"].str.match(p_AA)
        f_HLA = __ASSOC_RESULT__.loc[:, "SNP"].str.match(p_HLA)


        __ASSOC_RESULT_AA__ = __ASSOC_RESULT__.loc[f_AA, :]
        __ASSOC_RESULT_HLA__ = __ASSOC_RESULT__.loc[f_HLA, :]

        # print("AA : \n{}\n".format(__ASSOC_RESULT_AA__.head()))
        # print("HLA : \n{}\n".format(__ASSOC_RESULT_HLA__.head()))


        # Checking subsetted logistic regression result.
        if __ASSOC_RESULT_AA__.shape[0] == 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "Given association result file({}) doesn't contain any Amino Acid markers. Please check it again.\n".format(_p_assoc_result))
            return -1

        if __ASSOC_RESULT_HLA__.shape[0] == 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "Given association result file({}) doesn't contain any HLA allele markers. Please check it again.\n".format(_p_assoc_result))
            return -1





    if PREPROCESSING_MAPTABLE:

        ########## < [3] Preprocessing MAPTABLE. > ##########

        # print(std_MAIN_PROCESS_NAME + "[1] Preprocessing MAPTABLE. ('*.map.txt')")

        """
        Maptable info && Association test result.
        
        [1] (Row) Excluding the HLA alleles that don't appear in given association test result.
        [2] (Column) Excluding non-polymorphic positions. (unrelated to association test)
        [3] (Column) Excluding the positions that don't appear in Amino acid markers of given association test result.
        
        """


        ### [1] (Row) Excluding the HLA alleles that don't appear in given association test result.

        HLA_alleles = __ASSOC_RESULT_HLA__.loc[:, "SNP"].apply(lambda x : re.sub(pattern=r'HLA_', repl='', string=x))
        # print("\nGiven HLA_alleles in association result: \n{0}\n".format(HLA_alleles))

        f_given_HLA_alleles = H_MARKERS.index.to_series().isin(HLA_alleles)

        __MAPTABLE__ = H_MARKERS.loc[f_given_HLA_alleles]
        # print(__MAPTABLE__)

        # __MAPTABLE__.to_csv(_out+".1st.txt", sep='\t', header=True, index=True)



        ### [2] (Column) Excluding non-polymorphic positions. (unrelated to association test)

        f_isPolymorphic = __MAPTABLE__.apply(lambda x : len(set(x)), axis=0) > 1

        __MAPTABLE__ = __MAPTABLE__.loc[:, f_isPolymorphic]
        # print(__MAPTABLE__)

        # __MAPTABLE__.to_csv(_out+".2nd.txt", sep='\t', header=True, index=True)



        ### [3] (Column) Excluding the positions that don't appear in Amino acid markers of given association test result.

        t_maptable_columns = pd.Series([item[0] for item in __MAPTABLE__.columns.tolist()])
        t_assoc_AA_MARKERS = __ASSOC_RESULT_AA__.loc[:, "SNP"]
        # print(t_maptable_columns.head())

        p_relPOS = re.compile(r'^AA_{}_(-?\d+)_\d+_.+'.format(_hla_name))
        t_assoc_relPOS = t_assoc_AA_MARKERS.str.extract(p_relPOS, expand=False)

        flag_valid_relPOS = t_maptable_columns.isin(t_assoc_relPOS).tolist()
        # print(flag_valid_relPOS)

        __MAPTABLE__ = __MAPTABLE__.loc[:, flag_valid_relPOS]

        # print("\nPreprocessed MAPTABLE.\n")
        # print(__MAPTABLE__.head())

        # __MAPTABLE__.to_csv(_out+".maptable.{}.txt".format(_hla_name), sep='\t', header=True, index=True)






        ### Reindexing `__ASSOC_RESULT_AA__` DataFrame with "relative_position".
        __ASSOC_RESULT_AA__.index = pd.Index(t_assoc_relPOS, name="rel_pos")  # This re-indexing will be used in next code block.
        # print(__ASSOC_RESULT_AA__.head())




    ########## < Main job to process maptable and Making new association > ##########




    if MAKING_NEW_ASSOC:

        ########## < [4] Making new *.assoc file for AA markers. > ##########

        # print(std_MAIN_PROCESS_NAME + "[2] Making new *.assoc file for AA markers. ('*.assoc.txt')\n")

        """
        Above preprocessed `__MAPTABLE__` will be the main content of Heatmap plot.
        (i.e. amino acid characters in each cell)

        In this block, P-value information corresponding to those amino acid characters will be obtained here.
        (i.e. the color of each of those amino acid characters)
        Besides, the sign will be given to those p-value depending on whether its OR value is above 1 or not.
        (ex. if OR is 0.278(<1), then its p-value will be -7.792000e-117(negative. to represent it is risky but not protective).



        Rememeber that `__MAPTABLE__` was subsetted to have only AA markers of which the relative position appears in `__ASSOC_RESULT_AA__`.

        As iterating over AA relative poisitions of `__MAPTABLE__` (ex. for i in [-27, -21, -18, ..., 232, 233, 234] <= filtered above.)
        and checking how many corresponding markers are in `__ASSOC_RESULT__`, the process is considering next two major cases.

        (1) Bi-allelic,
        (2) more than Tri-allelic.

        If the markers of `__ASSOC_RESULT_AA__` at each relative position is given as `df_BroughtMarkers`, then it would looks
        like this.


        ### More than bi-allelic ###

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


        (cf) Only single characters will be dealt. (i.e. The markers like 'AA_DQB1_71_32632639_exon2_AD' will be ignored.


        """


        # The main process will be done iterating the column of `__MAPTABLE__` DataFrame.
        COLNAMES = __MAPTABLE__.columns.tolist()


        for_new_assoc = []

        for i in range(0, __MAPTABLE__.shape[1]):
        # for i in range(0, 10):

            # print("=================================%d iteration=================================" % i, end='\r')

            # t_col_maptable := `AAname` (variable in R)
            t_col_maptable = COLNAMES[i] # ex) ('-27', '32634368', 'exon1')
            # print("\nAAname : {}".format(t_col_maptable))


            ### Bringing markers in each relative position to `df_BroughtMarkers`.
            df_BroughtMarkers = __ASSOC_RESULT_AA__.loc[[t_col_maptable[0]]]
            # print("\ndf_BroughtMarkers : \n{}".format(df_BroughtMarkers))



            if df_BroughtMarkers.shape[0] > 0:

                ### 1. P-value
                # AAvar_assoc (P-value column)
                AAvar_assoc = df_BroughtMarkers.loc[:, "P"]
                # print("\nAAvar_assoc(P-value) : \n")
                # print(AAvar_assoc)

                ### 2. AA character in `__MAPTABLE__` - (variable in R : `AAs`)
                AAs = __MAPTABLE__.iloc[:, i]
                # print("\nAAs")
                # print(AAs)

                ### 3. OR
                t_OR = df_BroughtMarkers.loc[:, "OR"]
                # print("\nOR")
                # print(t_OR)


                ### Preparing index for `AAvar_assoc`. ###

                if df_BroughtMarkers.shape[0] > 1:      ##### More than bi-allelic #####

                    ### Newly index `AAvar_assoc`.

                    # (variable in R) - (3) : `AAvariants`
                    AAvariants = df_BroughtMarkers.loc[:, "SNP"]

                    AAvar_assoc.index = pd.Index([item.split('_')[-1] for item in AAvariants.tolist()])
                    # ex) ['AA_DRB1_-25_32665484_K', 'AA_DRB1_-25_32665484_R', 'AA_DRB1_-25_32665484_x']
                    # => ['K', 'R', 'x']

                    # print("\nNew AAvar_assoc : \n")
                    # print(AAvar_assoc)

                    ### Transforming AA character seq. of `__MAPTABLE__` to '-log(x)' value.
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


                    I want to use the `df_BroughtMarkers` with that index(Only AA character).
                    Because, if that DataFrame is prepared, when `__MAPTABLE__` is given like this,

                    genomic_position  32557506 32557503 32557482 32557479 32557437 32557434  \
                    relative_position      -25      -24      -17      -16       -2       -1
                    codon                  AAG      CTC      ACA      GCG      TTG      GCT
                    DRB1*01:01:01            K        L        T        A        L        A
                    DRB1*03:01:01:01         R        L        A        V        L        A
                    DRB1*04:01:01            K        F        A        A        L        A
                    DRB1*04:03:01            K        F        A        A        L        A
                    DRB1*04:04:01            K        F        A        A        L        A
                    DRB1*04:05:01            K        F        A        A        L        A

                    We can make `__NEW_ASSOC__` just by this single command.

                    > AAvar_assoc.loc[x])

                    It will be much easier with .loc[] operator.

                    """

                elif df_BroughtMarkers.shape[0] == 1:   ##### Bi-allelic #####

                    ### Newly indexing `AAvar_assoc`.

                    refA = df_BroughtMarkers.loc[:, "A1"].iat[0]
                    AAvar_assoc.index = pd.Index([refA])

                    # print("\nNew AAvar_assoc : \n")
                    # print(AAvar_assoc)

                    ### Transforming AA character seq. of `__MAPTABLE__` to '-log(x)' value.
                    AAs2 = AAs.apply(lambda x: -log10(AAvar_assoc.loc[refA]))
                    # Same p-value will be given to two amino acid characters.
                    # Remember that association test on two factors generates one p-value.

                    ### Flippng based on OR
                    t_OR.index = AAvar_assoc.index
                    AAs3 = AAs.apply(lambda x : (2*int(t_OR.loc[refA] > 1) - 1)*(2*int(x == refA) - 1))



                ### Finally processed AA character seq. of `__MAPTABLE__`.
                AAs4 = AAs2*AAs3
                # print("\nAAs4 : \n")
                # print(AAs4)

                for_new_assoc.append(AAs4)


        __NEW_ASSOC__ = pd.concat(for_new_assoc, axis=1)
        # print(__NEW_ASSOC__.head())
        # __NEW_ASSOC__.to_csv(_out+".assoc2.{}.txt".format(_hla_name), sep='\t', header=True, index=True)



    if MAKING_ASSOC_P:

        ########## < [5] Making Assoc_P file. > ##########

        # print(std_MAIN_PROCESS_NAME + "[3] Making Assoc_P file. ('*.alleleP.txt')\n")

        """
        Obtaining 'p-values of HLA alleles' in association test result. 
        These p-value will be used to represent the color of HLA allele (Right part of Heatmap).
        
        1-field HLA allele names should be filtered out.
        
        """

        ### Excluding 1-field HLA alleles in `__ASSOC_RESULT_HLA__`.

        t_HLA_alleles_in_assoc = __ASSOC_RESULT_HLA__.loc[:, ["SNP", "P", "OR"]]

        # print("\nHLA alleles in Association Test result: ")
        # print(t_HLA_alleles_in_assoc)

        t_HLA_alleles_in_MAPTABLE = __MAPTABLE__.index.to_frame().apply(lambda x : "HLA_"+x).reset_index(drop=False)
        t_HLA_alleles_in_MAPTABLE.columns = pd.Index(["HLA1", "HLA2"])
        # print("\nHLA alleles in `__MAPTABLE__`: ")
        # print(t_HLA_alleles_in_MAPTABLE)


        df_No1Field = t_HLA_alleles_in_MAPTABLE.merge(t_HLA_alleles_in_assoc, left_on="HLA2", right_on="SNP", how="left").loc[:, ["HLA1", "OR", "P"]]
        # print(df_No1Field)

        sr_log10P = df_No1Field.loc[:, "P"].apply(lambda x : -log10(x))
        # print("\nsr_log10P : \n{}\n".format(sr_log10P))

        sr_OR = df_No1Field.loc[:, "OR"].apply(lambda x : (2*int(x > 1) -1))
        # print("\nsr_OR : \n{}\n".format(sr_OR))


        __alleleP__ = sr_log10P*sr_OR
        __alleleP__.index = df_No1Field.loc[:, "HLA1"]

        # print(__alleleP__)


        """
        (2019. 04. 09.)
        New major bug found.

        In the process of removing 1-field allele in above `__alleleP__`, indexes are mismatched.

        In next block, 2-field alleles are extracted in `__MAPTABLE__` dataframe but using these extracted 2-field alleles as an index for
        `__alleleP__` doesn't make sense. Consequently, Processed -log10(P) values in `__alleleP__` are given to wrong HLA alleles.

        So, indexing should be done to `__alleleP__` here, and 2-field alleles should be extracted from this `__alleleP__` dataframe.
        
        => Solved by using merge() function.

        """


    if EXPORTING_OUTPUT:

        ########## < [6] Exporting output files. > ##########

        # print(std_MAIN_PROCESS_NAME + "[4] Exporting output files. (Forwarding them to Rscript.)")

        """
        Files to export.
        
        (1) *.map.txt (differnt to the ones used by plink.)
        (2) *.allleP.txt
        (3) *.assoc.txt
        

        (2018. 8. 23.)
        In this code block, the proper index should be made and assigned to `__MAPTABLE__`, `__alleleP__` and `__NEW_ASSOC__`.       
        Also, as professor Han requested, the HLA allele names in final index should be in the form of 2-field.
        So, i need the marker labels in the next forms
        
        (1) AA_A_-?\d+_
        (2) HLA_A_\d{2,3}\:\d{2,3}
        
        """

        ##### Processing HLA markers.

        # if not __as4field:
        #     p_HLA = re.compile(r"^(%s\*\d{2,3}:\d{2,3}[A-Z]?).*" % (_hla_name))
        #     t_hla_markers = __MAPTABLE__.index.to_series().str.extract(p_HLA, expand=False)
        # else:
        #     t_hla_markers = __MAPTABLE__.index

        t_hla_markers = __MAPTABLE__.index


        ##### Processing AA markers.
        t_aa_markers = __MAPTABLE__.columns.to_frame().loc[:, "AA_rel_pos"]

        # print("\nt_hla_markers :\n{0}\nt_aa_markers :\n{1}\n".format(t_hla_markers, t_aa_markers))

        ##### Assigning processed indexes.

        # `__MAPTABLE__`
        __MAPTABLE__.index = t_hla_markers
        __MAPTABLE__.columns = t_aa_markers
        __MAPTABLE__.columns.name = None

        # `__alleleP__`
        __alleleP__.index = t_hla_markers

        # `__NEW_ASSOC__`
        __NEW_ASSOC__.index = t_hla_markers
        __NEW_ASSOC__.columns = t_aa_markers
        __NEW_ASSOC__.columns.name = None



        # print("\n(1) __MAPTABLE__\n")
        # print(__MAPTABLE__.head())
        #
        # print("\n(2) __alleleP__\n")
        # print(__alleleP__.head())
        #
        # print("\n(3) __NEW_ASSOC__\n")
        # print(__NEW_ASSOC__.head())



        __MAPTABLE__.to_csv(_out+'.map.txt', sep='\t', header=True, index=True)
        __alleleP__.to_csv(_out+".alleleP.txt", sep='\t', header=True, index=True)
        __NEW_ASSOC__.to_csv(_out+".assoc.txt", sep='\t', header=True, index=True)





    if PLOT_HEATMAP:

        ########## < [7] Plotting heatmap. > ##########

        # print(std_MAIN_PROCESS_NAME + "[7] Plotting heatmap.")

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
                            _out,
                            '1>{LOG} 2>{LOG}'.format(LOG=(_out+'.log'))])

        # print(command)
        os.system(command)


    if not __save_intermediates:

        l_remove = [".map.txt", ".alleleP.txt", ".assoc.txt"]

        for item in l_remove:

            command = ' '.join(["rm", _out+item])
            # print(command)
            os.system(command)



    __RESULTS__ = _out + ".pdf"

    if os.path.exists(__RESULTS__):
        return _out + ".pdf"
    else:
        print(std_WARNING_MAIN_PROCESS_NAME + "Heatmap failed. ('{}')".format(_out+'.pdf'))
        return -1


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

    parser.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(ex. *.assoc.logistic).\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nSave intermediate files.\n\n", action='store_true')




    ##### < for Test > #####

    # (2018. 10. 29.) HATK Integration test
    # args = parser.parse_args(["--HLA", "DQB1",
    #                           "-mt", "/Users/wansun/Git_Projects/HLA_Heatmap/data/HLA_MAPTABLE_DQB1.hg19.imgt3320.txt",
    #                           "-o", "tests/T1D_DQB1_test",
    #                           "-ar", "/Users/wansun/Git_Projects/HLA_Heatmap/data/example/20190327_WTCCC_T1D.assoc.logistic",
    #                           "--save-intermediates"
    #                           ])

    # args = parser.parse_args(["--HLA", "DQB1",
    #                           "-mt", "/Users/wansun/Git_Projects/HLA_Heatmap/data/HLA_MAPTABLE_DQB1.hg19.imgt3320.txt",
    #                           "-o", "tests/T1D_DQB1_test",
    #                           "-ar", "/Users/wansun/Git_Projects/HLA_Heatmap/data/example/20190327_WTCCC_T1D.assoc.logistic",
    #                           "--as4field",
    #                           "--save-intermediates"
    #                           ])

    # args = parser.parse_args(["--HLA", "DRB1",
    #                       "-mt", "data/HLA_MAPTABLE_DRB1.hg19.imgt3320.txt",
    #                       "-o", "tests/WTCCC_RA_DRB1/WTCCC_RA_DRB1",
    #                       "-ar", "/Users/wansun/Git_Projects/HLA_Heatmap/data/example/20190327_WTCCC_RA.assoc.logistic",
    #                       "--save-intermediates"
    #                       ])

    # args = parser.parse_args(["--HLA", "A",
    #                           "-mt", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/d20190701/HLA_MAPTABLE_A.hg19.imgt3320.txt",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/d20190701/Heatmap_Test_A",
    #                           "-ar", "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/d20190701/20190701_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.assoc.logistic",
    #                           "--save-intermediates"
    #                           ])



    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    ##### < Additional Argument Processing. > #####

    # main function execution.
    HEATMAP(_hla_name=args.HLA, _out=args.o, _p_maptable=args.maptable, _p_assoc_result=args.assoc_result,
            __save_intermediates=args.save_intermediates)
