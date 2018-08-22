#-*- coding: utf-8 -*-


import os, sys, re
import argparse, textwrap
import pandas as pd
from mpmath import log10
from shutil import which

########## <Core Global Variable> ##########

# Paths
p_Rscript = which("Rscript")

def HEATMAP(_hla_name, _out, _input,
            _p_hla_dict, _p_assoc_logistic,
            _bim_HLA = "Not_given", _bim_AA = "Not_given", _bim_merged = "Not_given",
            _4field=False, _oldv=False,
            _p_Rscript=p_Rscript, _p_heatmapR="./src/HLA_Analysis/plot/heatmap/8b_plot_WS.R"):

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


    ########## < Additional Argument Checking. > ##########

    # When "_input" argument is given
    if bool(_input):

        # Override every prefix of dependent files.
        _p_hla_dict = _input+".txt"
        _p_assoc_logistic = _input+".assoc.logistic"
        _bim_merged = _input+".bim" # It is natural that "*.AA.bim" and "*.HLA.bim" are given as merged .bim file.






    ##### < Control Flags > #####

    LOADING_HLA_MARKERS_DICTIONARY = 1
    LOADING_BIM = 1
    LOADING_ASSOC = 1

    MAKING_NEW_ASSOC = 1
    MAKING_ASSOC_P = 1
    EXPORTING_OUTPUT = 0
    PLOT_HEATMAP = 0




    ########## < Loading and Preparing Core Data Files given by arguments. > ##########


    if LOADING_HLA_MARKERS_DICTIONARY:

        ### Loading HLA_MARKERS_DICTIONARY file.
        print(std_MAIN_PROCESS_NAME + "Loading HLA marker information.\n")

        H_MARKERS = pd.read_table(_p_hla_dict, sep=' |\t', engine='python', header=[0, 1, 2], index_col=0).filter(regex=_hla_name + "\*", axis=0)
        # filter() due to the case of "DRB2", "DRB3", etc.

        print(H_MARKERS.head())

        # So far, Only by loading HLA marker dictionary file, we finished preparing `maptable`.
        # (2018. 6. 12) Anyway, `H_MARKERS` corresponds to "maptable" in Professor Han's pipeline.

        print("\n")


    if LOADING_BIM:

        ### Loading "*.bim" file.

        print(std_MAIN_PROCESS_NAME + "Loading *.bim file.\n")

        # patterns which will be used.
        p_HLA = re.compile("^HLA_{0}".format(_hla_name))

        # Just classify the cases where bim files are given separately or not.
        if (_bim_HLA != "Not_given" and _bim_AA != "Not_given") and _bim_merged == "Not_given":

            # "_bim_HLA" and "_bim_AA" are given.
            print("\"_bim_HLA\" : {0}".format(_bim_HLA))

            # Loading each of bim file
            BIM_HLA = pd.read_table(_bim_HLA, sep='\t| ', engine='python', header=None,
                                    names=["Chr", "Label", "GD", "Pos", "Al1", "Al2"],
                                    usecols=[1,3,4,5]).set_index("Label", drop=False).filter(regex=p_HLA, axis=0).reset_index(drop=True)


            print("\nDataFrame \"BIM_HLA\"\n")
            print(BIM_HLA.head())
            print(BIM_HLA.tail())


        # *.bim file for HLA markers has been loaded.


        ### Subsetting HLA dictionary(`H_MARKERS`) DataFrame

        # (1st Filtering condition - Row) Find the markers of *.bim file in HLA dict(`H_MARKERS`) DataFrame.

        p = re.compile(r'HLA_')

        LABELS = BIM_HLA.loc[:, "Label"].apply(lambda x : p.sub(repl='', string=x)).tolist()
        print("\nLABELS : \n{0}\n".format(LABELS))

        flag_LABELSinDICT = H_MARKERS.index.to_series().isin(LABELS)

        ### Filtering `H_MARKERS` with 1st condition.
        sub_H_MARKERS = H_MARKERS.loc[flag_LABELSinDICT]

        # sub_H_MARKERS.to_csv(_out+".1st.txt", sep='\t', header=True, index=True)


        # (2nd Filtering Condition 2 - Column) Find the spots which have more than equal 2 kind of AA markers.

        flag_Marker_Count = sub_H_MARKERS.apply(lambda x : len(set(x)), axis=0) > 1

        ### Filtering `H_MARKERS` with 2nd condition.
        sub_H_MARKERS = sub_H_MARKERS.loc[:, flag_Marker_Count]

        # `sub_H_MARKERS` == "maptable" (***)


        print(sub_H_MARKERS.head())
        print(sub_H_MARKERS.tail())

        # sub_H_MARKERS.to_csv(_out+".2nd.txt", sep='\t', header=True, index=True)



    if LOADING_ASSOC:

        ### Loading *.assoc.logistic. (The reuslt of Logistic Regression on Amino Acid markers)

        print(std_MAIN_PROCESS_NAME + "Loading '*.assoc.logistc' file.\n\n")

        ASSOC_LOGISTIC = pd.read_table(_p_assoc_logistic, header=0, sep='\s+', usecols=["SNP", "P", "OR", "A1"]
                                       ).set_index("SNP", drop=False).filter(regex="_".join(["AA", _hla_name]), axis=0).reset_index(drop=True)

        print(ASSOC_LOGISTIC.head())

        # In general, the marker set of HLA_Dicitionary(`sub_H_MARKERS`) would be larger than that of `ASSOC_LOGISTIC`.
        # So, the markers of `sub_H_MARKERS` will be subsetted based on the markers in `ASSOC_LOGISTIC`.

        # Also, this code block might deal with the result of "Omnibus_Test" too. Ask professor about this.

        # ASSOC의 marker set에서 나타나는지...

        # (1) `variants`
        ASSOC_MARKER_SET = ASSOC_LOGISTIC.loc[:, "SNP"]
        print("\nASSOC_MARKER_SET")
        print(ASSOC_MARKER_SET.head())

        # (2) `assocrst`
        assocrst = ASSOC_LOGISTIC.loc[:, "P"]
        print("\nassocrst")
        print(assocrst.head())

        # # Making Row index which will be used in next 'MAKING_NEW_ASSOC' code block.
        # MARKER_CHARACTERS = [ASSOC_MARKER_SET.apply(lambda x : x.split('_')[-1]).iat[i] for i in range(0, ASSOC_MARKER_SET.shape[0])]
        # idx_asscrst = pd.MultiIndex.from_arrays([ASSOC_MARKER_SET.index.tolist(), MARKER_CHARACTERS], names=["LABELS", "AA_CHARACTER"]) # (Marker_Full_Name, Actual_Marker_Character)
        # assocrst.index = idx_asscrst

        # (3) `OR`
        OR = ASSOC_LOGISTIC.loc[:, "OR"]
        # OR.index = idx_asscrst
        print("\nOR")
        print(OR.head())

        # (4) `A1`
        A1 = ASSOC_LOGISTIC.loc[:, "A1"]
        print("\nA1")
        print(A1.head())



        ### Subsetting HLA dictionary(`H_MARKERS`) DataFrame

        sub_H_MARKERS_columns = pd.Series(["_".join(["AA", _hla_name, item[1]]) for item in sub_H_MARKERS.columns.tolist()]) # (*****) This part is core to subset
        # Using "relative_position" information to extract markers which appear in both "sub_H_MARKERS" and "ASSOC_LOGISTIC".

        # (3rd Filtering condition - Overlapping relative position)

        p_relPOS = re.compile("(AA_{0}_-?\d+)".format(_hla_name))

        flag_valid_relPOS = sub_H_MARKERS_columns.isin(ASSOC_MARKER_SET.str.extract(p_relPOS, expand=False).tolist()).tolist()


        ### Filtering `H_MARKERS` with 3rd condition.
        sub_H_MARKERS = sub_H_MARKERS.loc[:, flag_valid_relPOS]

        print("\nFinally subsetted HLA dictionary file(\"maptable\").\n")
        print(sub_H_MARKERS.head())

        # sub_H_MARKERS.to_csv(_out+".3rd.txt", sep='\t', header=True, index=True)

        # (deprecated) - 2018. 8. 20.
        # # idx_rel_pos = pd.Index([ '_'.join(["AA", _hla_name, str(item[1])]) for item in sub_H_MARKERS.columns.tolist()])
        # idx_rel_pos = pd.Index([str(item[1]) for item in sub_H_MARKERS.columns.tolist()])
        # # (2018. 6. 19) 헤더에다가 그냥 relative position만 박음. ("AA_DRB1_-24" => "-24")
        #
        #
        # # # (2018. 6. 14) polymorphic position test
        # # sub_H_MARKERS.to_csv("./checkpolymorphic.txt", sep='\t', header=True, index=True)


        ### Reindexing `ASSOC_LOGISTIC` DataFrame with "relative_position".

        p_relPOS = re.compile("AA_{0}_(-?\d+)".format(_hla_name))

        ASSOC_LOGISTIC.index = pd.Index(ASSOC_MARKER_SET.str.extract(p_relPOS, expand=False).tolist())

        # This re-indexing will be used in next code block.




    ########## < Main job to process maptable and Making new association > ##########




    if MAKING_NEW_ASSOC:

        ###### < Making new *.assoc file > #####

        print(std_MAIN_PROCESS_NAME + "Making new *.assoc file.")

        # The main process will be done iterating the column of `sub_H_MARKERS` DataFrame.
        COLNAMES = sub_H_MARKERS.columns.tolist() # cf) len(COLNAMES) == len(col of `sub_H_MARKERS`)

        # # Indexing `ASSOC_MARKER_SET` with "relative position" beforehand.
        # p_relPOS = re.compile("AA_{0}_(-?\d+)".format(_hla_name))
        # ASSOC_MARKER_SET.index = pd.Index(ASSOC_MARKER_SET.str.extract(p_relPOS).tolist())


        for_new_assoc = []

        for i in range(0, sub_H_MARKERS.shape[1]):
        # for i in range(0, 5):

            print("\n=================================\n%d iteration\n=================================" % i)

            ## (variable in R) - (1) : `AAs`
            AAs = sub_H_MARKERS.iloc[:, i]
            print("\nAAs")
            print(AAs)


            ## (variable in R) - (2) : `AAname`
            AAname = COLNAMES[i] # ex) ('32557506', '-25', 'AAG')
            print("\nAAname")
            print(AAname)


            """
            (2018. 8. 20.)
            Remember that there could be more than 1 marker in one relative position.
            That's why one more searching job should be done to `ASSOC_LOGISTIC` DataFrame.
            
            ex) HLA-A, rel:-15
            1. AA_A_-15_29910359_L
            2. AA_A_-15_29910359_V
            3. AA_A_-15_29910359_x

            => Three markers in -15 relative position of Amino Acids.
            
            """

            # (2018. 8. 20.) This temporary DataFrame `df_temp` will do a central role in each iteration.
            df_temp = ASSOC_LOGISTIC.loc[[AAname[1]]]
            print("\ndf_temp : \n")
            print(df_temp)


            ## (variable in R) - (3) : `AAvariants`
            AAvariants = df_temp.loc[:, "SNP"]
            print("\nAAvariants")
            print(AAvariants)



            if len(AAvariants) > 0:

                """
                Introducing this condition (len(AAvariants) > 0) might seem trivial, but it filters unexpected markers 
                which don't appear in "*.assoc.logistic(`ASSOC_LOGISTIC`)" but in "HLA dictionary(`sub_H_MARKERS`).
                
                cf) len(AAvariants) == df_temp.shape[0] == The number of rows of `df_temp`
                
                """

                """
                As iterating over AA poisition of "HLA dictionary(`sub_H_MARKERS`)" and checking how many markers of 
                "*.assoc.logistc(`ASSOC_LOGISTIC`)" there are, the process should do it considering next two major cases.
                
                (1) Bi-allelic,
                (2) more than Tri-allelic.
                
                If the markers of `ASSOC_LOGISTIC` at each relative position is given as `df_temp`, then it would looks
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


                # AAvar_assoc (P-value column)
                AAvar_assoc = df_temp.loc[:, "P"]
                print("\nAAvar_assoc")
                print(AAvar_assoc)


                # Preparing index

                if len(AAvariants) > 1:

                    ### Tri-allelic ###

                    # The case where the position has more than 2 kind of markers (=> more than Tri-allelic)

                    index_AAvar_assoc = pd.Index([item.split('_')[-1] for item in AAvariants.tolist()])
                    # ex) ['AA_DRB1_-25_32665484_K', 'AA_DRB1_-25_32665484_R', 'AA_DRB1_-25_32665484_x']
                    # => ['K', 'R', 'x']

                    """
                    So, in the previous `df_temp` of Bi or Tri -allelic case,

                    SNP
                    AA_DRB1_-25_32665484_K    0.3133
                    AA_DRB1_-25_32665484_R    0.3173
                    AA_DRB1_-25_32665484_x    0.9150
                    Name: P, dtype: float64

                    
                    I will make this to next form.
                    
                    K    0.3133
                    R    0.3173
                    x    0.9150
                    Name: P, dtype: float64
                    
                    
                    I want to use the `df_temp` like this with that index.
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

                elif len(AAvariants) == 1:

                    ### Bi-allelic ###

                    print("\n\n# of label is 1")
                    print(df_temp.loc[:, "SNP"])
                    print(df_temp.loc[:, "A1"])

                    # refA = A1.loc[index_AAvar_assoc[0]]
                    refA = df_temp.loc[:, "A1"].iat[0]

                    index_AAvar_assoc = pd.Index([refA])



                AAvar_assoc.index = index_AAvar_assoc

                print("\nNew AAvar_assoc")
                print(AAvar_assoc)

                """                
                Newly made `AAvar_assoc`.

                K    0.3133
                R    0.3173
                x    0.9150
                Name: P, dtype: float64
                
                """


                ### This part is core job to transform "AA Marker Character" to the -log10 values.

                try:
                    # Tri-allelic or more.
                    AAs2 = AAs.apply(lambda x : -log10(AAvar_assoc.loc[x]))
                except KeyError:
                    # Bi-allelic
                    AAs2 = AAs.apply(lambda x : -log10(AAvar_assoc.loc[refA]))

                print("\nAAs2")
                print(AAs2)



                ### Flippng based on OR

                t_OR = df_temp.loc[:, "OR"]
                t_OR.index = index_AAvar_assoc
                print("\nOR")
                print(t_OR)

                print("\nCondition")
                try:
                    AAs3 = AAs.apply(lambda x : (2*int(t_OR.loc[x] > 1) - 1))
                except KeyError:
                    AAs3 = AAs.apply(lambda x : (2*int(t_OR.loc[refA] > 1) - 1))
                    # (2018. 8. 21.) 여기 Bi-allelic한케이스 flipping이 덜 됐음.

                print(AAs3)

                AAs4 = AAs2*AAs3
                print("\nAAs4")
                print(AAs4)

                for_new_assoc.append(AAs4)


        NEW_ASSOC = pd.DataFrame(for_new_assoc)

        # idx_NEW_ASSOC = [ '_'.join(["AA", _hla_name, item[1]]) for item in NEW_ASSOC.index.tolist()]
        idx_NEW_ASSOC = [str(item[1]) for item in NEW_ASSOC.index.tolist()] # (?) check it later whether it is needed or not.
        # NEW_ASSOC.index = idx_NEW_ASSOC

        print("\nNEW_ASSOC : \n{0}\n\n".format(NEW_ASSOC.head()))


    if MAKING_ASSOC_P:

        ##### < Making Assoc_P file > #####

        print(std_MAIN_PROCESS_NAME + "Making Assoc_P file.\n\n")

        print("\nmaptable")
        print(sub_H_MARKERS.head())
        print(sub_H_MARKERS.index.tolist())
        print(len(sub_H_MARKERS))


        ### Preparing "HLA_DRB1_***" marker from "*.assoc.logistic" again.

        """
        앞에서도 "*.assoc.logistic" 파일을 로드해서 다뤘지만 구체적으로 앞에서는 "AA_DRB1_****" 이런형태의 marker들이었음.
        이번에도 똑같이 logistic regression의 결과가 필요하지만 이번에는 "HLA_DRB1_****"이런 형태의 marker들을 활용해야함.
        
        한 가지 recurrent하는 문제는 2,4-field Input/Output 문제임.
        Output의 경우 뭐 4-field로 작업해서 기본적으로 4-field로 만들어 뽑다가 필요하면 2-field로 잘라서 출력해주는건 쉬우니까 큰 문제는 없고.
        
        문제는 Input으로 받을때 
        
        (1) 2-field일 경우
            => "HLA_DRB1_\d{4}"
            (ex) "HLA_DRB1_0101"
            
        (2) 4-field일 경우
            => "HLA_DRB1_\d{2}\:"
            (ex) "HLA_DRB1_01:01:01:01}
            
        이렇게 받아올거임.
        
        또, 이렇게 새로 받아온 `assoc`에서 "P" 컬럼만 따로 때서 `assocsrt`를 준비해야함. 
        
        
        사실 생각해보니, 그래서 "HLA_DRB1_0101"이렇게 생긴애들을 모아다가
        
        저 레이블에서 0101추리고,
        0101 TEMPORARY_TRANSFORMATION에 집어넣어서 digit checking해서 DRB1*01:01 이렇게 만들고,
        DRB1*01:01이거를 TWOtoFOUR dictionary에 집어넣어서 DRB1*01:01:01 이렇게 만들면되는데.
        
        사실 생각해보면 방금까지 언급한 step들이 위에서 `filtered_HLA_MARKER_DICTIONARY` (진짜 딱 `filtered_HLA_MARKER_DICTIONARY`까지만 만드는 과정)
        까지 만드는 과정과 똑같음. 결국 `TWOtoFOUR` 딕셔너리에 집어넣어서 KeyError뜨는 레이블들까지만 필터링 해놓은 결과물.
        (cf. `filtered_HLA_MARKER_DICTIONARY2`는 relative position당 marker의 종류가 2개 이상 나타난 애들, 
             `sub_H_MARKERS`는 "*.assoc.logistic"의 marker set과의 교집합.)
        
        그냥 얘를 가져다 쓰는것도 괜춘할듯.
        
        
        야야야 지금 "HLA_DRB1_****"이런 마커 로드하는 작업까지 위에서 한번씩 다 했네 여기 싹다 갈아 엎어도 되겠다.
        멍청아 그건 maptable에 대해서 한거고 지금은 "*.assoc.logistic"에 대해서 "HLA_DRB1_****" 마커 캐는거라고
        
        """

        if _oldv:

            assoc = assocfile.filter(regex="_".join(["HLA", _hla_name, "\d{4,}"]), axis=0)
            # 앞서서는 "*.assoc.logistic"에서 "AA_DRB1_****"형태인 marker를 활용했고, 이번에는 다시 "HLA_DRB1_****"형태의 marker를 활용할 차례.

            print("\nfirst assoc")
            print(assoc.head(30))

            # _oldv 라는 가정하에 또 '0101' => '01:01" 이렇게 바꿔줘야함.
            # 그때 함수만들어놓은거 쓰면 될거같음.

            idx_assoc = assoc.index.tolist()

            p = re.compile("\d{4}")

            idx_assoc = pd.Series(idx_assoc).apply(lambda x : p.search(string=x).group()).apply(lambda x : single_2DIGIT_CHECK(_hla_name, x, H_MARKERS))
            print("\nold idx_assoc")
            print(idx_assoc)
            # 여기까지 "*.assoc.logistic" 파일의 "HLA_DRB1_****" 마커들을 추리고, 4-digit추리고, 이거 DRB1*01:01형태로 까지 만들었음
            # 여기까지 작업한 idx_assoc는 아직 "DRB1*01:60"같은 쭉정이 allele들을 포함하고 있음. 여기서 idx_assoc의 원소들 하나하나를 `TWOtoFOUR`에 집어넣어서 쭉정이 allele들만 걸러내면

            No_Dummy_allele = idx_assoc.apply(lambda x : True if x in TWOtoFOUR.keys() else False).tolist()
            print("\nNo_Dummy")
            print(No_Dummy_allele)

            idx_assoc = idx_assoc.loc[No_Dummy_allele]
            assoc = assoc.loc[No_Dummy_allele]

            print("\nNo JJookJungYee")
            print(idx_assoc.head())

            idx_assoc = idx_assoc.apply(lambda x : TWOtoFOUR[x]) # Now, Finally 4-field allele names for "HLA_DRB1_****" markers for "*.assoc.logistic" 파일
            print("\nnew idx_assoc")
            print(idx_assoc)
            print(len(idx_assoc))

            assoc.index = pd.Index(idx_assoc)

        else:
            assoc = assocfile.filter(regex="_".join(["HLA", _hla_name, "\d{2,}\:"]), axis=0)


            # (2018. 6. 18) 여기 input이 4-field일때 대응해서 처리하는 코드 나중에 추가할 것.



        print("\nold assoc")
        print(assoc)

        assocrst = assoc.loc[:, "P"]
        OR = assoc.loc[:, "OR"]
        print("\nassocrst")
        print(assocrst)


        # `filtered_HLA_MARKER_DICTIONARY`의 marker_label로 "*.assoc.logsitc"의 일부인 `asscrst`의 값에 접근하는 상황임.
        alleleP = [-log10(assocrst.loc[item]) for item in sub_H_MARKERS.index.tolist()]
        print("\nalleleP")
        print(alleleP)

        ### Flipping
        alleleP = [alleleP[i]*(2*(OR.iat[i]>1)-1) for i in range(0, len(alleleP))]
        print("\nflipped alleleP")
        print(alleleP)

        alleleP = pd.DataFrame(alleleP)
        alleleP.index = pd.Index(assoc.loc[:, "SNP"])

        print("\nFinal allelP")
        alleleP.index.name = None
        # alleleP.columns.name = "assocrst"
        print(alleleP)

        ##### OUTPUT 3. alleleP.txt file


    if EXPORTING_OUTPUT:


        ### Introducing new indexes(maptable, assoc)
        NEW_ASSOC = NEW_ASSOC.transpose()
        NEW_ASSOC.columns = idx_NEW_ASSOC
        # (2018. 6. 19) 에전에 작업하던거 왜 헤더가 n-1인지 알겠네. 첫번째 라인만 n-1개인데 헤더로 잡아주면 그 아랫줄 부터는 첫번째 컬럼은 모두 인덱스로 알아서 잡음. 다음단계에서 편하게 작업하려고 전부 첫번째 라인만 ncol-1 로 일부러 잡으신거같음 교수님께서.

        sub_H_MARKERS.columns = idx_rel_pos
        sub_H_MARKERS = sub_H_MARKERS.loc[:, NEW_ASSOC.columns.tolist()]

        # 혹시 몰라서...
        if sub_H_MARKERS.shape[1] != NEW_ASSOC.shape[1]:
            print("\n[heatmap]: The number of relative position(# of markers) to plot heatmap is different!\nSomething Wrong!\n")
            sys.exit()

        """
        (2018. 6. 19)
        그 앞에서도 한번 언급했는데, maptable에서 잘 필터링해도 "*.assoc.logistic"의 결과물에는 안나타나는 AA relative position이 몇개 있음.
        대표적인 예로 -2인 포지션이 maptable에는 계속 남아있는데 "*.assoc.logistic"에는 안나타남.
        
        이 때문에 그냥 `filtered_HLA_MARKERS_DICTIONARY`를 파일로 방출해버리면 다음 "8b_plot.R"에서 문제가 생겨서 그냥 최종적으로 "*.assoc.txt"의 relative position
        에 해당하는 애들만 한번더 추려서 maptable을 수정하고 export시키기로 함. 왜냐하면 어차피 heatmap을 그리는 본질적인 목표가 "*.assoc.logistic"에 있는 마커들에 관심이
        있는 거기 때문에 이게 낫겠음. 어차피 maptable을 이 로지스틱 리그레션을 거치고 나온 marker들을 그리기 위한 준비물들일 뿐임. 
        """

        sub_H_MARKERS.to_csv(_out+'.map.txt', sep='\t', header=True, index=True)
        alleleP.to_csv(_out+".alleleP.txt", sep='\t', header=True, index=True)
        NEW_ASSOC.to_csv(_out+".assoc.txt", sep='\t', header=True, index=True)

        # (2018. 6. 19.) 파일이름 관련해서 예전에는 _map.txt 에서 .map.txt이렇게 바꿈




    if PLOT_HEATMAP:

        print("\n[heatmap]: Plotting.\n\n")

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




    return 0


# Submodule for `HEATMAP()` function

def MAPPING_2to4FIELD(_hla_name, _HLA_MARKERS_DICTIONARY):

    """
    HLA Name (ex. "A", "C", "DRB1" 등)과 이에 해당하는 HLA_MARKERS_DICTIONARY를 가지고 2-field => 4-field mapping을 함.
    지금은 HLA_MARKERS_DICTIONARY를 파이썬 자료구조인 Dictionary로 받았는데, 추후에 확장할때 바로 파일로 받아오는것도 괜춘할듯.

    예상하다시피 4-field를 2-field로 다듬으면 비결정성이 발생함. 예를 들어, DRB1*01:01:01:01, DRB1*01:01:01:02, DRB1*01:01:02:01
    얘네들 모두 2-field로 자르면 똑같이 DRB1*01:01이 될거임. 이걸 역으로, 이렇게 만들어진 DRB1*01:01을 key로 해서 groupby를 해다가

        DRB1*01:01 => DRB1*01:01:01:01

    이렇게 매핑한거임(DRB1*01:01:01:01, DRB1*01:01:01:02, DRB1*01:01:02:01 이 세 개 중에 그냥 첫번째껄로. 나머지 2가지 경우의 수는 버려지는거지만,
    반대로 얘기해서 실제로 2-field까지만 주어진 정보로 저 3가지 중 어떤게 맞는지는 알 수 없음. 다시 한번 말하지만 이 framework는 2-field기반으로 작동하던 이전의
    framework들을 4-field기반과 호환되게 만드는 일환일 뿐임).

    최종적으로 crude mapping결과를 TWOtoFOUR 딕셔너리로 만들어 리턴하는 거임.

    """

    """
    (2018. 6. 11)
    As an extension of `COATING_PED.py`, other modules that uses output of MakeDictionary Framework need to have this sub-module.
    This module conduct temporary transformation
    
    On the second thought, it could be the core function of `COATING_PED.py`...
    Maybe `COATING_PED.py` could be re-written with this function(Seriously...)
    
    
    Anyway, this function will be moved to `COATING_PED.py` script.
    
    """

    if not isinstance(_HLA_MARKERS_DICTIONARY, pd.DataFrame):
        print("\n[heatmap/ApproximationBox]: Given _HLA_MARKERS_DICTIONARY is not a DataFrame. Please Check it again.\n\n")
        sys.exit()

    # No problem with given `_HLA_MARKERS_DICTIONARY`

    LABELS_4field = pd.DataFrame(_HLA_MARKERS_DICTIONARY.index.tolist())
    LABELS_4field.columns = ["Alleles4"]

    LABELS_4field = LABELS_4field.set_index("Alleles4", drop=False)

    # 다음의 정규표현식으로 4-field를 2-field로 잘라냄
    p = re.compile("\*".join([_hla_name, "\d{2}\:\d{2}"]))

    LABELS_2field = LABELS_4field.applymap(lambda x : p.match(string=x).group())
    LABELS_2field.columns = pd.Index(["Alleles2"])

    # (2018. 6. 10) Because filtering to include only "DRB1" not "DRB2,3,4" was done when loading `HLA_MARKERS_DICTIONARY` file.
    # No need to give if-else in lambda expression(I mean, no need to code like below.)
    #
    # lambda x : p.match(string=x).group() if p.match(string=x) else ""


    # print(LABELS_2field.head())


    LABELS = LABELS_2field.reset_index("Alleles4")
    groupby_2field = LABELS.groupby(["Alleles2"])


    TWOtoFOUR = {k1: group.iat[0,0] for k1, group in groupby_2field}

    print(TWOtoFOUR)



    return TWOtoFOUR






def single_2DIGIT_CHECK(_hla_name, _the_allele, _HLA_MARKERS_DICTIONARY):

    t_name = _the_allele

    if len(t_name) == 4:
        # print(t_name)
        # OUTPUT.append(''.join([_hla_name, "*", t_name[0:2], ":", t_name[2:4]]))
        return ''.join([_hla_name, "*", t_name[0:2], ":", t_name[2:4]])

    elif len(t_name) == 2:
        # print(t_name)
        # OUTPUT.append(''.join([_hla_name, "*", t_name[0:2]]))
        return ''.join([_hla_name, "*", t_name[0:2]])

    elif len(t_name) == 5:

        # print(t_name)

        try:
            t_name2 = ''.join([_hla_name, "*", t_name[0:2], ":", t_name[2:5]])
            _HLA_MARKERS_DICTIONARY.loc[t_name2, :]  # If 2:3 digit name is wrong, then the program will split an error here.

        except KeyError:
            # If error is spitted out above line, then choosing 3:2 digit name is enough.
            t_name2 = ''.join([_hla_name, "*", t_name[0:3], ":", t_name[3:5]])


        return t_name2



if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################
    
        heatmap.py
    
        - Originally, Plotting heatmap entailed too many steps based on bash, python and R. Now it is
        integrated to this python script. It will be used in "HLA_Analysis.py" script.
    
    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-hla", help="\nHLA gene name.\n\n", required=True,
                        choices = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'])
    parser.add_argument("-o", help="\nOutput file name prefix.\n\n", required=True)


    # for Merged Prefix use(".bim", ".assoc.logistc", ...).
    parser.add_argument("--input", "-i", help="\nPrefix of mereged input files(\".bed\", \".covar\", or \".phe\" etc).\n\n")

    # Respective .bim files.
    parser.add_argument("--bim-HLA", "-bh", help="\n\".bim\" file of HLA markers.\n\n", default="Not_given")

    # (2018. 8. 20.) Deprecate next two arguments later.
    parser.add_argument("--bim-AA", "-ba", help="\n\".bim\" file of AA markers.\n\n", default="Not_given")
    # Merged .bim file.
    parser.add_argument("--bim-merged", "-bm", help="\nMerged \".bim\" file.\n\n", default="Not_given")

    # Output from MakeDictionary.
    parser.add_argument("--hla-dict", "-hd", help="\nMarker Dictionary file generated by 'MakeDictionary'.\n"                                                    "(If not given, it will search '~/data/HLA-Analysis/plot/heatmap' in default.)\n\n")
    # Result of Association Test.
    parser.add_argument("--logistic-result", "-lr", help="\nOutput from logtistic regression(*.assoc.logistic) by plink.\n\n")


    # We might need "*.alleles.DRB1" but not for now. because it could be no longer needed in generalizing 2,4-field HLA allele naming system.





    ##### < for Test > #####

    # (2018. 8. 20.)
    # args = parser.parse_args(["-hla", "DRB1",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_in_HATK2",
    #                           "-hd", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/WRAPER_TEST_DRB1.AA.markers.trim.labeled.txt",
    #                           "-lr", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/forHEATMAP_test.assoc.logistic",
    #                           "-bm", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/merged"
    #                           ])

    # Cancer data example.
    args = parser.parse_args(["-hla", "A",
                              "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer",
                              "-hd", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/heatmap/forHATK_A.AA.markers.trim.labeled.txt",
                              "-lr", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.AA.CODED.assoc.logistic",
                              "-bh", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.HLA.bim",
                              "-ba", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/AssociationTest/CancerResearch_Example/MARKER_PANEL/data_Rev_merged.AA.CODED.bim"
                              ])


    ##### < for Publish > #####

    # args = parser.parse_args()
    print(args)
    # print(vars(args))


    ##### < Additional Argument Processing. > #####



    # main function execution.
    HEATMAP(_hla_name=args.hla, _out=args.o, _input=args.input,
            _p_hla_dict=args.hla_dict, _p_assoc_logistic=args.logistic_result,
            _bim_HLA=args.bim_HLA, _bim_AA=args.bim_AA, _bim_merged=args.bim_merged,
            _p_heatmapR="/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/src/HLA_Analysis/Plotting/heatmap/8b_plot_WS.R")
