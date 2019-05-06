# -*- coding: utf-8 -*-


import os, sys, re
import argparse, textwrap
import pandas as pd
from collections import OrderedDict

########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
header_ped = ["FID", "IID", "PID", "MID", "Sex", "Phe"]


def NomenCleaner(_hped, _hped_descriptor, _iat, _imgt, _out, _field_format, __f_NoCaption=False,
                 __leave_NotFound=False):

    """
    """

    """
    [Input]
    `_hped_descriptor` := 1 (standard 4-field allele),
                          2 (G-group allele),
                          3 (P-group allele)
                          
    [Output]
    `_field_format` := 
        1 (1-field)
        2 (2-field)
        3 (3-field)
        4 (4-field)
        5 (G-group)
    
    """



    ########## < Core Variables > ##########

    __HPED__ = None
    __IAT__ = None

    ### Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')

    if bool(os.path.dirname(_out)):
        INTERMEDIATE_PATH = os.path.dirname(_out)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        INTERMEDIATE_PATH = "./"




    ########## < [1] Loading "*.(h)ped" file > ##########

    __HPED__ = pd.read_table(_hped, sep='\t', header=None, dtype=str,
                        names=header_ped + [item + "_" + str(i) for item in HLA_names for i in range(1, 3)]).set_index((header_ped))

    # print(std_MAIN_PROCESS_NAME + "Loaded \"*.ped\" file.")
    # print(__HPED__.head())




    ########## < [2] Loading "*.iat" file > ##########

    __IAT__ = pd.read_table(_iat, sep='\t', header=0, dtype=str)
    # print(std_MAIN_PROCESS_NAME + "Loaded \"*.iat\" file.\n")
    # print(__IAT__.head())


    ## Dividing `__IAT__` by HLA

    df_AllelebyHLA = __IAT__.loc[:, "Allele"].str.split(r'\*', expand=True)
    df_AllelebyHLA.columns = ["HLA", "Allele"]

    __IAT__ = pd.concat([df_AllelebyHLA, __IAT__.loc[:, ["G_group", "P_group"]]], axis=1).set_index('HLA')
    # print(__IAT__.head())


    __IAT_dict__ = {HLA_names[i]: __IAT__.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}

    # print("\nIAT divided by HLA gene names.\n")
    # print(__IAT_dict__["C"].head(10))




    ########## < [3] Searching given HLA alleles in *.iat file. > ##########


    count = 0
    l_rows_Searched = []
    l_rows_Check = []

    for t_row in __HPED__.itertuples():

        """
        t_row[0] := index (ex. ('WTCCC125636', 'WTCCC125636', '0', '0', '0', '-9'))
        t_row[1] := A_1 (ex. '0101');   t_row[2] := A_2 (ex. '0101');
        t_row[3] := B_1 (ex. '4402');   t_row[4] := B_1 (ex. '5701')        
        ...
        
        """
        # print(t_row)

        l_t_HLAs = []
        l_t_Check = []

        for j in range(0, len(HLA_names)):

            t_al1 = t_row[2*j+1]
            t_al2 = t_row[2*j+2]

            ### Searching given alleles in *.iat file.

            # print("{} and {}".format(t_al1, t_al2))

            ## al1
            if t_al1 != '0':
                t_new_allele1 = SearchIAT(t_al1, HLA_names[j], __IAT_dict__[HLA_names[j]], _hped_descriptor)
            else:
                t_new_allele1 = t_al1   # Just let it be.


            ## al2
            if t_al2 != '0':
                t_new_allele2 = SearchIAT(t_al2, HLA_names[j], __IAT_dict__[HLA_names[j]], _hped_descriptor)
            else:
                t_new_allele2 = t_al2


            ### Searched HLA alleles
            l_t_HLAs.append(t_new_allele1)
            l_t_HLAs.append(t_new_allele2)


            ### Check whether Searched or not.

            ## al1
            if t_new_allele1 == "-1":
                # Not found in *.iat file
                l_t_Check.append(t_al1)
            else:
                l_t_Check.append('F') # 'F'ound in *.iat

            ## al2
            if t_new_allele2 == "-1":
                # Not found in *.iat file
                l_t_Check.append(t_al2)
            else:
                l_t_Check.append('F') # 'F'ound in *.iat


        # print(l_t_HLAs)
        l_rows_Searched.append(l_t_HLAs)
        l_rows_Check.append(l_t_Check)

        count += 1
        # if count > 5 : break


    ### DataFrame of Searched alleles.
    df_IAT_Searched = pd.DataFrame(l_rows_Searched)
    # print(df_IAT_Searched.head())
    # df_IAT_Searched.to_csv(_out+'.justSearched.chped', sep='\t', header=False, index=False)


    ### DataFrame of marks of searched alleles.
    df_Check_Searched = pd.DataFrame(l_rows_Check)
    # print(df_Check_Searched.head())



    ########## < [4] Main Transformation > ##########

    ### Preparing Hashtable of mapping rule

    __IAT_dict2__ = {}

    col_from = -1
    col_to = -1

    if _hped_descriptor == 1:

        ## Input hped : Standard 4-field
        col_from = 0

        if _field_format <= 4:
            # to 4-field
            pass   # No need to do transformation (4-field to 4-field)
        elif _field_format == 5:
            # to G-group
            col_to = 1
        elif _field_format == 6:
            # to P-group
            col_to = 2


    elif _hped_descriptor == 2:

        ## Input hped : G-group
        col_from = 1

        if _field_format <= 4:
            col_to = 0
        elif _field_format == 5:
            # to G-group
            pass   # No need to do transformation (G-group to G-group)
        elif _field_format == 6:
            # to P-group
            col_to = 2


    elif _hped_descriptor == 3:

        ## Input hped : P-group
        col_from = 2

        if _field_format <= 4:
            col_to = 0
        elif _field_format == 5:
            # to G-group
            col_to = 1
        elif _field_format == 6:
            # to P-group
            pass   # No need to do transformation (P-group to P-group)


    if col_from != -1 and col_to != -1:

        for i in range(0, len(HLA_names)):

            from_which = "Allele" if _hped_descriptor == 1 else "G_group" if _hped_descriptor == 2 else "P_group"

            df_temp = __IAT_dict__[HLA_names[i]].drop_duplicates(from_which)
            # print(df_temp.head())

            sr_from = df_temp.iloc[:, col_from]
            sr_to = df_temp.iloc[:, col_to]

            sr_to.index = sr_from

            __IAT_dict2__[HLA_names[i]] = sr_to.to_dict()

    # print("__IAT_dict2__")
    # print(__IAT_dict2__)


    ### Transforming HLA alleles searched in "*.iat" (`df_IAT_Searched`) to (1) standard 4-field, (2) G-group, or (2) P-group allele.

    df_Transformed = Main_Transformation(df_IAT_Searched, __IAT_dict2__, _f_NoCaption=__f_NoCaption)
    df_Transformed.index = __HPED__.index
    # print("Completely Transformed dataframe.")
    # print(df_Transformed.head())



    ########## < [5] Final Trimming > ##########

    ### Trimming out `df_Transformed`(standard 4-field) (1-field, 2-field, 3-field, or old-format(MakeReference)).

    if _field_format < 4:

        if _field_format == 0:
            ## Old-format (for compatibility with original version of MakeReference(by S. Jia))
            p_old_format = re.compile(r'^(\w+\*)?(\d{2,3}:){0,1}\d{2,3}[A-Z]?')

            df_Transformed = df_Transformed.applymap(lambda x: p_old_format.match(x).group() if bool(p_old_format.match(x)) else x).applymap(lambda x : re.sub(pattern=r'\*', string=x, repl=':'))


        elif _field_format == 1:
            ## 1-field
            p_1field = re.compile(r'^(\w+\*)?\d{2,3}')

            df_Transformed = df_Transformed.applymap(lambda x: p_1field.match(x).group() if bool(p_1field.match(x)) else x)


        elif _field_format == 2:
            ## 2-field
            p_2field = re.compile(r'^(\w+\*)?(\d{2,3}:){0,1}\d{2,3}[A-Z]?')

            df_Transformed = df_Transformed.applymap(lambda x : p_2field.match(x).group() if bool(p_2field.match(x)) else x)


        elif _field_format == 3:
            ## 3-field
            p_3field = re.compile(r'^(\w+\*)?(\d{2,3}:){0,2}\d{2,3}[A-Z]?')

            df_Transformed = df_Transformed.applymap(lambda x : p_3field.match(x).group() if bool(p_3field.match(x)) else x)


        # print(df_Transformed.head())
        # df_Transformed.to_csv(_out+".tEst.{}.chped".format(MakeSuffixInfo(_imgt, _field_format)), sep='\t', header=False, index=True)




    ### Integrating Not found allele(novel allele) ... Yang requested.

    count = 0

    l_rows_merged = []

    for row in df_Check_Searched.itertuples():

        l_row = []

        for j in range(0, df_Transformed.shape[1]):

            idx = j + 1

            t_allele = None

            if row[idx] == 'F':
                t_allele = df_Transformed.iat[count, j]
            else:
                if __leave_NotFound:
                    t_allele = row[idx]
                else:
                    t_allele = "0"

            l_row.append(t_allele)

        l_rows_merged.append(l_row)

        count += 1

    __RETURN__ = pd.DataFrame(l_rows_merged, index=__HPED__.index)
    # print(__RETURN__.head())


    ### Writing final output(*.chped)
    __RETURN__.to_csv(_out+".{}.chped".format(MakeSuffixInfo(_imgt, _field_format)), sep='\t', header=False, index=True)



    return 0



def CHECK_DIGITS(_hla_name, _the_allele, _IAT_Allelelist):

    """

    Perform this function assuming given `_the_allele` neither is captioned or has double-colon.

    (2018. 6. 29.)
    More classification conditions are added.

    (2018. 7. 23.)
    I found the way to deal with a Tag has not been introduced yet.
    Now it is introduced.

    """

    ### Step1. Checking suffix character (some single character tagged along with the given allele).
    p_tag = re.compile(r'.+[A-Z]$')
    hasTag = p_tag.match(_the_allele)

    if hasTag:
        t_name = _the_allele[:-1]
        t_name_tag = _the_allele[-1]
    else:
        t_name = _the_allele
        t_name_tag = -1



    ### Step2. main digit checking part.

    if len(t_name) == 2:
        # 1-field / 2-digits
        return t_name[0:2] + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 3:
        # 1-field / 3-digits
        return t_name[0:3] + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 4:
        # 2 + 2 (+C) = 4
        return ':'.join([t_name[0:2], t_name[2:4]]) + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 5:
        # (1) 2 + 3 (+C) = 5
        # (2) 3 + 2 = 5

        p = re.compile('\:'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any = _IAT_Allelelist.str.match(p).any()

        return (':'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
            (':'.join([t_name[0:3], t_name[3:5]]) + (t_name_tag if t_name_tag != -1 else ""))

    elif len(t_name) == 6:
        # (1) 2 + 2 + 2 (+C)
        # (2) 3 + 3 (+C)

        p = re.compile('\:'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any = _IAT_Allelelist.str.match(p).any()

        return (':'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (
            t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
            (':'.join([t_name[0:3], t_name[3:6]]) + (t_name_tag if t_name_tag != -1 else ""))

    elif len(t_name) == 7:
        # (1) 2 + 3 + 2 (+C) = 7
        # (2) 2 + 2 + 3 (+C) = 7
        # (3) 3 + 2 + 2 (+C) = 7

        p1 = re.compile('\:'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any1 = _IAT_Allelelist.str.match(p1).any()

        if not Found_Any1:
            # Not the case of "(1) 2 + 3 + 2 = 7"
            p2 = re.compile(
                '\:'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any2 = _IAT_Allelelist.str.match(p2).any()

            if not Found_Any2:
                # Not the case of "(2) 2 + 2 + 3 = 7"
                p3 = re.compile(
                    '\:'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
                Found_Any3 = _IAT_Allelelist.str.match(p3).any()

                if not Found_Any3:
                    # Not the case of "(3) 3 + 2 + 2 = 7"
                    # Not Found
                    return "-1"
                else:
                    # The case of "(3) 3 + 2 + 2 = 7"
                    return ':'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")
            else:
                # The case of "(2) 2 + 2 + 3 = 7"
                return ':'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else "")
        else:
            # The case of "(1) 2 + 3 + 2 = 7"
            return ':'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")

        # return ':'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) if Found_Any else \
        #     ':'.join([t_name[0:3], t_name[3:5], t_name[5:7]])

    elif len(t_name) == 8:
        # (1) 2 + 2 + 2 + 2
        # (2) 3 + 2 + 3
        # (3) 3 + 3 + 2

        p1 = re.compile(
            '\:'.join([t_name[0:2], t_name[2:4], t_name[4:6], t_name[6:8]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any1 = _IAT_Allelelist.str.match(p1).any()

        if not Found_Any1:
            # Not the case of "(1) 2 + 2 + 2 + 2"
            p2 = re.compile(
                '\:'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any2 = _IAT_Allelelist.str.match(p2).any()

            if not Found_Any2:
                # Not the case of "(2) 3 + 2 + 3"
                p3 = re.compile(
                    '\:'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else ""))
                Found_Any3 = _IAT_Allelelist.str.match(p3).any()

                if not Found_Any3:
                    return str(-1)

                else:
                    return ':'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else "")
            else:
                return ':'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else "")
        else:
            return ':'.join([t_name[0:2], t_name[2:4], t_name[4:6], t_name[6:8]]) + (
                t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 9:
        # (1) 2 + 3 + 2 + 2
        # (2) 3 + 2 + 2 + 2

        p = re.compile(
            '\:'.join([t_name[0:2], t_name[2:5], t_name[5:7], t_name[7:9]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any = _IAT_Allelelist.str.match(p).any()

        return (':'.join([t_name[0:2], t_name[2:5], t_name[5:7], t_name[7:9]]) + (
            t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
            (':'.join([t_name[0:3], t_name[3:5], t_name[5:7], t_name[7:9]]) + (t_name_tag if t_name_tag != -1 else ""))

    elif len(t_name) == 10:
        # (1) 3 + 3 + 2 + 2

        return ':'.join([t_name[0:3], t_name[3:6], t_name[6:8], t_name[8:10]]) + (
            t_name_tag if t_name_tag != -1 else "")


def CHECK_DIGITS_PorGgroup(_hla_name, _the_allele, _IAT_Allelelist, _ped_descriptor):
    """

    ### < G-group > ###

    # 2 + 2 (+C) = 4
    - 46:51Q

    # 2 + 3 (+C) = 5
    - 08:124, 07:534
    - 01:247N, 01:250N

    # 3 + 2 = 5
    - 255:01, 527:01, 632:01, 137:01, 338:01

    # 2 + 2 + 2 (+C) = 6
    - 01:01:02
    - 01:01:01G

    # 2 + 3 + 2 = 7
    - 01:146:01
    - 06:127:01G

    # 2 + 2 + 3 = 7
    - 02:01:100, 02:01:101

    # 3 + 2 + 2 = 7
    - 155:01:01


    ### < P-group > ###

    # 2 + 2 = 4
    - 10:24, 11:04
    - 01:01P, 02:01P

    # 2 + 3 = 5
    - 14:110, 14:139
    - 01:146P, 02:101P, 02:610P

    # 3 + 2 = 5
    - 100:01, 119:01
    - 155:01P, 279:01P

    # 2 + 2 + 2 = 6
    - 30:12:01 (2018. 7. 5) I found this case lol!


    """


    # check whether some single character is tagged along with the given allele.
    p_tag = re.compile(r'.+[A-Z]$')

    hasTag = p_tag.match(_the_allele)

    if hasTag:
        t_name = _the_allele[:-1]
        t_name_tag = _the_allele[-1]
    else:
        t_name = _the_allele
        t_name_tag = -1

    # print("\n[CHECK_DIGITS_PorG] : {0} and {1}".format(t_name, t_name_tag))

    ### In case of "G-group"
    if _ped_descriptor == 2:

        # Main matching job

        if len(t_name) == 4:
            # (1) 2 + 2 (+C) = 4
            return ':'.join([t_name[0:2], t_name[2:4]]) + (t_name_tag if t_name_tag != -1 else "")

        elif len(t_name) == 5:
            # (1) 2 + 3 (+C) = 5
            # (2) 3 + 2 = 5

            p = re.compile('\:'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any = _IAT_Allelelist.str.match(p).any()

            return (':'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
                (':'.join([t_name[0:3], t_name[3:5]]) + (t_name_tag if t_name_tag != -1 else ""))

        elif len(t_name) == 6:
            # (1) 2 + 2 + 2 (+C) = 6

            return ':'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else "")

        elif len(t_name) == 7:
            # (1) 2 + 3 + 2 = 7
            # (2) 2 + 2 + 3 = 7
            # (3) 3 + 2 + 2 = 7

            p1 = re.compile(
                '\:'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any1 = _IAT_Allelelist.str.match(p1).any()

            if not Found_Any1:
                # Not the case of "(1) 2 + 3 + 2 = 7"
                p2 = re.compile(
                    '\:'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else ""))
                Found_Any2 = _IAT_Allelelist.str.match(p2).any()

                if not Found_Any2:
                    # Not the case of "(2) 2 + 2 + 3 = 7"
                    p3 = re.compile(
                        '\:'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
                    Found_Any3 = _IAT_Allelelist.str.match(p3).any()

                    if not Found_Any3:
                        return "-1"

                    else:
                        return ':'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (
                            t_name_tag if t_name_tag != -1 else "")
                else:
                    return ':'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else "")
            else:
                return ':'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")

        else:
            # print("\nNo match for G_group.\n")
            return "-1"


    ### In case of "P-group"
    elif _ped_descriptor == 3:

        # Main Matching Job

        if len(t_name) == 4:
            # 2 + 2 = 4
            return ':'.join([t_name[0:2], t_name[2:4]]) + (t_name_tag if t_name_tag != -1 else "")

        elif len(t_name) == 5:
            # 2 + 3 = 5
            # 3 + 2 = 5

            p = re.compile('\:'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any = _IAT_Allelelist.str.match(p).any()

            return (':'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else (
                    ':'.join([t_name[0:3], t_name[3:5]]) + (t_name_tag if t_name_tag != -1 else ""))

        elif len(t_name) == 6:
            # 2 + 2 + 2
            # 3 + 3 ... Not Found.

            return ':'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else "")


        else:
            # print("\nNo match for P_group.\n")
            return "-1"


def Find_1st_Allele(_the_allele, _IAT_Allelelist):

    """
    In case where given allele has double-colon(implying its digits are determined) but not in a form of complete allele name in "Allelelist.txt" file,
    then i will classify this allele as "Non-deterministic 4-field allele".

    So, that incomplete allele("Non-deterministic 4-field allele") will be matched by "re.match" function and return 1st element.
    """


    if 2 <= len(_the_allele) <= 5:

        """
        (2018. 7. 23.)
        This classification block was introduced as I found some unproper exceptions.
        For example, When the allele "B*1501" is given, the transformed result is "15:170" not "B*15:17:01:01".
        It means digit information isn't fully considered in this function, i.e. this part which makes regular expression pattern.

        As a solution, I decided to add more classification related to `len(_the_allele)`.
        """

        Flag_Found = _IAT_Allelelist.str.match(str(_the_allele) + ":")

        if Flag_Found.any():
            return _IAT_Allelelist.loc[Flag_Found].iat[0]
        else:
            # One more time
            Flag_Found = _IAT_Allelelist.str.match(str(_the_allele))

            if Flag_Found.any():
                return _IAT_Allelelist.loc[Flag_Found].iat[0]
            else:
                return 0


    else:
        p = re.compile(str(_the_allele))
        # Flag_Found = _IAT_Allelelist.apply(lambda x: bool(p.match(string=x)))
        Flag_Found = _IAT_Allelelist.str.match(str(_the_allele))

        if Flag_Found.any():
            return _IAT_Allelelist.loc[Flag_Found].iat[0]
        else:
            return 0


def Find_1st_Allele_PorGgroup(_the_allele, _IAT_Allelelist, _ped_descriptor):


    # in case the allele is not found in the table, return its original value
    Flag_Found = _IAT_Allelelist.str.match(_the_allele)

    if Flag_Found.any():
        return _IAT_Allelelist.loc[Flag_Found].iat[0]
    else:
        return 0

    # # (2018. 7. 4.) Found result `Found_df` could be 'str' object. So, you shouldn't use ".iat" function.
    #
    # # return str(Found_df.iat[0]) if len(Found_df) > 0 else "0"
    #
    # if isinstance(Found_df, pd.Series):
    #     # ex. "01:01:01G" => ["01:01:01:01", "01:01:01:02N", "01:01:01:03", ..., "01:253"]
    #     # Multiple results retured as "Series" object.
    #     return str(Found_df.iat[0])
    # elif isinstance(Found_df, str):
    #     return Found_df
    # else:
    #     return "0"


# def Main_Transformation(_single_allele, _hla_name, _df_IAT_Allelelist, _ped_descriptor):
#     """
#     Main Processes which were located in main for loop. They are now moved to this function.
#     From now on, those main processes will be preformed to single allele.
#     """
#
#     _IAT_Allelelist = _df_IAT_Allelelist.index.tolist()
#
#     ### [1] Checking whether given allele is captioned or not.
#
#     _al_filetered1 = isCaptioned(_single_allele, _hla_name)
#
#     ### [2] Checking whether given allele has double-colon or not.
#
#     _al_filetered2 = ""
#
#     if _ped_descriptor == 1:
#         # When given ped file is standard 4-field allele(where `_ped_descriptor` == 1).
#         _al_filetered2 = CHECK_DIGITS(_hla_name, _al_filetered1, _IAT_Allelelist) if not hasDoubleColon(
#             _al_filetered1) else _al_filetered1
#
#     elif _ped_descriptor == 2:
#         # When given ped file is P or G-group allele(where `_ped_descriptor` == 2 or 3).
#         _al_filetered2 = CHECK_DIGITS_PorGgroup(_hla_name, _al_filetered1, _df_IAT_Allelelist.loc[:, "G_group"],
#                                                 _ped_descriptor) if not hasDoubleColon(
#             _al_filetered1) else _al_filetered1
#
#     elif _ped_descriptor == 3:
#         _al_filetered2 = CHECK_DIGITS_PorGgroup(_hla_name, _al_filetered1, _df_IAT_Allelelist.loc[:, "P_group"],
#                                                 _ped_descriptor) if not hasDoubleColon(
#             _al_filetered1) else _al_filetered1
#
#     ### [3] Digit-Checking to determine proper number of digits per field.
#
#     _al_filetered3 = ""
#
#     if _ped_descriptor < 2:
#         _al_filetered3 = Find_1st_Allele(_al_filetered2, _IAT_Allelelist)
#     elif _ped_descriptor >= 2:
#         _al_filetered3 = Find_1st_Allele_PorGgroup(_al_filetered2, _df_IAT_Allelelist, _ped_descriptor)
#
#     return str(_al_filetered3)



def SearchIAT(_single_allele, _hla, _IAT, _hped_descriptor):

    t_allele = _single_allele

    ### step 1. Removing HLA gene caption part (ex. "A*", "DRB1*")

    p_GeneCaption = re.compile(r'^\w+\*')

    if p_GeneCaption.match(t_allele):
        t_allele = p_GeneCaption.sub(string=t_allele, repl='')



    ### step 2. Check Digits

    if not bool(re.search(r':', t_allele)):

        if _hped_descriptor == 1:
            t_allele = CHECK_DIGITS(_hla, t_allele, _IAT.loc[:, "Allele"])
        elif _hped_descriptor == 2:
            t_allele = CHECK_DIGITS_PorGgroup(_hla, t_allele, _IAT.loc[:, "G_group"], _hped_descriptor)
        elif _hped_descriptor == 3:
            t_allele = CHECK_DIGITS_PorGgroup(_hla, t_allele, _IAT.loc[:, "P_group"], _hped_descriptor)

    # print("Digit checked allele : {}".format(t_allele))


    ### step 3. Search IAT (Main job).

    __RETURN_allele__ = ""

    if _hped_descriptor == 1:
        __RETURN_allele__ = Find_1st_Allele(t_allele, _IAT.loc[:, "Allele"])
    elif _hped_descriptor == 2:
        __RETURN_allele__ = Find_1st_Allele_PorGgroup(t_allele, _IAT.loc[:, "G_group"], _hped_descriptor)
    elif _hped_descriptor == 3:
        __RETURN_allele__ = Find_1st_Allele_PorGgroup(t_allele, _IAT.loc[:, "P_group"], _hped_descriptor)



    ### Return

    if __RETURN_allele__ == 0:
        # Nothing found in *.iat file. (maybe novel HLA allele.)
        # Then, just return original value (requested by Yang.)
        return "-1"
    else:
        return __RETURN_allele__



def Main_Transformation(_df_IAT_Searched, _IAT_dict2, _f_NoCaption=False):


    _df_Transformed = None

    if bool(_IAT_dict2):

        l_temp = []

        for i in range(0, len(HLA_names)):

            idx1 = 2*i
            idx2 = idx1 + 1

            t_df_IAT_Searched = \
                _df_IAT_Searched.iloc[:, [idx1, idx2]].applymap(lambda x : _IAT_dict2[HLA_names[i]][x] if (x != '-1' and x != '0') else x)

            l_temp.append(t_df_IAT_Searched)
            # print(t_df_IAT_Searched)

        _df_Transformed = pd.concat(l_temp, axis=1)

    else:
        _df_Transformed = _df_IAT_Searched


    ### Check gene caption.
    if not _f_NoCaption:
        _df_Transformed = pd.concat([_df_Transformed.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if (x != '-1' and x != '0') else x) for i in range(0, len(HLA_names))], axis=1)


    return _df_Transformed



# def Main_Transformation2(_hped):
#     #### [Warning] Temporary Hard coding. This funciton will be either adopted to the main transformation function or deprecated.
#
#     with open(_hped, 'r') as f_hped:
#         count = 0
#
#         for l in f_hped:
#             t_line = re.split(r'\s+', l.rstrip('\n'))
#
#             """
#             [0,1,2,3,4,5] := ped file information
#             [6,7] := HLA-A,
#             [8,9] := HLA-B,
#             ...,
#             [20, 21] := HLA-DRB1
#             """
#             __ped_info__ = '\t'.join(t_line[:6])
#             __genomic_info__ = '\t'.join(
#                 [Old_Transformation(t_line[2 * i + 6], t_line[2 * i + 7], HLA_names[i]) for i in
#                  range(0, len(HLA_names))])
#
#             yield '\t'.join([__ped_info__, __genomic_info__]) + "\n"
#
#             count += 1
#             # if count > 5: break


# def Old_Transformation(_HLA_allele1, _HLA_allele2, _hla):
#     if (_HLA_allele1 == "0" and _HLA_allele2 == "0"):
#         return '\t'.join(["0", "0"])
#     else:
#
#         p = re.compile(r'\d{4}[A-Z]?$')  # Only 4-digit with or without single suffix is accepted in old transformation.
#         p_1field = re.compile(r'\d{2}')  # decided to also process 1-field HLA alleles.
#
#         # Allele1
#         if p.match(_HLA_allele1):
#             _new_al1 = ':'.join([_hla, _HLA_allele1[:2], _HLA_allele1[2:4]])
#         elif p_1field.match(_HLA_allele1):
#             _new_al1 = ':'.join([_hla, _HLA_allele1])
#         else:
#             _new_al1 = "0" if _HLA_allele1 == "0" else "-1"
#
#         # Allele2
#         if p.match(_HLA_allele2):
#             _new_al2 = ':'.join([_hla, _HLA_allele2[:2], _HLA_allele2[2:4]])
#         elif p_1field.match(_HLA_allele2):
#             _new_al2 = ':'.join([_hla, _HLA_allele2])
#         else:
#             _new_al2 = "0" if _HLA_allele2 == "0" else "-1"
#
#         return '\t'.join([_new_al1, _new_al2])



def MakeSuffixInfo(_imgt, _field_format):

    s_field = ""

    if _field_format == 1:
        s_field = "1field"
    elif _field_format == 2:
        s_field = "2field"
    elif _field_format == 3:
        s_field = "3field"
    elif _field_format == 4:
        s_field = "4field"
    elif _field_format == 5:
        s_field = "Ggroup"
    elif _field_format == 6:
        s_field = "Pgroup"
    elif _field_format == 0:
        s_field = "OLD"


    return "imgt{}.{}".format(_imgt, s_field)





if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        NomenCleaner.py
        
         
        


    #########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    # Input (1) : *.ped file
    PED_TYPE = parser.add_mutually_exclusive_group(required=True)
    PED_TYPE.add_argument("--hped", help="\nHLA Type Data with raw HLA allele(ex. 0101).\n\n", dest="hped")
    PED_TYPE.add_argument("--hped-Ggroup", help="\nHLA Type Data with raw G-group allele(ex. 020601G).\n\n", dest="hped_G")
    PED_TYPE.add_argument("--hped-Pgroup", help="\nHLA Type Data with raw P-group allele (ex. 0102P).\n\n", dest="hped_P")

    # Input (2) : *.iat file
    parser.add_argument("--iat", help="\nIntegrated Allele Table file(*.iat).\n\n", required=True)
    parser.add_argument("--imgt", help="\nSpecifying the IMGT-HLA version.\n\n", required=True)

    # Ouptut Prefix
    parser.add_argument("--out", "-o", help="\nOutput file prefix.\n\n", required=True)
    parser.add_argument("--leave-NotFound", help="\nLeaving HLA alleles which can't be found in given *.iat file(Novel or Erroneous allele) intact.\n\n", action='store_true')

    # Output format
    output_digit_selection = parser.add_mutually_exclusive_group(required=True)
    output_digit_selection.add_argument("--1field", help="\nOutput ped file as '1-field' format.\n\n",
                                        action="store_true", dest="oneF")
    output_digit_selection.add_argument("--2field", help="\nOutput ped file as '2-field' format.\n\n",
                                        action="store_true", dest="twoF")
    output_digit_selection.add_argument("--3field", help="\nOutput ped file as '3-field' format.\n\n",
                                        action="store_true", dest="threeF")
    output_digit_selection.add_argument("--4field",
                                        help="\nOutput ped file as '4-field(Current Standard Names)' format.\n\n",
                                        action="store_true", dest="fourF")
    output_digit_selection.add_argument("--G-group", help="\nOutput ped file as 'G-group' format.\n\n",
                                        action="store_true")
    output_digit_selection.add_argument("--P-group", help="\nOutput ped file as 'P-group' format.\n\n",
                                        action="store_true")

    output_digit_selection.add_argument("--old-format",
                                        help="\nOutput ped file as previous version format. (ex. A:01:01)\n\n",
                                        action="store_true")

    # Flag to remove HLA gene caption.
    parser.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*\").\n\n", action='store_true')




    ##### <for Test> #####

    ### < Raw 4-field *.hped >

    # # raw 4-field to standard 4-field / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.before.hped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.after",
    #                           "--4field",
    #                           "--leave-NotFound"
    #                           ])

    # # raw 4-field to standard 4-field / with Novel Allele / NotFound : 0
    # args = parser.parse_args(["-ped", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.before.hped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.after",
    #                           "--4field"
    #                           ])


    # # raw 4-field to G-group / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.before.hped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.after",
    #                           "--G-group",
    #                           "--leave-NotFound"
    #                           ])

    # # raw 4-field to P-group / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.before.hped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/wtccc_filtered_58C_NBS.raw4field.novelallele.after",
    #                           "--P-group",
    #                           "--leave-NotFound"
    #                           ])



    ### < raw G-group >

    # # raw G-group to 4-field / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--4field",
    #                           "--leave-NotFound"
    #                           ])

    # # raw G-group to 3-field / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--3field",
    #                           "--leave-NotFound"
    #                           ])

    # # raw G-group to 2-field / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--2field",
    #                           "--leave-NotFound"
    #                           ])

    # # raw G-group to 1-field / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--1field",
    #                           "--leave-NotFound"
    #                           ])

    # # raw G-group to old-format / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--old-format",
    #                           "--leave-NotFound"
    #                           ])



    # # raw G-group to G-group / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--G-group",
    #                           "--leave-NotFound"
    #                           ])


    # # raw G-group to P-group / with Novel Allele / --leave-NotFound
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190123_hg18_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg18.imgt3320.iat",
    #                           "-imgt", "3320",
    #                           "-o", "/Users/wansun/Projects/20190222_JSH_GWAS/NomenCleanerTest/DummyPED.Ggroup.Ncap.Ndc.10.novelallele",
    #                           "--P-group",
    #                           "--leave-NotFound"
    #                           ])


    ### < raw P-group to >

    # 4-field
    # args = parser.parse_args(["--hped-Pgroup", "/Users/wansun/Git_Projects/MakeReference_v2/data/NomenCleaner/v2/BEFORE_DummyHPED.rawPgroup.Ncap.Ndc.10.novelallele.hped",
    #                           "--iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190303_hg19_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg19.imgt3320.iat",
    #                           "--imgt", "3320",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/data/NomenCleaner/v2/from_Pgroup_to/EMERGENCY_AFTER_DummyHPED.rawPgroup.Ncap.Ndc.10.novelallele.imgt3320.Pgroup.chped",
    #                           "--4field",
    #                           "--leave-NotFound"
    #                           ])

    # P-group
    # args = parser.parse_args(["--hped-Pgroup", "/Users/wansun/Git_Projects/MakeReference_v2/data/NomenCleaner/v2/BEFORE_DummyHPED.rawPgroup.Ncap.Ndc.10.novelallele.hped",
    #                           "--iat", "/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190303_hg19_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg19.imgt3320.iat",
    #                           "--imgt", "3320",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/data/NomenCleaner/v2/from_Pgroup_to/EMERGENCY_AFTER_DummyHPED.rawPgroup.Ncap.Ndc.10.novelallele.imgt3320.Pgroup.chped",
    #                           "--P-group",
    #                           "--leave-NotFound"
    #                           ])





    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)



    ### Additional Argument processing


    ## Output Format Flags

    if args.oneF:
        FIELD_FORMAT = 1
    elif args.twoF:
        FIELD_FORMAT = 2
    elif args.threeF:
        FIELD_FORMAT = 3
    elif args.fourF:
        FIELD_FORMAT = 4
    elif args.G_group:
        FIELD_FORMAT = 5
    elif args.P_group:
        FIELD_FORMAT = 6
    elif args.old_format:
        FIELD_FORMAT = 0 # Old format will be assigned 0 for `FIELD_FORMAT`.
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Argument for Field format is not appropriate. Please check them again.\n")
        sys.exit()

    ## Which type of ped file given?
    _p_hped = -1
    _p_hped_descriptor = -1

    if args.hped:
        # Standard 4-field *.ped file given
        _p_hped = args.hped
        _p_hped_descriptor = 1
    elif args.hped_G:
        # G-group *.ped file given
        _p_hped = args.hped_G
        _p_hped_descriptor = 2
    elif args.hped_P:
        # P-group *.ped file given
        _p_hped = args.hped_P
        _p_hped_descriptor = 3
    else:
        # Assuming at least three of them given, there won't be the case which comes to here.
        print(std_ERROR_MAIN_PROCESS_NAME + "Argument for input *.ped file is not appropriate. Please check it again.\n")
        sys.exit()


    NomenCleaner(_p_hped, _p_hped_descriptor, args.iat, args.imgt, args.out, FIELD_FORMAT, __f_NoCaption=args.NoCaption,
                 __leave_NotFound=args.leave_NotFound)

