# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
header_ped = ["FID", "IID", "PID", "MID", "Sex", "Phe"]

# Patterns
p_Suffix = re.compile(r'.+[A-Z]$')


def NomenCleaner(_hped, _hat, _imgt, _out, __f_NoCaption=False, __leave_NotFound=False, **kwargs):


    ##### < Output Format > #####

    OUTPUT_FORMAT = -1

    if kwargs['__oneF']:
        OUTPUT_FORMAT = 1
    elif kwargs['__twoF']:
        OUTPUT_FORMAT = 2
    elif kwargs['__threeF']:
        OUTPUT_FORMAT = 3
    elif kwargs['__fourF']:
        OUTPUT_FORMAT = 4
    elif kwargs['__Ggroup']:
        OUTPUT_FORMAT = 5
    elif kwargs['__Pgroup']:
        OUTPUT_FORMAT = 6




    ##### < Loading Data > #####

    # *.hped
    __HPED__ = pd.read_csv(_hped, sep='\s+', header=None, dtype=str)
    # print("__HPED__ :\n{}\n".format(__HPED__.head()))

    # *.hat
    __HAT__ = pd.read_csv(_hat, sep='\s+', header=0, dtype=str, index_col=0)
    # print("__HAT__ :\n{}\n".format(__HAT__.head()))

    d__HAT__ = {HLA_names[i]: __HAT__.loc[HLA_names[i], :] for i in range(len(HLA_names))}

    # for k,v in d__HAT__.items():
    #     print("HLA : {}\n{}\n".format(k, v.head()))



    ##### < Main iteration > #####

    f_chped = open(_out+'.chped', 'w')
    f_chped_log = open(_out+'.chped.log', 'w')

    count = 0

    for eachRow in __HPED__.itertuples():

        print(eachRow)

        [t_FID, t_IID] = eachRow[1:3]

        t_alleles = eachRow[7:]
        # print(t_alleles)

        ### [1] Allele conversion

        l_alleles_row = []

        for i in range(len(HLA_names)):

            idx1 = 2*i
            idx2 = 2*i + 1

            if t_alleles[idx1] != '0':

                [t_converted_allele1, LOG_MESSAGE1] = getConvertedAllele(HLA_names[i], t_alleles[idx1], d__HAT__[HLA_names[i]], OUTPUT_FORMAT, __leave_NotFound)

                print("Converted Alleles : {}".format(t_converted_allele1))
                print("LOG_MESSAGE :{}\n".format(LOG_MESSAGE1))

            else:

                t_converted_allele1 = t_alleles[idx1]
                LOG_MESSAGE1 = 'No conversion ({})'.format(t_alleles[idx1])


            if t_alleles[idx2] != '0':

                [t_converted_allele2, LOG_MESSAGE2] = getConvertedAllele(HLA_names[i], t_alleles[idx2], d__HAT__[HLA_names[i]], OUTPUT_FORMAT, __leave_NotFound)

                print("Converted Alleles : {}".format(t_converted_allele2))
                print("LOG_MESSAGE :{}\n".format(LOG_MESSAGE2))

            else:
                t_converted_allele2 = t_alleles[idx2]
                LOG_MESSAGE2 = 'No conversion ({})'.format(t_alleles[idx2])




            l_alleles_row.append(('%s*'%(HLA_names[i]) if (not __f_NoCaption and t_converted_allele1 != '0') else '') + t_converted_allele1)
            l_alleles_row.append(('%s*'%(HLA_names[i]) if (not __f_NoCaption and t_converted_allele2 != '0') else '') + t_converted_allele2)

            ### [2] logging the result
            f_chped_log.write('\t'.join([t_FID, t_IID, HLA_names[i]+'_1', LOG_MESSAGE1])+'\n')
            f_chped_log.write('\t'.join([t_FID, t_IID, HLA_names[i]+'_2', LOG_MESSAGE2])+'\n')



            """
            G-group 2 UPDATED
            P-group 2 UPDATED
            이거 기능 두 개 아직 완성 못함.
            
            """


        ### Generating chped
        f_chped.write('\t'.join(list(eachRow[1:7]) + l_alleles_row) + '\n')


        count += 1
        # if count > 1 : break


    f_chped.close()
    f_chped_log.close()

    return _out+'.chped'




def getConvertedAllele(_hla, _allele, _d__HAT__, _OUTPUT_FORMAT, __leave_NotFound):

    """

    """

    __RETURN__ = None # Converted allele
    LOG_MESSAGE = None



    ## OUTPUT format classification.

    isGgroup = _allele.endswith('G')
    isPgroup = _allele.endswith('P')


    if isGgroup:

        if _OUTPUT_FORMAT == 5:
            # No need conversion.
            __RETURN__ = _allele
            LOG_MESSAGE = "[Need not conversion] {} -> {}".format(_allele, _allele)
            return [__RETURN__, LOG_MESSAGE]


        [__RETURN__, LOG_MESSAGE] = get1stAllele(_allele, _d__HAT__, 'Ggroup', _OUTPUT_FORMAT, __leave_NotFound)



    elif isPgroup:

        if _OUTPUT_FORMAT == 6:
            # No need conversion.
            __RETURN__ = _allele
            LOG_MESSAGE = "[Need not conversion] {} -> {}".format(_allele, _allele)
            return [__RETURN__, LOG_MESSAGE]


        [__RETURN__, LOG_MESSAGE] = get1stAllele(_allele, _d__HAT__, 'Pgroup', _OUTPUT_FORMAT, __leave_NotFound)



    else:

        if ':' in _allele:

            ### Must be in UPDATED nomenclature.

            # (1) Exact as it is (ex.
            # (2) as a Prefix of some HLA alleles
            # ex) 03:28 (Possible candidates : ['03:28', '03:280'])

            [__RETURN__, LOG_MESSAGE] = get1stAllele(_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)


        else:

            ### Possibly OLD nomenclature[__RETURN__, LOG_MESSAGE]

            if len(_allele) <= 3:



                # Exception handling for 1field allele (ex. DPB1*03 -> 0302(-> 102:01; OLD) / 03:01:01:01(UPDATED))

                flag_matched_UPDATED = _d__HAT__.loc[:, 'STANDARD'].str.match(_allele)
                # Check both first.

                if flag_matched_UPDATED.any():

                    [__RETURN__, LOG_MESSAGE] = get1stAllele(_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)

                else:

                    [__RETURN__, LOG_MESSAGE] = get1stAllele(_allele, _d__HAT__, 'OLD', _OUTPUT_FORMAT, __leave_NotFound)




            else:

                ### Possibly OLD nomenclature

                flag_matched_OLD = _d__HAT__.loc[:, 'OLD'].str.match(_allele)

                if flag_matched_OLD.any():

                    [__RETURN__, LOG_MESSAGE] = get1stAllele(_allele, _d__HAT__, 'OLD', _OUTPUT_FORMAT, __leave_NotFound)


                else:

                    ### Must be in UPDATED nomenclature.


                    ## Checking Suffix
                    hasSuffix = p_Suffix.match(_allele)

                    if hasSuffix:
                        _allele2 = _allele[:-1]
                        _allele_Suffix = _allele[-1]
                    else:
                        _allele2 = _allele
                        _allele_Suffix = -1


                    ### Main Digit-checking

                    sr_UPDATED = _d__HAT__.loc[:, 'STANDARD'].str

                    if len(_allele2) == 4:
                        # 2 + 2 (+C) = 4

                        trial1_allele = ':'.join([_allele2[:2], _allele2[2:4]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)




                    elif len(_allele2) == 5:

                        # (1) 2 + 3 (+C) = 5
                        # (2) 3 + 2 (+C) = 5

                        trial1_allele = ':'.join([_allele2[:2], _allele2[2:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            trial2_allele = ':'.join([_allele2[:3], _allele2[3:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial2 = sr_UPDATED.match(trial2_allele)

                            if Flag_trial2.any():
                                [__RETURN__, LOG_MESSAGE] = get1stAllele(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                            else:
                                __RETURN__ = _allele if __leave_NotFound else "0"
                                if __leave_NotFound:
                                    LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                else:
                                    LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

                    elif len(_allele2) == 6:

                        # (1) 2 + 2 + 2 (+C)
                        # (2) 3 + 3 (+C)


                        # (1) 2 + 2 + 2 (+C)
                        trial1_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            # (2) 3 + 3 (+C)
                            trial2_allele = ':'.join([_allele2[:3], _allele2[3:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial2 = sr_UPDATED.match(trial2_allele)

                            if Flag_trial2:
                                [__RETURN__, LOG_MESSAGE] = get1stAllele(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                            else:
                                __RETURN__ = _allele if __leave_NotFound else "0"
                                if __leave_NotFound:
                                    LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                else:
                                    LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                    elif len(_allele2) == 7:

                        # (1) 2 + 3 + 2 (+C) = 7
                        # (2) 2 + 2 + 3 (+C) = 7
                        # (3) 3 + 2 + 2 (+C) = 7


                        # (1) 2 + 3 + 2 (+C) = 7
                        trial1_allele = ':'.join([_allele2[:2], _allele2[2:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            # (2) 2 + 2 + 3 (+C) = 7
                            trial2_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial2 = sr_UPDATED.match(trial2_allele)

                            if Flag_trial2:
                                [__RETURN__, LOG_MESSAGE] = get1stAllele(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                            else:
                                # (3) 3 + 2 + 2 (+C) = 7
                                trial3_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                                Flag_trial3 = sr_UPDATED.match(trial3_allele)

                                if Flag_trial3.any():
                                    [__RETURN__, LOG_MESSAGE] = get1stAllele(trial3_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                                else:
                                    __RETURN__ = _allele if __leave_NotFound else "0"
                                    if __leave_NotFound:
                                        LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                    else:
                                        LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                    elif len(_allele2) == 8:

                        # (1) 2 + 2 + 2 + 2 (+C)
                        # (2) 3 + 2 + 3 (+C)
                        # (3) 3 + 3 + 2 (+C)


                        # (1) 2 + 2 + 2 + 2
                        trial1_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:6], _allele2[6:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            # (2) 3 + 2 + 3
                            trial2_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial2 = sr_UPDATED.match(trial2_allele)

                            if Flag_trial2:
                                [__RETURN__, LOG_MESSAGE] = get1stAllele(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                            else:
                                # (3) 3 + 3 + 2
                                trial3_allele = ':'.join([_allele2[:3], _allele2[3:6], _allele2[6:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                                Flag_trial3 = sr_UPDATED.match(trial3_allele)

                                if Flag_trial3.any():
                                    [__RETURN__, LOG_MESSAGE] = get1stAllele(trial3_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                                else:
                                    __RETURN__ = _allele if __leave_NotFound else "0"
                                    if __leave_NotFound:
                                        LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                    else:
                                        LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                    elif len(_allele2) == 9:

                        # (1) 2 + 3 + 2 + 2 (+C)
                        # (2) 3 + 2 + 2 + 2 (+C)

                        # (1) 2 + 3 + 2 + 2
                        trial1_allele = ':'.join([_allele2[:2], _allele2[2:5], _allele2[5:7], _allele2[7:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            # (2) 3 + 2 + 2 + 2
                            trial2_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:7], _allele2[7:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial2 = sr_UPDATED.match(trial2_allele)

                            if Flag_trial2.any():
                                [__RETURN__, LOG_MESSAGE] = get1stAllele(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                            else:
                                __RETURN__ = _allele if __leave_NotFound else "0"
                                if __leave_NotFound:
                                    LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                else:
                                    LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                    elif len(_allele2) == 10:

                        # (1) 3 + 3 + 2 + 2

                        trial1_allele = ':'.join([_allele2[0:3], _allele2[3:6], _allele2[6:8], _allele2[8:10]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial1 = sr_UPDATED.match(trial1_allele)

                        if Flag_trial1.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound)
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)




        if _OUTPUT_FORMAT <= 3:

            if _OUTPUT_FORMAT == 1:
                p_OUTPUT_FORMAT = re.compile(r'\d{2,3}')
            elif _OUTPUT_FORMAT == 2:
                p_OUTPUT_FORMAT = re.compile(r'\d{2,3}:\d{2,3}[A-Z]?')
            elif _OUTPUT_FORMAT == 3:
                p_OUTPUT_FORMAT = re.compile(r'\d{2,3}:\d{2,3}:\d{2,3}[A-Z]?')


            m_field = p_OUTPUT_FORMAT.match(__RETURN__)

            __RETURN__ = m_field.group() if bool(m_field) else __RETURN__





    return [__RETURN__, LOG_MESSAGE]




def get1stAllele(_hped_allele, _df, _from, _OUTPUT_FORMAT, __leave_NotFound=False):

    """
    Only for UPDATED nomenclature;
    """

    ### OUTPUT field format
    _to = 1 if _OUTPUT_FORMAT <= 4 else 3 if _OUTPUT_FORMAT == 5 else 4 if _OUTPUT_FORMAT == 6 else -1


    __1st_Allele__ = None
    LOG_MESSAGE = None


    sr_from = _df.loc[:, _from]


    ### Matching job.

    ## Exact match
    Flag_exact_match = sr_from.str.match(_hped_allele+'$')
    f_exact = False
    f_single_candidate = False

    if Flag_exact_match.any():

        # (ex. '03:28' -> '03:28')

        df_candidate = _df.loc[Flag_exact_match, :]
        __1st_Allele__ = df_candidate.iat[0, _to]

        f_exact = True
        f_single_candidate = (df_candidate.shape[0] == 1)

    else:

        ## Just match
        Flag_match = sr_from.str.match(_hped_allele)

        if Flag_match.any():

            # (ex. '03:28' -> '03:280')

            df_candidate = _df.loc[Flag_match, :]
            __1st_Allele__ = df_candidate.iat[0, _to]

            f_single_candidate = (df_candidate.shape[0] == 1)

        else:

            __1st_Allele__ = _hped_allele if __leave_NotFound else '0'

            if __leave_NotFound:
                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_hped_allele, _hped_allele)
            else:
                LOG_MESSAGE = "[Not Found] {} -> 0".format(_hped_allele)

            return [__1st_Allele__, LOG_MESSAGE]



    ### Log message
    if _from == 'OLD':

        if _to == 1:
            # to Standard
            LOG_MESSAGE = "[{}] {} -> {}(OLD) -> {}".format(
                'Exact match' if f_exact else 'Approximated', _hped_allele, df_candidate.iat[0, 2], __1st_Allele__)



            if not f_exact:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [2, 1]].apply(' -> '.join, axis=1).tolist())])

        elif _to == 3 or _to == 4:
            # to Ggroup
            LOG_MESSAGE = "[{}] {} -> {}(OLD) -> {} -> {}".format(
                'Exact match' if f_exact else 'Approximated', _hped_allele, df_candidate.iat[0, 2], df_candidate.iat[0, 1], __1st_Allele__)

            if not f_exact:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [2, 1, _to]].apply(' -> '.join, axis=1).tolist())])


    else:

        if _to == 1:
            # to Standard
            LOG_MESSAGE = "[{}] {} -> {}".format(
                'Exact match' if f_exact else 'Approximated', _hped_allele, __1st_Allele__)

            if not f_exact:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, 1].tolist())])

        elif _to == 3 or _to == 4:
            # to Ggroup
            LOG_MESSAGE = "[{}] {} -> {} -> {}".format(
                'Exact match' if f_exact else 'Approximated', _hped_allele, df_candidate.iat[0, 1], __1st_Allele__)

            if not f_exact:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [1, _to]].apply(' -> '.join, axis=1).tolist())])



    return [__1st_Allele__, LOG_MESSAGE]




if __name__ == '__main__':


    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        < NomenCleaner.v3.py >
        
        - Transforms *.hped file to *.chped file.
        - *.hat file must be given as input file.






    #########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    # Input (1) : *.ped file
    PED_TYPE = parser.add_mutually_exclusive_group(required=True)
    PED_TYPE.add_argument("--hped", help="\nHLA Type Data with raw HLA allele(ex. 0101).\n\n", dest="hped")
    # PED_TYPE.add_argument("--hped-Ggroup", help="\nHLA Type Data with raw G-group allele(ex. 020601G).\n\n",
    #                       dest="hped_G")
    # PED_TYPE.add_argument("--hped-Pgroup", help="\nHLA Type Data with raw P-group allele (ex. 0102P).\n\n",
    #                       dest="hped_P")

    # Input (2) : *.iat file
    parser.add_argument("-hat", help="\nHLA Allele Table file(*.hat).\n\n", required=True)
    parser.add_argument("-imgt", help="\nSpecifying the IMGT-HLA version.\n\n", required=True)

    # Ouptut Prefix
    parser.add_argument("--out", "-o", help="\nOutput file prefix.\n\n", required=True)
    parser.add_argument("--leave-NotFound",
                        help="\nLeaving HLA alleles which can't be found in given *.iat file(Novel or Erroneous allele) intact.\n\n",
                        action='store_true')

    # Output format
    output_digit_selection = parser.add_mutually_exclusive_group(required=True)
    output_digit_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n",
                                        action="store_true", dest="oneF")
    output_digit_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n",
                                        action="store_true", dest="twoF")
    output_digit_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n",
                                        action="store_true", dest="threeF")
    output_digit_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields\n\n",
                                        action="store_true", dest="fourF")
    output_digit_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n",
                                        action="store_true")
    output_digit_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n",
                                        action="store_true")

    # Flag to remove HLA gene caption.
    parser.add_argument("--NoCaption", help="\nConverted HLA alleles NOT to have HLA gene prefix(ex. \"A*\").\n\n", action='store_true')



    ##### <for Test> #####

    # in OS X
    _hped = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.Ggroup.chped'
    _hat = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/HLA_ALLELE_TABLE.imgt3320.hat'
    _out = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.Ggroup.20190924'

    # in Ubuntu
    # _hped = '/home/wanson/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.ruined.hped'
    # _hat = '/home/wanson/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/HLA_ALLELE_TABLE.imgt3320.hat'
    # _out = '/home/wanson/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.ruined.20190922'

    ## OLD or UPDATED alleles.
    # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out+'.1field', '--1field', '-imgt', '3320'])
    # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out+'.2field', '--2field', '-imgt', '3320'])
    # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out+'.3field', '--3field', '-imgt', '3320'])
    # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out+'.4field', '--4field', '-imgt', '3320'])

    ## OLD or UPDATED alleles (no colon).
    # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '--4field', '-imgt', '3320'])
    # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out+'.Ggroup', '--Ggroup', '-imgt', '3320'])
    args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out+'.Pgroup', '--Pgroup', '-imgt', '3320'])


    ##### <for Publication> #####

    # args = parser.parse_args()
    print(args)


    NomenCleaner(args.hped, args.hat, args.imgt, args.out,
                 __f_NoCaption=args.NoCaption, __leave_NotFound=args.leave_NotFound,
                 __oneF=args.oneF, __twoF=args.twoF, __threeF=args.threeF, __fourF=args.fourF, __Ggroup=args.Ggroup, __Pgroup=args.Pgroup)