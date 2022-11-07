# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np
import pandas as pd

from src.util import printDF, printDict


def genHPED(_hat, _N, _out=None, _HLA_req=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"),
            _field_format='2-field', _prop_zero=0.1, _seed=-1, _f_hasGenePrefix=True, _f_hasFieldSep=True):

    """
    Generate a HPED file based on the HAT file. (N rows(samples))

    _field_format ==
        '1-field', '2-field', '3-field', '4-field', 'STANDARD', 'Ggroup', 'Pgroup'
    """

    ### Main variable ###
    which_column = None
    if _field_format == '1-field' or _field_format == '2-field' or _field_format == '3-field' or _field_format == '4-field':
        which_column = 'STANDARD'
    elif _field_format == 'STANDARD' or _field_format == 'Ggroup' or _field_format == 'Pgroup': pass
    else: print("Invalid output field format('{}').".format(_field_format)); return -1

    dict_HAT = HAT2Dict(_hat, which_column)
    HLA_avail = list(dict_HAT.keys())
    HLA_target = None


    ### Main actions ###

    ## `HLA_target`
    f_HLA = np.isin(_HLA_req, HLA_avail)
    HLA_target = np.array(_HLA_req)[f_HLA]
    HLA_discard = np.array(_HLA_req)[~f_HLA]

    print("Target HLAs: ", HLA_target)
    print("Discarded HLAs: ", HLA_discard)


    ## `df_Right` - Main sampling.
    if _seed >=0: np.random.seed(_seed)

    dict_temp = {"{}_{}".format(hla, chrN): None for hla in HLA_target for chrN in np.arange(1, 3)}

    for hla in HLA_target:
        for chrN in np.arange(1, 3):

            k_temp = "{}_{}".format(hla, chrN)

            dict_temp[k_temp] = sampleAllele(hla, dict_HAT[hla], _N, _prop_zero)


    df_Right = pd.DataFrame(dict_temp)
    printDF("df_Right", df_Right)


    ## Tuning the `df_Right`.
    # N-field cutting.
    if _field_format == '1-field' or _field_format == '2-field' or _field_format == '3-field':
        p_Nfield = re.compile(r'\d+[A-Z]?' if _field_format == '1-field' else
                              r'\d+:\d+[A-Z]?' if _field_format == '2-field' else
                              r'\d+:\d+:\d+[A-Z]?')

        func = np.vectorize(lambda x : (p_Nfield.match(x).group() if p_Nfield.match(x) else '0') if x != '0' else x)

        df_Right = pd.DataFrame(func(df_Right.values), index=df_Right.index, columns=df_Right.columns)
        printDF("df_Right('{}')".format(_field_format), df_Right)

    # has the field separator?
    if not _f_hasFieldSep:
        func = np.vectorize(lambda x: re.sub(r':', '', x))
        df_Right = pd.DataFrame(func(df_Right.values), index=df_Right.index, columns=df_Right.columns)
        printDF("df_Right(No field separator", df_Right)

    # Gene Prefix
    if _f_hasGenePrefix:
        df_Right = setGenePrefix(df_Right, HLA_target)
        printDF("df_Right(with Gene prefix)", df_Right)


    ## `df_Left` - Meta Info.
    df_Left = pd.DataFrame({
        'FID': np.arange(_N) + 1,
        'IID': np.arange(_N) + 1,
        'PID': np.repeat('0', _N),
        'MID': np.repeat('0', _N),
        'Sex': np.repeat('-9', _N),
        'Phe': np.repeat('-9', _N)
    })


    df_RETURN = pd.concat([df_Left, df_Right], axis=1)

    if _out:
        df_RETURN.to_csv(_out, sep='\t', header=True, index=False)
        return _out
    else:
        printDF("df_RETURN", df_RETURN)
        return df_RETURN


def HAT2Dict(_hat, _field_format='STANDARD'):

    df_hat = pd.read_csv(_hat, sep='\s+', header=0, dtype=str)
    # printDF("df_hat", df_hat)

    # print(df_hat[['HLA', _field_format]])

    # dict_RETURN = dict(list(df_hat[['HLA', _field_format]].groupby('HLA')[_field_format]))
    dict_RETURN = {hla: sr.values for hla, sr in df_hat[['HLA', _field_format]].groupby('HLA')[_field_format]}
    # printDict(dict_RETURN)
    # print(type(dict_RETURN['A']))

    return dict_RETURN


def sampleAllele(_hla, _arr_HLA, _N, _prop_zero):

    flag_Allele = np.random.choice(2, _N, p=[_prop_zero, 1-_prop_zero]) # either allele or '0'.
    # print(flag_Allele)

    func = np.vectorize(lambda x: np.random.choice(_arr_HLA, 1)[0] if x else '0')
    arr_RETURN = func(flag_Allele)

    # print(arr_RETURN)
    return arr_RETURN


def setGenePrefix(_df_Right, _HLA_target):

    func = np.vectorize(lambda _hla, _allele: '{}*'.format(_hla) + str(_allele) if _allele != '0' else _allele)
    dict_temp = {col: None for col in _df_Right.columns}

    for i in np.arange(len(_HLA_target)):
        hla = _HLA_target[i]
        idx1 = 2 * i
        idx2 = idx1 + 1

        arr_chr1 = func(hla, _df_Right.iloc[:, idx1].values)
        arr_chr2 = func(hla, _df_Right.iloc[:, idx2].values)

        dict_temp['{}_1'.format(hla)] = arr_chr1
        dict_temp['{}_2'.format(hla)] = arr_chr2

    df_RETURN = pd.DataFrame(dict_temp)
    return df_RETURN



if __name__ == '__main__':

    """
    genHPED.py
    
    - Randomly generate `N` samples with the list of given HLA genes.
    - Choose the `field_format` among '1-field', '2-field', '3-field', '4-field', 'Ggroup', or 'Pgroup'
    """

    hat = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220221/HLA_ALLELE_TABLE.imgt3460.hat"
    out = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_NC_20220221/Dummy.N10.hped"
    N = 50
    seed = 123
    field_format = '2-field'
    f_hasGenePrefix = False
    f_hasFieldSep = False
    HLA_req = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DOA", "DOB", "E", "F", "G")

    # HAT2Dict(hat)

    # [hat, N, field_format, out, seed] = sys.argv[1:6]
    # HLA_req = sys.argv[6:]

    t = genHPED(hat, N, out, HLA_req, _seed=seed, _field_format=field_format,
                _f_hasGenePrefix=f_hasGenePrefix, _f_hasFieldSep=f_hasFieldSep)
