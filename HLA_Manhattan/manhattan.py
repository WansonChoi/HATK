# -*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join
import subprocess as sbp

import numpy as np
import pandas as pd
from shutil import which
from math import log10

from src.HATK_Error import HATK_R_Execution_Error
from src.util import Exists, getColumn

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def Manhattan(_assoc_result, _out_prefix, _hg, _pointcol="#778899", _topcol="#FF0000", _min_pos="29.60E6",
              _max_pos="33.2E6", _pointsize="15", _yaxis_unit="10", _p_src='HLA_Manhattan/src/',
              _p_data='HLA_Manhattan/data', _Rscript=which("Rscript"), _f_save_intermediates=False):
    """
    Manhattan plot for PLINK association test result. (ex. *.assoc.{linear,logistic})
    """

    ### Main Variables ###
    # Point color
    _pointcol = re.escape(_pointcol if _pointcol.startswith('#') else ("#"+_pointcol))
    _topcol = re.escape(_topcol if _topcol.startswith('#') else ("#"+_topcol))

    # hg (Human Genome)
    _knownGene = join(_p_data, 'known_genes/known_genes_chr6.hg{0}.txt'.format(_hg))


    ### Main Actions ###

    l_TOP_LABEL = []
    l_yaxis = []

    for i in range(0, len(_assoc_result)):

        # print("\n[{0}] : {1}\n".format(i, _assoc_result[i]))

        t_ar = pd.read_csv(_assoc_result[i], sep='\s+', header=0, usecols=["SNP", "P"]).dropna().sort_values("P")

        # MARKER_set = t_ar.iloc[:, 0].tolist()
        # print(std_MAIN_PROCESS_NAME + "Marker Labels are \n{0}".format(MARKER_set))

        # print(t_ar.head())

        TOP_LABEL = t_ar.iat[0, 0]
        TOP_LABEL_value = t_ar.iat[0, 1]

        # print("\nTop signal is \"{0}(P-value : {1})\"".format(TOP_LABEL, TOP_LABEL_value))


        temp = -log10(TOP_LABEL_value)

        if temp <= 10:
            _yaxis = 10
        else:

            if temp % 5 == 0:
                _yaxis = temp
            else:
                _yaxis = (int(temp/5) + 1)*5


        # print("maximum y-axis is {0}".format(_yaxis))

        l_TOP_LABEL.append(TOP_LABEL)
        l_yaxis.append(str(_yaxis))


    export_ar = ','.join(_assoc_result)
    export_TOP_LABEL = ','.join(l_TOP_LABEL)
    export_yaxis = ','.join(l_yaxis)


    ########## < Plotting Manhattan > ##########

    # command = [_Rscript, join(_p_src, "manhattan_HLA_HATK.R"),
    #            export_ar, _out,
    #            _pointcol, _pointsize, _topcol,
    #            _min_pos, _max_pos, export_TOP_LABEL, export_yaxis, _yaxis_unit,
    #            _knownGene, _p_src, '> {}'.format(_out+'.log')]

    command = \
        "{Rscript} {src} {assoc} {out_prefix} {pointcol} {pointsize} {topcol} {min_pos} {max_pos} {top_label} " \
        "{yaxis} {yaxis_unit} {knownGene} {p_src}" \
        .format(Rscript=_Rscript, src=join(_p_src, "manhattan_HLA_HATK.R"), assoc=export_ar, out_prefix=_out_prefix,
                pointcol=_pointcol, pointsize=_pointsize, topcol=_topcol, min_pos=_min_pos, max_pos=_max_pos,
                top_label=export_TOP_LABEL, yaxis=export_yaxis, yaxis_unit=_yaxis_unit, knownGene=_knownGene, p_src=_p_src)
    # print(command)

    __RETURN__ = plotManhattan(command, _out_prefix, _f_save_intermediates) + ".pdf"

    return __RETURN__



def Manhattan_OM(_assoc_OM_result, _out_prefix, _hg, _pointcol="#778899", _topcol="#FF0000", _min_pos="29.60E6",
                 _max_pos="33.2E6", _pointsize="15", _yaxis_unit="10", _p_src='HLA_Manhattan/src/',
                 _p_data='HLA_Manhattan/data', _Rscript=which("Rscript"), _f_save_intermediates=False):
    """
    Manhattan plot for Omnibus test result. (ex. *.omnibus)
    """

    ### Main Variables ###
    # Point color
    _pointcol = re.escape(_pointcol if _pointcol.startswith('#') else ("#"+_pointcol))
    _topcol = re.escape(_topcol if _topcol.startswith('#') else ("#"+_topcol))

    # Patterns to use
    p_Omnibus_AA = re.compile(r'AA_([A-Z0-9]+)_(-?\d+)')
    p_Omnibus_INS = re.compile(r'INS_AA_([A-Z0-9]+)_(-?\d+)x(-?\d+)')

    # HLA target
    _HLA_target = getTargetHLA(_assoc_OM_result, p_Omnibus_AA, p_Omnibus_INS)

    __RETURN__ = []


    ### Main Actions ###
    for i in range(0, len(_assoc_OM_result)):

        # Header : [Variant, deltaDeviance, deltaDF, N, log10_P, Residues]
        t_ar = pd.read_csv(_assoc_OM_result[i], sep='\s+', header=0) \
                    .dropna()

        # print(t_ar['log10_P'].map(lambda x : 10**x))
        t_ar['P'] = t_ar['log10_P'].map(lambda x : 10**x)


        # Extracting 'HLA' and 'Relative position' of Amino acid marker
        f_AA = t_ar['Variant'].str.match(p_Omnibus_AA)
        df_AA = t_ar['Variant'][f_AA].str.extract(p_Omnibus_AA, expand=True)
        df_AA.columns = ['HLA', 'REL_POS']
        # print("df_AA :\n{}\n".format(df_AA))


        # Extracting 'HLA' and 'Relative position' of INS a.a. marker
        f_INS = t_ar['Variant'].str.match(p_Omnibus_INS)

        if f_INS.any():
            df_INS = t_ar['Variant'][f_INS].str.extract(p_Omnibus_INS, expand=True)
            df_INS = pd.concat([
                df_INS.iloc[:, 0],
                df_INS.iloc[:, [1,2]].apply(lambda x : str(x.astype(int).mean()), axis=1)
            ], axis=1)
        else:
            df_INS = pd.DataFrame([]) # Null DataFrame

        df_INS.columns = ['HLA', 'REL_POS']

        # print("df_INS :\n{}\n".format(df_INS))

        df_HLA_RelPos = pd.concat([df_AA, df_INS], axis=0).sort_index()
        df_HLA_RelPos['REL_POS'] = df_HLA_RelPos['REL_POS'].map(lambda x : float(x)) # as float

        t_ar = pd.concat([t_ar, df_HLA_RelPos], axis=1)
        # print(t_ar)


        for i in range(len(_HLA_target)):

            ### Plotting by each HLA gene.

            t_ar_byHLA = t_ar[t_ar['HLA'] == _HLA_target[i]]
            # print(t_ar_byHLA.head())

            t_ar_byHLA_2 = t_ar_byHLA[['Variant', 'P', 'REL_POS']].sort_values('P')
            # print(t_ar_byHLA_2.head())


            TOP_LABEL = t_ar_byHLA_2['Variant'].iat[0]
            TOP_LABEL_value = t_ar_byHLA_2['P'].iat[0]

            # print("\nTop signal is \"{0}(P-value : {1})\"".format(TOP_LABEL, TOP_LABEL_value))

            temp = -log10(TOP_LABEL_value)

            if temp <= 10:
                _yaxis = 10
            else:

                if temp % 5 == 0:
                    _yaxis = temp
                else:
                    _yaxis = (int(temp / 5) + 1) * 5

            # print("maximum y-axis is {0}".format(_yaxis))


            _forManhattn = _out_prefix + '.{}.forManhattan.txt'.format(_HLA_target[i])
            # print(t_ar_byHLA.sort_values('REL_POS'))
            t_ar_byHLA.sort_values('REL_POS').to_csv(_forManhattn, sep='\t', header=True, index=False)

            export_ar = _forManhattn
            export_TOP_LABEL = TOP_LABEL
            export_yaxis = str(_yaxis)

            _min_pos = str(t_ar_byHLA['REL_POS'].min())
            _max_pos = str(t_ar_byHLA['REL_POS'].max())



            ########## < Plotting Manhattan > ##########

            # print(std_MAIN_PROCESS_NAME + "Plotting Manhattan.\n")

            OUT_prefix = _out_prefix + '.{}'.format(_HLA_target[i])

            # command = [_Rscript, os.path.join(_p_src, "manhattan_HLA_HATK.Omnibus.R"),
            #            export_ar, OUT_prefix,
            #            _pointcol, _pointsize, _topcol,
            #            _min_pos, _max_pos, export_TOP_LABEL, export_yaxis, _yaxis_unit,
            #            _p_src, '> {}'.format(OUT_prefix + '.log')]

            command = \
                "{Rscript} {src} {assoc} {out_prefix} {pointcol} {pointsize} {topcol} {min_pos} {max_pos} {top_label} " \
                "{yaxis} {yaxis_unit} {p_src}" \
                .format(Rscript=_Rscript, src=join(_p_src, "manhattan_HLA_HATK.Omnibus.R"), assoc=export_ar, out_prefix=OUT_prefix,
                        pointcol=_pointcol, pointsize=_pointsize, topcol=_topcol, min_pos=_min_pos, max_pos=_max_pos,
                        top_label=export_TOP_LABEL, yaxis=export_yaxis, yaxis_unit=_yaxis_unit, p_src=_p_src)
            # print(command)

            r = plotManhattan_OM(command, OUT_prefix, _f_save_intermediates, [_forManhattn])
            __RETURN__.append(r+'.pdf')


    return __RETURN__



def plotManhattan(_command, _out_prefix, _f_save_intermediates=False, _l_to_remove:list=()):
    try:
        with open(_out_prefix + '.log', 'w') as f_out_log:
            sbp.run(_command.split(), check=True, stdout=f_out_log, stderr=f_out_log)

    except sbp.CalledProcessError:
        raise HATK_R_Execution_Error(
            "\nFailed to plot the Manhattan plot. ('{}')".format(_command) + \
            ("\nPlease refer to its log file. ('{}')".format(_out_prefix + '.log') if Exists(_out_prefix + '.log') else "")
        )
    else:
        # remove
        if not _f_save_intermediates:
            if Exists(_out_prefix + '.log'): os.remove(_out_prefix + '.log')
            for item in _l_to_remove:
                if Exists(item): os.remove(item)

        return _out_prefix



def plotManhattan_OM(_command, _out_prefix, _f_save_intermediates=False, _l_to_remove:list=()):
    try:
        with open(_out_prefix + '.log', 'w') as f_out_log:
            sbp.run(_command.split(), check=True, stdout=f_out_log, stderr=f_out_log)

    except sbp.CalledProcessError:
        # raise HATK_R_Execution_Error(
        #     "\nFailed to plot the Omnibus Manhattan plot. ('{}')".format(_command) + \
        #     ("\nPlease refer to its log file. ('{}')".format(_out_prefix + '.log') if Exists(_out_prefix + '.log') else "")
        # )
        print(std_ERROR + "\nFailed to plot the Omnibus Manhattan plot. ('{}')".format(_command) + \
            ("\nPlease refer to its log file. ('{}')".format(_out_prefix + '.log') if Exists(_out_prefix + '.log') else ""))
    else:
        # remove
        if not _f_save_intermediates:
            if Exists(_out_prefix + '.log'): os.remove(_out_prefix + '.log')
            for item in _l_to_remove:
                if Exists(item): os.remove(item)

        return _out_prefix



def getUniqueHLA(_assoc_OM, _p_Omnibus_AA, _p_Omnibus_INS):

    l_RETURN = []

    for item in getColumn(_assoc_OM, 0):
        m_Omnibus_AA = _p_Omnibus_AA.match(item)
        m_Omnibus_INS = _p_Omnibus_INS.match(item)

        if bool(m_Omnibus_AA): l_RETURN.append(m_Omnibus_AA.group(1))
        elif bool(m_Omnibus_INS): l_RETURN.append(m_Omnibus_INS.group(1))
        else: pass

    return l_RETURN

def getTargetHLA(_assoc_OM_result, _p_Omnibus_AA, _p_Omnibus_INS):
    arr_whole_HLA = np.array(
        [getUniqueHLA(_assoc_OM, _p_Omnibus_AA, _p_Omnibus_INS) for _assoc_OM in _assoc_OM_result]
    ).flatten()

    return list(np.unique(arr_whole_HLA))
