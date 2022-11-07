#-*- coding: utf-8 -*-
import os, sys, re
from os.path import join, basename, dirname
# import numpy as np
# import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from HLA_Heatmap.__main__ import HATK_Heatmap
from src.PLINK_ASSOC import ASSOC
from bMarkerGenerator.bMarker import bMarker
from IMGT2Seq.src.IMGT2Seq_Output import IMGT2Seq_Output
from src.util import Exists

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def doHeatmap(_ASSOC:ASSOC, _bmarker:bMarker, _imgt_output:IMGT2Seq_Output, _out_prefix=None, _f_save_intermediates=False):

    if not _ASSOC:
        print(std_WARNING + "No assoc file was found. Please perform association test first.".format(_ASSOC))
        return None

    if not _imgt_output:
        print(std_WARNING + "No maptable is given. Please specify IMGT2Seq output directory.")
        return None


    _HLAs = _bmarker.l_HLA_target
    _maptable = _imgt_output.HLA_MAPTABLE

    dict_Heatmap = {hla: None for hla in _HLAs} # To return.

    for hla in _HLAs:

        if hla not in _maptable.keys():
            print(std_WARNING + "'HLA-{}' maptable file is not available in IMGT2Seq output directory('{}'). "
                                "Skipping this Heatmap.".format(hla, _imgt_output.out_dir))
            dict_Heatmap[hla] = None
            continue

        if not Exists(_maptable[hla]):
            print(std_WARNING + "'HLA-{}' maptable file can't be found in IMGT2Seq output directory('{}'). "
                                "Skipping this Heatmap".format(hla, _imgt_output.out_dir))
            dict_Heatmap[hla] = None
            continue

        out_prefix = _out_prefix if _out_prefix else _ASSOC.assoc+'.heatmap.HLA-{}'.format(hla)
        dict_Heatmap[hla] = HATK_Heatmap([hla], out_prefix, _maptable[hla], _ASSOC.assoc, _f_save_intermediates)

    dict_Heatmap_pdfs = {hla: _HATK_Heatmap.Heatmap if _HATK_Heatmap else None for hla, _HATK_Heatmap in dict_Heatmap.items()} # To return.

    return dict_Heatmap, dict_Heatmap_pdfs