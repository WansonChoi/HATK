#-*- coding: utf-8 -*-
import os, sys, re
from os.path import join, basename, dirname
# import numpy as np
# import pandas as pd

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from HLA_Heatmap.__main__ import HATK_Heatmap
from src.PLINK_ASSOC import ASSOC
from src.util import Exists

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def doHeatmap(_ASSOC:ASSOC, _HLA, _maptable, _out_prefix=None, _f_save_intermediates=False):

    if not _ASSOC:
        print(std_WARNING + "No assoc file was found. Please perform association test first.".format(_ASSOC))
        return None

    if not Exists(_maptable):
        print(std_WARNING + "Maptable file for HLA-{}('{}') can't be found.".format(_HLA, _maptable))
        return None

    out_prefix = _out_prefix if _out_prefix else _ASSOC.assoc+'.heatmap.HLA-{}'.format(_HLA)

    return HATK_Heatmap([_HLA], out_prefix, _maptable, _ASSOC.assoc, _f_save_intermediates)