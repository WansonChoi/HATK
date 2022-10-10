# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join, isdir

from src.HATK_Error import HATK_InputPreparation_Error, HATK_PLINK_Execution_Error, RaiseError
from src.PLINK_ASSOC import ASSOC
from src.util import Exists
from HLA_Manhattan.__main__ import HATK_Manhattan

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def doManhattanPlot_assoc(_ASSOC:ASSOC, _hg, _f_save_intermediates=False):

    if not _ASSOC:
        print(std_WARNING + "No assoc file was found. Please perform association test first.".format(_ASSOC))
        return None

    return HATK_Manhattan([_ASSOC.assoc], _ASSOC.assoc + '.manhattan', _hg, _f_save_intermediates=_f_save_intermediates)

