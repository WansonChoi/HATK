# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join, isdir

from src.HATK_Error import HATK_InputPreparation_Error, HATK_PLINK_Execution_Error, RaiseError
from OmnibusTest.__main__ import HATK_OmibusTest
from src.util import Exists
from HLA_Manhattan.__main__ import HATK_Manhattan

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def doManhattanPlot_OM(_OMNIBUS:HATK_OmibusTest, _hg, _out_prefix=None, _f_save_intermediates=False):

    if not _OMNIBUS:
        print(std_WARNING + "No omnibus test result was found. Please perform the Omnibus test first. Skipping Manhattan plot.")
        return 0

    out_prefix = _out_prefix if _out_prefix else _OMNIBUS.OMNIBUS+".manhattan"

    return HATK_Manhattan([_OMNIBUS.OMNIBUS], out_prefix, _hg, _f_save_intermediates=_f_save_intermediates)
