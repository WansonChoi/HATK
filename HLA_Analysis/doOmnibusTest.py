#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join

from OmnibusTest.__main__ import HATK_OmibusTest
from bMarkerGenerator.bMarker import bMarker
from src.PLINK_PHENO import PHENO
from src.PLINK_COVAR import COVAR
from src.PLINK_CONDITION import CONDITION

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def doOmnibusTest(_out_prefix, _bmarker:bMarker, _PHENO:PHENO, _phe_name, _COVAR:COVAR, _CONDITION:CONDITION,
                  _bgl_phased=None, _aa=None, _f_save_intermediates=False):

    _pheno = _PHENO.phe

    if _COVAR:
        _covar = _COVAR.file
        _covar_name =_COVAR.covar_name_target
    else:
        _covar = None
        _covar_name = None
        # Redundant, but can't bring "_COVAR.file" if _COVAR is None.

    if _CONDITION:
        _condition = _CONDITION.condition
        _condition_list = _CONDITION.condition_list
    else:
        _condition = None
        _condition_list = None


    return HATK_OmibusTest(_out_prefix, _bmarker.file_prefix, _pheno, _phe_name, _aa, _bgl_phased,
                           _covar, _covar_name, _condition, _condition_list, _f_save_intermediates)


