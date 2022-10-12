#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, join

from Phasing.__main__ import HATK_Phasing
from src.util import Exists

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


def doPhasing(_hStudy, _out_prefix=None):

    bgl_phased = _hStudy.bgl_phased

    if Exists(bgl_phased): # Already phased.
        print(std_WARNING + "There is already phased output('{}'). Skipping phasing.".format(bgl_phased))
        return bgl_phased

    _out_prefix = _out_prefix if _out_prefix else _hStudy.out_prefix+'.PH'

    Phased = HATK_Phasing(_out_prefix, _bfile=_hStudy.bmarker.file_prefix,
                          _java_mem=_hStudy.java_mem, _nthreads=_hStudy.nthreads,
                          _f_save_intermediates=_hStudy.f_save_intermediates)


    return Phased, Phased.bgl_phased