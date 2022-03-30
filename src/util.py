#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError

"""
Collection of small utility functions
"""

Exists = lambda x, func=exists: x and func(x)
checkFile = lambda x, _msg="" : x if Exists(x) else RaiseError(HATK_InputPreparation_Error, _msg)


def printDF(_df_name, _df):
    print("{}:\n{}".format(_df_name, _df))


def printDict(_dict):
    for k, v in _dict.items():
        print("{}: {}".format(k, v))


def hasHeader(_file, _pattern=re.compile(r'^FID\s+IID\s+')) -> bool:
    """
    To check whether given file(ex. pheno, covar, etc.) has the header or not.
    """
    with open(_file, 'r') as f_file:
        line_1st = f_file.readline()
        return bool(_pattern.match(line_1st))


def getHeader(_file) -> list:
    """
    Get the 1st line, i.e. header line.
    """
    with open(_file, 'r') as f_file:
        return f_file.readline().split()


def get_N_row(_file) -> int:
    with open(_file, 'r') as f_file:
        count = 0
        for _ in f_file: count +=1
        return count


def getColumn(_file, _N_col):
    with open(_file, 'r') as f_input:
        for line in f_input:
            l = line.split()
            yield l[_N_col]