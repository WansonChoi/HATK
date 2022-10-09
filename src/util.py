#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join
from shutil import which

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError

"""
Collection of small utility functions
"""

Exists = lambda x, func=exists: bool(x) and func(x)
checkFile = lambda x, _msg="" : x if Exists(x) else RaiseError(HATK_InputPreparation_Error, _msg)

def checkFile2(_file, _msg):
    """
    To Accept 'None' value.

    In the 'COVAR' and 'CONDITION" class, '_file = None' should be acceptable.
    """
    if _file == None or Exists(_file):
        return _file
    else:
        RaiseError(HATK_InputPreparation_Error, _msg)


def printDF(_df_name, _df):
    print("{}:\n{}".format(_df_name, _df))


def printDict(_dict, _N=-1):
    count = 0
    for k, v in _dict.items():
        print("{}: {}".format(k, v))

        if _N > 0:
            count += 1
            if count >= _N: break


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


def FieldFormat2Label(_which_format):
    if _which_format == 1: return "1-field"
    elif _which_format == 2: return "2-field"
    elif _which_format == 3: return "3-field"
    elif _which_format == 4: return "4-field"
    elif _which_format == 5: return "G-group"
    elif _which_format == 6: return "P-group"
    else: return None


def findExec(_cmd, _msg):
    if Exists(which(_cmd)):
        return which(_cmd)
    elif Exists(join("./dependency", _cmd)):
        return join("./dependency", _cmd)
    elif Exists(join("../dependency", _cmd)):
        return join("../dependency", _cmd)
    else:
        RaiseError(HATK_InputPreparation_Error, _msg)