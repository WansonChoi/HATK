#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join, exists

from src.util import getColumn

def exclude_AA_SNPS_INS(_bim, _out):

    p_AA = re.compile(r'^AA_([A-Z]+\d?)')
    p_SNPS = re.compile(r'^SNPS?_([A-Z]+\d?)_')
    p_INS = re.compile(r'^INS_(AA_|SNPS_)?([A-Z]+\d?)')
    # Above patterns are the same ones in the 'checkBIM()' function in the 'bMarker.py'.

    with open(_out, 'w') as f_out:
        for label in getColumn(_bim, 1):
            if p_AA.match(label) or p_SNPS.match(label) or p_INS.match(label):
                f_out.write(label + "\n")

    return _out