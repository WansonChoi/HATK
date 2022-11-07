#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join, exists
import pandas as pd
import numpy as np

from src.util import getColumn


def isBIMwithOldLabel(_bim):

    p_HLA = re.compile(r'^HLA_')
    p_HLA_old = re.compile(r'^HLA_[A-Z]+\d?_\d+[A-Z]?$')
    # p_HLA_HATK = re.compile(r'^HLA_[A-Z]+\d?\*\d{2,3}(:\d{2,3})*[A-Z]?$')

    l_is_old_label = [bool(p_HLA_old.match(_label)) for _label in getColumn(_bim, 1) if p_HLA.match(_label)]

    return (float(sum(l_is_old_label)) / len(l_is_old_label)) > 0.8



def updateNameANDAlleles(_bim, _out_update_name, _out_update_alleles):

    p_HLA = re.compile(r'^HLA_')

    with open(_bim, 'r') as f_bim, \
            open(_out_update_name, 'w') as f_out_update_name, \
            open(_out_update_alleles, 'w') as f_out_update_alleles:

        for line in f_bim:
            """
            6	HLA_F*01:01	0	29799220	p	a
            6	HLA_F*01:02	0	29799220	p	a
            (=>)
            6	HLA_F_0101	0	29799220	P	A
            6	HLA_F_0102	0	29799220	P	A
            
            
            """

            l = line.split()

            if p_HLA.match(l[1]):
                ## update_name
                f_out_update_name.write('\t'.join([l[1], toOldLabel(l[1])]) + "\n")

                ## update_allele
                f_out_update_alleles.write(
                    # '\t'.join([l[1], l[4], l[5], l[4].capitalize(), l[5].capitalize()]) + "\n"
                    '\t'.join([toOldLabel(l[1]), l[4], l[5], l[4].capitalize(), l[5].capitalize()]) + "\n"
                ) # New names are needed, not the former ones.

    return _out_update_name, _out_update_alleles



def toOldLabel(_label:str) -> str:
    # replace asterisk to underscore
    str_temp1 = re.sub(r'\*', '_', _label)

    # remove colon(':')
    str_temp2 = re.sub(r':', '', str_temp1)

    return str_temp2





if __name__ == '__main__':

    bim_old = "/home/wansonchoi/sf_VirtualBox_Share/UC-CD-HLA/data/Merged/merged.bim"
    bim_HATK = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATKv2_wholeImple2_20221012/wtccc_58C_RA.hatk.300+300.hg18.chr6.29-34mb.ALL.bim"

    print(isBIMwithOldLabel(bim_old))
    print(isBIMwithOldLabel(bim_HATK))
