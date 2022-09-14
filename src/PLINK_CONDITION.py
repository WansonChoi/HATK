#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists
import pandas as pd
import numpy as np

# from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
# from src.util import Exists, checkFile, hasHeader, getHeader

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)



class CONDITION(object):

    def __init__(self, _condition:str=None, _condition_list:str=None, _l_markers:list=()):
        ### Main Variables ###
        self.condition = _condition # comma-separated string (ex. "AA_A_11_exon1,rs1234")
        self.condition_list = _condition_list # file path.

        self.l_condition = []
        self.N_condition = 0


        ### Main Actions ###
        # Condition itme check.
        if self.condition:
            self.l_condition = self.condition.split(',')
        else:
            with open(self.condition_list, 'r') as f_condlist:
                for line in f_condlist:
                    self.l_condition.append(line.rstrip('\n'))
        self.N_condition = len(self.l_condition)

        # All given conditions are in target marekrs?
        if len(_l_markers) > 0:
            l_excluded = list(np.setdiff1d(self.l_condition, _l_markers)) # supposed to [].
            if len(l_excluded) > 0:
                print(std_WARNING + "Some conditions are not in the target marker set: {}".format(l_excluded))


    def __bool__(self): return bool(self.condition) != bool(self.condition_list) # xor
    def isGivenAsFile(self): return (not self.condition) and self.condition_list

    def __repr__(self):
        # str_main = "< PLINK Condition file summary >\n"
        str_file = \
            "- Condition file: {}\n".format(self.condition_list) if self.isGivenAsFile() else \
            "- Condition string: {}\n".format(self.condition)
        str_conditions = \
            "- Conditions: {}\n".format(self.l_condition)
        str_N_conditions = \
            "- # of Conditions: {}\n".format(self.N_condition)

        str_summary = ''.join([str_file, str_conditions, str_N_conditions]).rstrip('\n')
        return str_summary



if __name__ == '__main__':

    ### Condition class
    # cond = CONDITION(_condition_list="/home/wansonchoi/sf_VirtualBox_Share/UC-CD-HLA/analysis/02-conditional-5/ToCond_CD.txt")
    # print(cond)

    cond = CONDITION(_condition="AA_DRB1_37_32660037_S,AA_DRB1_37_32660037_N")
    print(cond)
