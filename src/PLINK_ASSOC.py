#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join
import pandas as pd

from src.util import Exists

class ASSOC(object):

    def __init__(self, _assoc, _regr_type, _out=None, _top_N=5):
        """
        A container class for association test(linear/logistic regression) result.
        """

        ### Main Variables ###
        self.assoc = _assoc
        self.regr_type = _regr_type # must be either "linear" or "logistic".
        self.top_N = _top_N

        self.assoc_sort, \
        self.df_assoc_topN = sort_assoc(self.assoc, _out if _out else self.assoc+'.sort', self.top_N)
        self.l_topN = self.df_assoc_topN.iloc[:, 1].tolist()



    def __repr__(self):

        str_assoc = \
            "- Assoc file: {}\n".format(self.assoc)
        str_regr_type = \
            "- Regression type: {}\n".format(self.regr_type)
        str_assoc_sort = \
            "- Sorted Assoc file: {}\n".format(self.assoc_sort)
        str_top_N = \
            "- Top {} SNPs:\n{}\n".format(self.top_N, self.df_assoc_topN)

        str_summary = ''.join([
            str_assoc, str_regr_type, str_assoc_sort, str_top_N
        ])

        return str_summary


    def __bool__(self): return Exists(self.assoc)


def sort_assoc(_assoc, _out, _top_N=5):

    df_assoc = pd.read_csv(_assoc, sep='\s+', header=0)
    # print("df_assoc:\n{}\n".format(df_assoc))

    """
     CHR                                  SNP         BP   A1       TEST    NMISS         OR       SE      L95      U95         STAT            P 
       6                            rs2894066   29001906    G        ADD     3128     0.8604  0.04966   0.7806   0.9484       -3.028     0.002464
       6                            rs3763338   29002290    T        ADD     3128     0.6916  0.05478   0.6212   0.7699       -6.733    1.661e-11
       6                            rs1237485   29002323    C        ADD     3128     0.8228  0.05062   0.7451   0.9087       -3.852    0.0001171
       6                            rs2015436   29004893    T        ADD     3128     0.8922  0.04943   0.8098   0.9829       -2.308      0.02097
       6                             rs932776   29009013    A        ADD     3128     0.8228  0.05062   0.7451   0.9087       -3.852    0.0001171
   """

    df_assoc_sort = df_assoc.sort_values('P')
    # print("df_assoc_sort:\n{}\n".format(df_assoc_sort))

    df_assoc_sort_topN = df_assoc_sort.iloc[:_top_N, :]
    # print("df_assoc_sort_topN:\n{}\n".format(df_assoc_sort_topN))


    # write sorted assoc.
    df_assoc_sort.to_csv(_out, sep='\t', header=True, index=False)

    return _out, df_assoc_sort_topN



if __name__ == '__main__':

    ## LINUX

    ## Mac
    _assoc = "/Users/wansonchoi/Git_Projects/HATK/tests/20221007_bMG/asdf.assoc.logistic"

    # sort_assoc(_assoc, _assoc+'.sort')
    r = ASSOC(_assoc, "logistic")
    print(r)