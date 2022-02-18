#-*- coding: utf-8 -*-

from os.path import basename, dirname, exists, isdir, join

"""
Collection of small utility functions
"""

Exists = lambda x, func=exists: x and func(x)


def printDF(_df_name, _df):
    print("{}:\n{}".format(_df_name, _df))


def printDict(_dict):
    for k, v in _dict.items():
        print("{}: {}".format(k, v))