#-*- coding: utf-8 -*-

import os, sys, re


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


class HLA_Study(object):


    def __init__(self, _args):

        """

        Main implementation body of HATK.
        HATK will work mainly by creating `HLA_Study` instance.

        """

        print(_args)
