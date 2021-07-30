# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname

std_MAIN_PROCESS_NAME = "\n[%s]: " % (basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (basename(__file__))



class IMGT2SeqError(Exception):

    def __init__(self, _msg):
        print(_msg)



class ArgumentCheck(object):

    def __init__(self, f):
        self.f = f


    def __call__(self, _imgt, _hg, _out, *args, **kwargs):

        """
        Argument Check.
        """

        if not _imgt:
            raise IMGT2SeqError(
                std_ERROR_MAIN_PROCESS_NAME +
                "IMGT version information wasn't given. Please check the '--imgt' argument('{}') again.".format(_imgt)
            )


        if not _hg:
            raise IMGT2SeqError(
                std_ERROR_MAIN_PROCESS_NAME +
                "HG(Human Genome) version information wasn't given. Please check the '--hg' argument('{}') again.".format(_hg)
            )


        if not kwargs['_imgt_dir']:
            raise IMGT2SeqError(
                std_ERROR_MAIN_PROCESS_NAME +
                "Directory of the IMGT-HLA database wasn't given. Please check the '--imgt-dir' argument('{}') again.".format(kwargs['_imgt_dir'])
            )
        else:
            if not os.path.isdir(kwargs['_imgt_dir']):
                raise IMGT2SeqError(
                    std_ERROR_MAIN_PROCESS_NAME +
                    "Given directory for the IMGT-HLA database is not a directory. Please check the '--imgt-dir' argument('{}') again.".format(kwargs['_imgt_dir'])
                )


