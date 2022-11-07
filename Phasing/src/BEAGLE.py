# -*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, exists, join

from src.util import Exists, checkFile
from src.HATK_Error import HATK_InputPreparation_Error, RaiseError

std_MAIN_PROCESS_NAME = "\n[%s]: " % basename(__file__)
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % basename(__file__)


class BEAGLE(object):

    def __init__(self, _bgl, _markers, _f_isPhased=False):
        ### Main Variables ###
        self.bgl = checkFile(_bgl, "Given Beagle file('{}') can't be found.".format(_bgl))
        self.markers = checkFile(_markers, "Given Beagle file('{}') can't be found.".format(_markers))

        self.f_isPhased = _f_isPhased

        ### Main Actions ###



    def GCtrick(self):
        pass

    def refineBP(self):
        pass

    # def __repr__(self):
    #     return 0


