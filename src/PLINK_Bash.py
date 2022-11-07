#-*- coding: utf-8 -*-
import os, sys, re
from os.path import exists, join
import subprocess as sbp
from shutil import which

from src.HATK_Error import HATK_PLINK_Execution_Error, HATK_InputPreparation_Error, RaiseError
from src.util import Exists


def getPLINKexec():

    if exists(which("plink")):
        return which("plink")
    elif exists("./dependency/plink"):
        return "./dependency/plink"
    else:
        RaiseError(
            HATK_InputPreparation_Error,
            "PLINK(v1.9b) binary executable file can't be found. "
            "Please check 'dependency/' folder or Install PLINK in your system manually."
        )



def Bash_RUN_PLINK(_command, _out_prefix, _f_save_intermediates=False, _f_save_log=False):
    try:
        # print(_command)
        # print(_command.split())
        sbp.run(_command.split(), check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)

    except sbp.CalledProcessError:
        raise HATK_PLINK_Execution_Error(
            "Next PLINK execution failed. ('{}')".format(_command) +
            (" Please refer to its log('{}') file.".format(_out_prefix + '.log') if Exists(_out_prefix + '.log') else "")
        )
    else:
        if not _f_save_intermediates:
            if exists(_out_prefix+'.nosex'): os.remove(_out_prefix+'.nosex')
            if exists(_out_prefix+'.log') and not _f_save_log: os.remove(_out_prefix + '.log')
        return _out_prefix


