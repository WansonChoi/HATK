#-*- coding: utf-8 -*-
import os, sys, re
from os.path import exists, join
import subprocess as sbp

from src.HATK_Error import HATK_PLINK_Execution_Error

def Bash_RUN_PLINK(_command, _out_prefix, _f_remove_log=True):
    try:
        # print(_command)
        # print(_command.split())
        sbp.run(_command.split(), check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)

    except sbp.CalledProcessError:
        raise HATK_PLINK_Execution_Error(
            "PLINK execution for the command('{}') failed. "
            "Please refer to its log('{}') file.".format(_command, _out_prefix + '.log')
        )
    else:
        if exists(_out_prefix+'.nosex'): os.remove(_out_prefix+'.nosex')
        if exists(_out_prefix+'.log') and _f_remove_log: os.remove(_out_prefix+'.log')
        return _out_prefix


