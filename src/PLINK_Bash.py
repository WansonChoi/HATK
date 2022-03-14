#-*- coding: utf-8 -*-
import os, sys, re
from os.path import exists, join
import subprocess as sbp

from src.HATK_Error import HATK_PLINK_Execution_Error
from src.util import Exists

def Bash_RUN_PLINK(_command, _out_prefix, _f_save_log=False):
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
        if exists(_out_prefix+'.nosex'): os.remove(_out_prefix+'.nosex')
        if exists(_out_prefix+'.log') and not _f_save_log: os.remove(_out_prefix + '.log')
        return _out_prefix


