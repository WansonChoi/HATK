#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, join
import subprocess as sbp

from src.HATK_Error import HATK_BEAGLE_Execution_Error
from src.util import Exists

def Bash_BEAGLE(_command, _out, _f_save_log=False):
    try:
        with open(_out, 'w') as f_out, open(_out + '.log', 'w') as f_out_log:
            sbp.run(_command.split(), check=True, stdout=f_out, stderr=f_out_log)

    except sbp.CalledProcessError:
        raise HATK_BEAGLE_Execution_Error(
            "Next BEAGLE execution failed. ('{}')".format(_command) + \
            (" Please refer to its log('{}') file.".format(_out + '.log') if Exists(_out + '.log') else "")
        )
    else:
        if Exists(_out + '.log') and not _f_save_log: os.remove(_out + '.log')
        return _out



def Bash_BEAGLE_Phase(_command, _out, _f_save_log=True):
    try:
        sbp.run(_command.split(), check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)

    except sbp.CalledProcessError:
        raise HATK_BEAGLE_Execution_Error(
            "Next BEAGLE Phasing failed. ('{}')".format(_command) + \
            (" Please refer to its log('{}') file.".format(_out + '.log') if Exists(_out + '.log') else "")
        )
    else:
        if Exists(_out + '.log') and not _f_save_log: os.remove(_out + '.log')
        return _out