#-*- coding: utf-8 -*-

class HATK_Error(Exception):
    """
    Base Error for HATK
    """
    pass

class HATK_InputPreparation_Error(HATK_Error):
    """
    Error related to input preparation
    """
    def __init__(self, _msg):
        pass

class HATK_PLINK_Execution_Error(HATK_Error):
    """
    Error related to PLINK bash execution.
    """
    def __init__(self, _msg):
        pass

class HATK_BEAGLE_Execution_Error(HATK_Error):
    """
    Error related to BEAGLE bash execution.
    """
    def __init__(self, _msg):
        pass

class HATK_R_Execution_Error(HATK_Error):
    """
    Error related to R(Rscript) execution.
    """
    def __init__(self, _msg):
        pass

def RaiseError(_Error, _msg):
    """
    To raise Error with Lambda expression.
    """
    raise _Error(_msg)