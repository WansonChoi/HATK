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
    def __init__(self, _msg, _func=None):

        if _func:
            _func(_msg)
        else:
            print(_msg)
    # pass
