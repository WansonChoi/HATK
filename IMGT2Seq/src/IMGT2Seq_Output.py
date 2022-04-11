#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np

from src.HATK_Error import RaiseError, HATK_InputPreparation_Error
from src.util import Exists, checkFile, getColumn

std_MAIN = "\n[%s]: " % basename(__file__)
std_ERROR = "\n[%s::ERROR]: " % basename(__file__)
std_WARNING = "\n[%s::WARNING]: " % basename(__file__)


class IMGT2Seq_Output(object):
    """
    - Wrapper class for output of the IMGT2Seq.

    1. HAT.
    2. HLA Dictionary
    3. Maptable
    """
    def __init__(self, _out_dir, _HLA_req=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

        self.hg = getPatterninFiles(_out_dir, re.compile(r'hg(\d+)'))
        self.imgt = getPatterninFiles(_out_dir, re.compile(r'imgt(\d+)'))
        self.HLA_req = list(_HLA_req)

        self.HAT = join(_out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(self.imgt)) \
            if Exists(join(_out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(self.imgt))) else None

        self.HLA_DICTIONARY = HLA_DICTIONARY(
            join(_out_dir, "HLA_DICTIONARY_AA.hg{}.imgt{}.txt".format(self.hg, self.imgt)),
            join(_out_dir, "HLA_DICTIONARY_AA.hg{}.imgt{}.map".format(self.hg, self.imgt)),
            join(_out_dir, "HLA_DICTIONARY_SNPS.hg{}.imgt{}.txt".format(self.hg, self.imgt)),
            join(_out_dir, "HLA_DICTIONARY_SNPS.hg{}.imgt{}.map".format(self.hg, self.imgt)),
            _HLA_req=self.HLA_req
        )

        self.HLA_MAPTABLE = {
            hla: join(_out_dir, "HLA_MAPTABLE_{}.hg{}.imgt{}.txt".format(hla, self.hg, self.imgt)) \
                if Exists(join(_out_dir, "HLA_MAPTABLE_{}.hg{}.imgt{}.txt".format(hla, self.hg, self.imgt))) else None
            for hla in self.HLA_DICTIONARY.HLA_avail
        }


    def __bool__(self):

        f_HAT = Exists(self.HAT)
        # print(f_HAT)
        f_HLA_DICTIONARY = bool(self.HLA_DICTIONARY)
        # print(f_HLA_DICTIONARY)
        f_Maptable = bool(np.all([Exists(mtable) for mtable in self.HLA_MAPTABLE.values()])) # bool() -> TypeError: __bool__ should return bool, returned numpy.bool_
        # print(f_Maptable)

        return f_HAT and f_HLA_DICTIONARY and f_Maptable



    def __repr__(self):
        str_imgt_version = \
            "- IMGT version: {}\n".format("v{}".format(self.imgt) if bool(self.imgt) else None)
        str_hg = \
            "- Human Genome version: {}\n".format("hg{}".format(self.hg) if bool(self.hg) else None)
        str_HLA = \
            "- HLA:\n" \
            "   (requested): {}\n" \
            "   (available,   AA): {}\n" \
            "   (available, SNPS): {}\n" \
            "   (available, Both): {}\n" \
            "   (Is all the requested are in the available?): {}\n" \
            "   (Not available): {}\n" \
                .format(self.HLA_req,
                        self.HLA_DICTIONARY.HLA_avail_AA, self.HLA_DICTIONARY.HLA_avail_SNPS, self.HLA_DICTIONARY.HLA_avail,
                        self.HLA_DICTIONARY.f_all_avail, self.HLA_DICTIONARY.HLA_Not_avail)

        str_HLA_Allele_Table = \
            "- HLA Allele Table:\n" \
            "   {}\n".format(self.HAT)

        str_dict_AA = \
            "- Dictionary(Amino acids)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.HLA_DICTIONARY.HLA_dict_AA_seq, self.HLA_DICTIONARY.HLA_dict_AA_map)

        str_dict_SNPS = \
            "- Dictionary(Intragenic SNPs)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.HLA_DICTIONARY.HLA_dict_SNPS_seq, self.HLA_DICTIONARY.HLA_dict_SNPS_map)

        str_maptable = \
            "- HLA Maptables for Heatmap plot: \n" \
            + ''.join([
                "   {hla}   : {maptable}\n".format(hla=hla, maptable=mtable) for hla, mtable in self.HLA_MAPTABLE.items()
            ])


        str_summary = ''.join([
            str_imgt_version, str_hg, str_HLA,
            str_HLA_Allele_Table, str_dict_AA, str_dict_SNPS,
            str_maptable
        ]).rstrip('\n')

        return str_summary



class HLA_DICTIONARY(object):

    def __init__(self, _HLA_dict_AA_seq, _HLA_dict_AA_map, _HLA_dict_SNPS_seq, _HLA_dict_SNPS_map,
                 _HLA_req=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

        self.HLA_dict_AA_seq = _HLA_dict_AA_seq if Exists(_HLA_dict_AA_seq) else None
        self.HLA_dict_AA_map = _HLA_dict_AA_map if Exists(_HLA_dict_AA_map) else None
        self.HLA_avail_AA = getHLAs(self.HLA_dict_AA_seq) if bool(self.HLA_dict_AA_seq) else []

        self.HLA_dict_SNPS_seq = _HLA_dict_SNPS_seq if Exists(_HLA_dict_SNPS_seq) else None
        self.HLA_dict_SNPS_map = _HLA_dict_SNPS_map if Exists(_HLA_dict_SNPS_map) else None
        self.HLA_avail_SNPS = getHLAs(self.HLA_dict_SNPS_seq) if bool(self.HLA_dict_SNPS_seq) else []

        self.HLA_avail = list(np.intersect1d(self.HLA_avail_AA, self.HLA_avail_SNPS))
        self.HLA_Not_avail = list(np.setdiff1d(_HLA_req, self.HLA_avail))
        self.f_all_avail = all(np.isin(_HLA_req, self.HLA_avail))


    def __bool__(self):
        f_HLA_dict_AA = Exists(self.HLA_dict_AA_seq) and Exists(self.HLA_dict_AA_map)
        f_HLA_dict_SNPS = Exists(self.HLA_dict_SNPS_seq) and Exists(self.HLA_dict_SNPS_map)

        return f_HLA_dict_AA and f_HLA_dict_SNPS and self.f_all_avail



def getHLAs(_HLA_dict_seq) -> list:

    l_temp_HLA = [HLA_allele.split('*')[0] for HLA_allele in getColumn(_HLA_dict_seq, 0)]

    return list(np.unique(l_temp_HLA))



def getPatterninFiles(_out_dir, _p):

    l_files = os.listdir(_out_dir)

    for file in l_files:
        s = _p.search(file)
        if bool(s):
            return s.group(1)



if __name__ == '__main__':

    HLA_Dictionary_AA_prefix = "/Users/wansonchoi/Git_Projects/HATK/tests/20220409_IMGT2Seq_output/HLA_DICTIONARY_AA.hg18.imgt3470"
    HLA_Dictionary_SNPS_prefix = "/Users/wansonchoi/Git_Projects/HATK/tests/20220409_IMGT2Seq_output/HLA_DICTIONARY_SNPS.hg18.imgt3470"

    # r_AA = getHLAs(HLA_Dictionary_AA_prefix+'.txt')
    # print(r_AA)
    # r_SNPS = getHLAs(HLA_Dictionary_SNPS_prefix+'.txt')
    # print(r_SNPS)

    # hla_dict = HLA_Dictionary(HLA_Dictionary_AA_prefix+'.txt', HLA_Dictionary_AA_prefix+'.map',
    #                           HLA_Dictionary_SNPS_prefix+'.txt', HLA_Dictionary_SNPS_prefix+'.map')
    # print(hla_dict.HLA_avail_AA)
    # print(hla_dict.HLA_avail_SNPS)
    # print(hla_dict.HLA_avail)


    out_dir = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_IMGT2Seq_20220216"
    r = IMGT2Seq_Output(out_dir)
    print(r)
    print(bool(r))


    pass