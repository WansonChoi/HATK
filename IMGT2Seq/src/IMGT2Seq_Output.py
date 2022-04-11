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
    def __init__(self, _out_dir, _hg, _imgt):

        self.hg = _hg
        self.imgt = _imgt

        self.HAT = join(_out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(_imgt)) \
            if Exists(join(_out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(_imgt))) else None

        self.HLA_DICTIONARY = HLA_Dictionary(
            join(_out_dir, "HLA_DICTIONARY_AA.hg{}.imgt{}.txt".format(_hg, _imgt)),
            join(_out_dir, "HLA_DICTIONARY_AA.hg{}.imgt{}.map".format(_hg, _imgt)),
            join(_out_dir, "HLA_DICTIONARY_SNPS.hg{}.imgt{}.txt".format(_hg, _imgt)),
            join(_out_dir, "HLA_DICTIONARY_SNPS.hg{}.imgt{}.map".format(_hg, _imgt))
        )

        self.HLA_MAPTABLE = {
            hla: join(_out_dir, "HLA_MAPTABLE_{}.hg{}.imgt{}.txt".format(hla, _hg, _imgt)) \
                if Exists(join(_out_dir, "HLA_MAPTABLE_{}.hg{}.imgt{}.txt".format(hla, _hg, _imgt))) else None
            for hla in self.HLA_DICTIONARY.HLA_avail
        }


    def __bool__(self):
        # print(Exists(self.HAT))
        # print(self.HLA_DICTIONARY.__bool__())
        # print(np.all([Exists(mtable) for mtable in self.HLA_MAPTABLE.values()]))
        # print(Exists(self.HAT) and self.HLA_DICTIONARY.__bool__() and \
        #        np.all([Exists(mtable) for mtable in self.HLA_MAPTABLE.values()]))
        # for mtable in self.HLA_MAPTABLE.values():
        #     print(mtable)
        #     print(Exists(mtable))

        return Exists(self.HAT) and self.HLA_DICTIONARY.__bool__() and \
               bool(np.all([Exists(mtable) for mtable in self.HLA_MAPTABLE.values()]))
        # bool() -> TypeError: __bool__ should return bool, returned numpy.bool_


    def __repr__(self):
        str_imgt_version = \
            "- IMGT version: v{}\n".format(self.imgt)
        str_hg = \
            "- Human Genome version: hg{}\n".format(self.hg)
        str_HLA = \
            "- HLA:\n" \
            "   (available,   AA): {}\n" \
            "   (available, SNPS): {}\n"\
                .format(self.HLA_DICTIONARY.HLA_avail_AA, self.HLA_DICTIONARY.HLA_avail_SNPS)


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




class HLA_Dictionary(object):

    def __init__(self, _HLA_dict_AA_seq, _HLA_dict_AA_map, _HLA_dict_SNPS_seq, _HLA_dict_SNPS_map):

        self.HLA_dict_AA_seq = checkFile(_HLA_dict_AA_seq, std_ERROR+"Given HLA Amino acid sequence dictionary('{}') can't be found.".format(_HLA_dict_AA_seq))
        self.HLA_dict_AA_map = checkFile(_HLA_dict_AA_map, std_ERROR+"Given HLA Amino acid sequence map dictionary('{}') can't be found.".format(_HLA_dict_AA_map))
        self.HLA_avail_AA = getHLAs(self.HLA_dict_AA_seq)

        self.HLA_dict_SNPS_seq = checkFile(_HLA_dict_SNPS_seq, std_ERROR+"Given HLA Intragenic DNA sequence dictionary('{}') can't be found.".format(_HLA_dict_SNPS_seq))
        self.HLA_dict_SNPS_map = checkFile(_HLA_dict_SNPS_map, std_ERROR+"Given HLA Intragenic DNA sequence map dictionary('{}') can't be found.".format(_HLA_dict_SNPS_map))
        self.HLA_avail_SNPS = getHLAs(self.HLA_dict_SNPS_seq)

        self.HLA_avail = list(np.intersect1d(self.HLA_avail_AA, self.HLA_avail_SNPS))


    def __bool__(self):
        return Exists(self.HLA_dict_AA_seq) and Exists(self.HLA_dict_AA_map) and \
                Exists(self.HLA_dict_SNPS_seq) and Exists(self.HLA_dict_SNPS_map)


def getHLAs(_HLA_dict_seq) -> list:

    l_temp_HLA = [HLA_allele.split('*')[0] for HLA_allele in getColumn(_HLA_dict_seq, 0)]

    return list(np.unique(l_temp_HLA))



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

    r = IMGT2Seq_Output("/Users/wansonchoi/Git_Projects/HATK/tests/20220409_IMGT2Seq_output", 18, 3470)
    print(r)
    print(bool(r))


    pass