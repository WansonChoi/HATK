#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np

from src.HATK_Error import RaiseError, HATK_InputPreparation_Error
from src.util import Exists, checkFile, getColumn, FieldFormat2Label

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
    def __init__(self, _out_dir, _HLA_target=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

        self.out_dir = _out_dir

        self.hg = getPatterninFiles(_out_dir, re.compile(r'hg(\d+)'))
        self.imgt = getPatterninFiles(_out_dir, re.compile(r'imgt(\d+)'))
        self.HLA_target = list(_HLA_target)

        self.HAT = join(_out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(self.imgt)) \
            if Exists(join(_out_dir, "HLA_ALLELE_TABLE.imgt{}.hat".format(self.imgt))) else None

        self.HLA_DICTIONARY = HLA_DICTIONARY(
            join(_out_dir, "HLA_DICTIONARY_AA.hg{}.imgt{}.txt".format(self.hg, self.imgt)),
            join(_out_dir, "HLA_DICTIONARY_AA.hg{}.imgt{}.map".format(self.hg, self.imgt)),
            join(_out_dir, "HLA_DICTIONARY_SNPS.hg{}.imgt{}.txt".format(self.hg, self.imgt)),
            join(_out_dir, "HLA_DICTIONARY_SNPS.hg{}.imgt{}.map".format(self.hg, self.imgt)),
            _HLA_req=self.HLA_target)

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
        f_Maptable = all([Exists(mtable) for mtable in self.HLA_MAPTABLE.values()])
        # print(f_Maptable)

        return f_HAT and f_HLA_DICTIONARY and f_Maptable


    def __repr__(self):
        str_imgt_out_dir = \
            "- IMGT output directory: {}\n".format(self.out_dir)
        str_imgt_version = \
            "- IMGT version: {}\n".format("v{}".format(self.imgt) if bool(self.imgt) else None)
        str_hg = \
            "- Human Genome version: {}\n".format("hg{}".format(self.hg) if bool(self.hg) else None)

        str_HLA_Allele_Table = \
            "- HLA Allele Table:\n" \
            "   {}\n".format(self.HAT)

        str_HLA_DICTIONARY = "{}\n".format(self.HLA_DICTIONARY)

        str_maptable = \
            "- HLA Maptables for Heatmap plot: \n" \
            + ''.join([
                "   {hla}   : {maptable}\n".format(hla=hla, maptable=mtable) for hla, mtable in self.HLA_MAPTABLE.items()
            ])


        str_summary = ''.join([
            str_imgt_out_dir,
            str_imgt_version, str_hg, str_HLA_Allele_Table,
            str_HLA_DICTIONARY, str_maptable
        ]).rstrip('\n')

        return str_summary



class HLA_DICTIONARY(object):

    def __init__(self, _HLA_dict_AA_seq, _HLA_dict_AA_map, _HLA_dict_SNPS_seq, _HLA_dict_SNPS_map,
                 _HLA_req=("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")):

        self.HLA_req = _HLA_req

        ## AA dictionary
        self.HLA_dict_AA_seq = _HLA_dict_AA_seq if Exists(_HLA_dict_AA_seq) else None
        self.HLA_dict_AA_map = _HLA_dict_AA_map if Exists(_HLA_dict_AA_map) else None

        self.HLA_avail_AA_seq = getHLAs_seq(self.HLA_dict_AA_seq) if bool(self.HLA_dict_AA_seq) else []
        self.HLA_avail_AA_map = getHLAs_map(self.HLA_dict_AA_map) if bool(self.HLA_dict_AA_map) else []
        self.HLA_avail_AA = list(np.intersect1d(self.HLA_avail_AA_seq, self.HLA_avail_AA_map))

        self.field_format_AA = guessFieldFormat(self.HLA_dict_AA_seq) if bool(self.HLA_avail_AA_seq) else None

        ## SNPS dictionary
        self.HLA_dict_SNPS_seq = _HLA_dict_SNPS_seq if Exists(_HLA_dict_SNPS_seq) else None
        self.HLA_dict_SNPS_map = _HLA_dict_SNPS_map if Exists(_HLA_dict_SNPS_map) else None

        self.HLA_avail_SNPS_seq = getHLAs_seq(self.HLA_dict_SNPS_seq) if bool(self.HLA_dict_SNPS_seq) else []
        self.HLA_avail_SNPS_map = getHLAs_map(self.HLA_dict_SNPS_map) if bool(self.HLA_dict_SNPS_map) else []
        self.HLA_avail_SNPS = list(np.intersect1d(self.HLA_avail_SNPS_seq, self.HLA_avail_SNPS_map))

        self.field_format_SNPS = guessFieldFormat(self.HLA_dict_SNPS_seq) if bool(self.HLA_avail_SNPS_seq) else None


        # HLA available
        self.HLA_avail = list(np.intersect1d(self.HLA_avail_AA, self.HLA_avail_SNPS))
        self.f_all_avail = all(np.isin(_HLA_req, self.HLA_avail))

        self.HLA_target = list(np.intersect1d(self.HLA_req, self.HLA_avail)) # HLAs that can be used for bMG.

        self.HLA_excluded_from_req = list(np.setdiff1d(_HLA_req, self.HLA_target)) # => To exclude in CHPED
        self.HLA_excluded_from_avail = list(np.setdiff1d(self.HLA_avail, self.HLA_target)) # => To subset the map of HLA dictionary.

        # N-field information
        self.field_format = None
        self.f_field_format_match = self.field_format_AA == self.field_format_SNPS

        if self.f_field_format_match:
            self.field_format = self.field_format_AA
        else:
            self.field_format = "None"
            RaiseError(HATK_InputPreparation_Error,
                  std_ERROR +
                  "Each field format of Amino acid('{}') and Intragenic('{}') dictionaries doesn't match.\n" \
                  "- (AA)  : {}\n" \
                  "- (SNPS): {}".format(self.HLA_dict_AA_seq, self.HLA_dict_SNPS_seq,
                                        FieldFormat2Label(self.field_format_AA), FieldFormat2Label(self.field_format_SNPS))
            )


    def __bool__(self):
        f_HLA_dict_AA = Exists(self.HLA_dict_AA_seq) and Exists(self.HLA_dict_AA_map)
        f_HLA_dict_SNPS = Exists(self.HLA_dict_SNPS_seq) and Exists(self.HLA_dict_SNPS_map)

        if not Exists(self.HLA_dict_AA_seq):
            print(std_WARNING + "Given Amino acid sequence dictionary file('{}') can't be found.".format(self.HLA_dict_AA_seq))
        if not Exists(self.HLA_dict_AA_map):
            print(std_WARNING + "Given Amino acid sequence dictionary map file('{}') can't be found.".format(self.HLA_dict_AA_map))
        if not Exists(self.HLA_dict_SNPS_seq):
            print(std_WARNING + "Given Intragenic DNA sequence dictionary file('{}') can't be found.".format(self.HLA_dict_SNPS_seq))
        if not Exists(self.HLA_dict_SNPS_map):
            print(std_WARNING + "Given Intragenic DNA sequence dictionary map file('{}') can't be found.".format(self.HLA_dict_SNPS_map))

        return f_HLA_dict_AA and f_HLA_dict_SNPS and self.f_all_avail and self.f_field_format_match


    def __repr__(self):
        str_HLA = \
            "- HLA:\n" \
            "   (requested): {}\n\n" \
            "   (available, AA, seq): {}\n" \
            "   (available, AA, map): {}\n" \
            "   (available, AA): {}\n" \
            "   (available, SNPS, seq): {}\n" \
            "   (available, SNPS, map): {}\n" \
            "   (available, SNPS): {}\n" \
            "   (available, Both): {}\n\n" \
            "   (Are all the requested genes in the available?): {}\n" \
            "   (excluded from requested): {}\n" \
            "   (excluded from available): {}\n\n" \
            "   (target): {}\n" \
                .format(list(self.HLA_req),
                        self.HLA_avail_AA_seq, self.HLA_avail_AA_map, self.HLA_avail_AA,
                        self.HLA_avail_SNPS_seq, self.HLA_avail_SNPS_map, self.HLA_avail_SNPS,
                        self.HLA_avail,
                        self.f_all_avail, self.HLA_excluded_from_req, self.HLA_excluded_from_avail,
                        self.HLA_target)

        str_dict_AA = \
            "- Dictionary(Amino acids)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.HLA_dict_AA_seq, self.HLA_dict_AA_map)

        str_dict_SNPS = \
            "- Dictionary(Intragenic SNPs)\n" \
            "   (seq): {}\n" \
            "   (map): {}\n" \
            .format(self.HLA_dict_SNPS_seq, self.HLA_dict_SNPS_map)

        str_field_format = \
            "- Field format: {}\n".format(FieldFormat2Label(self.field_format))

        str_summary = ''.join([
            str_HLA, str_dict_AA, str_dict_SNPS, str_field_format
        ]).rstrip('\n')

        return str_summary


    def subsetHLA_DICTIONARY(self, _out_dir, _HLA_ToExtract:list, _HLA_ToExclude=()):

        _HLA_ToExtract = list(np.setdiff1d(_HLA_ToExtract, _HLA_ToExclude))

        dict_AA_seq_subset = subsetSeqDict(self.HLA_dict_AA_seq, join(_out_dir, basename(self.HLA_dict_AA_seq)+'.subset'),
                                           _HLA_ToExtract)
        dict_AA_map_subset = subsetDictMap(self.HLA_dict_AA_map, join(_out_dir, basename(self.HLA_dict_AA_map)+'.subset'),
                                           _HLA_ToExtract)

        dict_SNPS_seq_subset = subsetSeqDict(self.HLA_dict_SNPS_seq, join(_out_dir, basename(self.HLA_dict_SNPS_seq)+'.subset'),
                                             _HLA_ToExtract)
        dict_SNPS_map_subset = subsetDictMap(self.HLA_dict_SNPS_map, join(_out_dir, basename(self.HLA_dict_SNPS_map)+'.subset'),
                                             _HLA_ToExtract)

        return HLA_DICTIONARY(
            dict_AA_seq_subset, dict_AA_map_subset, dict_SNPS_seq_subset, dict_SNPS_map_subset,
            _HLA_ToExtract
        )



def getHLAs_seq(_HLA_dict_seq) -> list:

    l_temp_HLA = [HLA_allele.split('*')[0] for HLA_allele in getColumn(_HLA_dict_seq, 0)]

    return list(np.unique(l_temp_HLA))



def getHLAs_map(_HLA_dict_map) -> list:

    func = lambda x: 2 if x.startswith('INS') else 1
    """
    - AA_'A'_-6_30018365_exon1, SNPS_'A'_73_30018382_exon1 => idx: 1
    - INS_AA_'A'_-6x-5_30018366, INS_SNPS_'A'_73x74_30018382 => idx: 2
    """

    l_temp_HLA = [label.split('_')[func(label)] for label in getColumn(_HLA_dict_map, 1)]

    return list(np.unique(l_temp_HLA))



def getPatterninFiles(_out_dir, _p):

    l_files = os.listdir(_out_dir)

    for file in l_files:
        s = _p.search(file)
        if bool(s):
            return s.group(1)



def guessFieldFormat(_dict_seq, _N_sample=-1, _seed=0):

    l_allele = [HLA_allele.split('*')[1] for HLA_allele in getColumn(_dict_seq, 0)]

    if _N_sample > 0:
        np.random.seed(_seed)
        _N_sample = len(l_allele) if _N_sample > len(l_allele) else _N_sample
        l_allele = np.random.choice(l_allele, _N_sample)

        # print(l_allele[:20])
        # print(len(l_allele))

    ## Main iteration
    Nfield_max = -1
    for allele in l_allele:
        # G-group
        if allele.endswith('G'):
            return 5
        # P-group
        if allele.endswith('P'):
            return 6

        temp_Nfield = len(allele.split(':'))
        if temp_Nfield > Nfield_max:
            Nfield_max = temp_Nfield

    return Nfield_max



def subsetSeqDict(_HLA_dict_seq, _out, _HLA_ToExtract):

    with open(_HLA_dict_seq, 'r') as f_HLA_dict_seq, open(_out, 'w') as f_out:
        for line in f_HLA_dict_seq:
            l = line.split()
            hla = l[0].split('*')[0]

            if hla in _HLA_ToExtract:
                f_out.write(line)

    return _out

def subsetDictMap(_HLA_dict_map, _out, _HLA_ToExtract):

    with open(_HLA_dict_map, 'r') as f_HLA_dict_seq, open(_out, 'w') as f_out:
        for line in f_HLA_dict_seq:
            l = line.split()

            temp = l[1].split('_')

            type = temp[0]

            if type == "AA" or type == "SNPS":
                hla = temp[1] # "AA_A_-6_30018365_exon1", "SNPS_A_73_30018382_exon1"
            else:
                hla = temp[2] # "INS_AA_A_-6x-5_30018366", "INS_SNPS_A_73x74_30018382"

            if hla in _HLA_ToExtract:
                f_out.write(line)

    return _out



if __name__ == '__main__':

    # HLA_Dictionary_AA_prefix = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_IMGT2Seq_20220411/HLA_DICTIONARY_AA.hg18.imgt3470"
    # HLA_Dictionary_SNPS_prefix = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_IMGT2Seq_20220411/HLA_DICTIONARY_SNPS.hg18.imgt3470"

    # r_AA = getHLAs(HLA_Dictionary_AA_prefix+'.txt')
    # print(r_AA)
    # r_SNPS = getHLAs(HLA_Dictionary_SNPS_prefix+'.txt')
    # print(r_SNPS)

    # hla_dict = HLA_DICTIONARY(HLA_Dictionary_AA_prefix+'.txt', HLA_Dictionary_AA_prefix+'.map',
    #                           HLA_Dictionary_SNPS_prefix+'.txt', HLA_Dictionary_SNPS_prefix+'.map')
    # print(hla_dict)
    #
    # hla_dict2 = \
    #     hla_dict.subsetHLA_DICTIONARY('/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_IMGT2Seq_20220411',
    #                                   ['A', 'B', 'C'], ['DRB1'])
    # print(hla_dict2)

    # out_dir = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_IMGT2Seq_20220216"
    # r = IMGT2Seq_Output(out_dir)
    # print(r)
    # print(bool(r))

    # Fieldformat_guess = guessFieldFormat(HLA_Dictionary_AA_seq)
    # print(Fieldformat_guess)


    # r = subsetSeqDict(HLA_Dictionary_AA_prefix+'.txt', HLA_Dictionary_AA_prefix+'.txt.subset', ['A', 'B', 'C'])
    # print(r)
    # r = subsetDictMap(HLA_Dictionary_AA_prefix+'.map', HLA_Dictionary_AA_prefix+'.map.subset', ['A', 'B', 'C'])
    # print(r)
    # r = subsetSeqDict(HLA_Dictionary_SNPS_prefix+'.txt', HLA_Dictionary_SNPS_prefix+'.txt.subset', ['A', 'B', 'C'])
    # print(r)
    # r = subsetDictMap(HLA_Dictionary_SNPS_prefix+'.map', HLA_Dictionary_SNPS_prefix+'.map.subset', ['A', 'B', 'C'])
    # print(r)

    ###########
    # (2022.10.11.)
    out_dir = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/HATK_rearchitect_IMGT2Seq_20220913_2field"
    imgt_output = IMGT2Seq_Output(out_dir, ['A', 'B', 'C', 'DMA', 'DMB', 'DOA', 'DOB', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'E', 'F', 'G', 'MICB'])
    print(imgt_output)

    pass