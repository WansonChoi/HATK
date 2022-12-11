# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

from src.util import printDict
import IMGT2Seq.src.processModules as modules


def processNuc(_nuc):

    """
    - Nuc file is just used for generating AA map.
    - Only ref seq is needed.

    (2022.12.11.)
    planned to jointly use the nuc's ref seq. in generating map file, but I postponed it again.
    Using nuc's ref seq to generate map file with AA and SNPS might cause more trouble.

    For now, likewise previous version, the relative position start information will be acquired. (ex. HLA-A: -24)

    """

    ### Load gen raw seqs
    dict_chunks_gen = modules.chopAlignmentFile(_nuc)
    # printDict(dict_chunks_gen)

    """

    """

    temp_chunk = dict_chunks_gen[1] # 1st chunk

    # line1 = temp_chunk[0].split()
    # print(line1)
    line2 = temp_chunk[1].split()
    # print(line2)

    AA_rel_pos_start = line2[2] # ex. '-32'

    return AA_rel_pos_start



if __name__ == '__main__':

    nuc = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3500/alignments/DQB1_nuc.txt"

    processNuc(nuc)

    pass