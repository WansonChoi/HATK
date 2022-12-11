# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

from src.util import printDict
import IMGT2Seq.src.processModules as modules

def processProt(_prot):

    print("\n\n< Load raw seq. >")
    dict_chunks = modules.chopAlignmentFile(_prot)

    print("\n\n< Join seqs within chunks >")
    l_joined_within = modules.join_within_chunks_Prot(dict_chunks)
    for i in range(len(l_joined_within)):
        print("{}th chunk:".format(i))
        printDict(l_joined_within[i], 5)

    print("\n\n< Join seqs between chunks >")
    dict_joined_between = modules.join_between_chunks_Prot(l_joined_within)
    printDict(dict_joined_between, 10)
    # print(dict_joined_between['A*03:437Q'])
    # print("{} vs. {}".format(len(dict_joined_between['A*01:01:01:01']), len(dict_joined_between['A*03:437Q'])))

    print("\n\n< Ligate the right end. >")
    dict_ligated = modules.ligateRightEnd(dict_joined_between)
    printDict(dict_ligated, 10)

    print("\n\n< Substitute bases. >")
    dict_substituted = modules.substituteBase(dict_ligated)
    printDict(dict_substituted, 10)

    print("\n\n< Process Indel. >")
    dict_IndexProcessed = modules.processIndel(dict_substituted)
    printDict(dict_IndexProcessed, 10)

    return 0



if __name__ == '__main__':

    # A, prot
    # prot = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3480/alignments/A_prot.txt"
    prot = "/Users/wansonchoi/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3500/alignments/A_prot.txt"

    processProt(prot)

    """
    (Prot)
    - Load raw sequences.
    - Join seqs within chunks.
    - Join seqs between chunks.

    - Ligate right end.
    - Substitute '-'(dash), '*', and 'x'.
    - Process indel.

    - Extract the ref. seq for MAP file generation.


    """

    pass