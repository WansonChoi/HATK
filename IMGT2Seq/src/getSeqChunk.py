#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, exists, isdir, join
import numpy as np

from src.util import printDict

"""
*_prot.txt => Prot
*_nuc.txt => cDNA
*_gen.txt => gDNA
"""



def getSeqChunk(_input):



    return 0



def getChunkOffsets(_input, _type="Prot"):

    l_offset = []

    with open(_input, 'r') as f_input:
        N_line = 0
        for line in f_input:
            if line.startswith(' {}'.format(_type)): l_offset.append(N_line)
            # if line.startswith('Please'): l_offset.append(N_line)
            N_line += 1

        l_offset.append(N_line) # The last line. (i.e. # of total line.)

    if len(l_offset) < 2: print("There could be only one chunk.")

    return l_offset



def getChunk(_input, _off_start, _off_end, _type):

    p_max_4field = re.compile(r'\w+\*\d{2,3}(:\d{2,3})*')

    with open(_input) as f_input:

        N_line = 0

        for line in f_input:

            # print("{}: {}".format(N_line, line))

            if N_line < _off_start:
                pass
            elif _off_start <= N_line and N_line < _off_end:

                l = line.split()
                # print(l)

                if len(l) >= 2: # ex. the last line is empty line.

                    m = p_max_4field.match(l[0])

                    if bool(m):
                        # __RETURN__[l[0]] = l[1:]
                        yield l[0], l[1:]

            else:
                break

            N_line += 1

    return 0


def getUniqHLA(_input, _type="Prot") -> list:

    l_offsets = getChunkOffsets(_input, _type)
    # print(l_offsets)

    # just the 1st chunk
    chunk_1st = getChunk(_input, l_offsets[0], l_offsets[1], _type)

    l_HLA_alleles = [hla_allele for hla_allele, _ in chunk_1st]

    __RETURN__ = list(np.unique([hla_allele.split('*')[0] for hla_allele in l_HLA_alleles]))
    # print(__RETURN__)

    return __RETURN__


if __name__ == '__main__':

    """
    
    """
    # input = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3470/alignments/DRB_prot.txt"
    input = "/home/wansonchoi/sf_VirtualBox_Share/IMGTHLA3470/alignments/DRB_nuc.txt"

    # r1 = getChunkOffsets(input) # [7, 3898, 7789, 11680, 11686]
    # print("Chunks: ", r1)

    # r2 = getChunk(input, 7, 3898, "Prot")
    # r2_sub = [l for n, l in enumerate(r2) if n < 10]
    # for _ in range(10): print(r2_sub[_])

    # r2 = getChunk(input, 11680, 11686, "Prot")
    # r2_sub = [l for n, l in enumerate(r2) if n < 10]
    # for item in r2_sub: print(item)

    ## *_nuc.txt
    r1 = getChunkOffsets(input, _type="cDNA")
    print(r1)
    r2 = getChunk(input, r1[0], r1[1], "cDNA")
    r2_sub = [l for n, l in enumerate(r2) if n < 10]
    for item in r2_sub: print(item)
    print(getUniqHLA(input, "cDNA"))

    # getUniqHLA(input)
