#-*- coding: utf-8 -*-

import sys, os


def Bgl2GC(_bglfile, _markerfile, _bgloutfile, _markeroutfile):

    alleles={}
    with open(_markerfile) as mrk, open(_markeroutfile, 'w') as mrkout:
        for l in mrk:
            [ID, pos, A1, A2]=l.split()
            alleles[ID]=(A1, A2)
            newl=' '.join([ID, pos, "G", "C"])+'\n'
            mrkout.write(newl)

    with open(_bglfile) as bgl, open(_bgloutfile, 'w') as bglout:
        for l in bgl:
            c=l.split()
            if c[0] != 'M':
                bglout.write(l)
            else: ## If 'M', we replace alleles.
                ID=c[1]
                assert (ID in alleles), "Marker ID %s is not in the marker file"%ID
                (A1, A2)=alleles[ID]
                data=c[2:]
                newdata=["G" if x==A1 else ("C" if x==A2 else "0") for x in data]
                newl=' '.join(c[:2]+newdata)+'\n'
                bglout.write(newl)

    return [_bgloutfile, _markeroutfile]





if __name__ == '__main__':

    ## This code temporarilly replaces Allele1 and 2 with G and C, because
    ## beagle does not understand allele names such as P(presence) and A(absence)
    ## Usage: python [this file] [beagle file] [marker file] [beagle output file] [marker output file]
    ## v1.2: 2017-07-05. python porting. Buhm Han.

    [bglfile, markerfile, bgloutfile, markeroutfile] = sys.argv[1:]

    Bgl2GC(bglfile, markerfile, bgloutfile, markeroutfile)