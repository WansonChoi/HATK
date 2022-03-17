#-*- coding: utf-8 -*-

import os
from shutil import which

from src.PLINK_Bash import Bash_RUN_PLINK
from Phasing.src.BEAGLE_Bash import Bash_BEAGLE

def PLINK2BEAGLE(_bfile, _out_prefix, plink, _linkage2beagle, _java=which("java"), _java_mem='1G',
                 _f_save_intermediates=False):

    ### Main Variables ###
    _bim = _bfile+'.bim'

    _map = None
    _ped = None

    _dat = None
    _ped_nopheno = None

    _markers = _out_prefix+'.markers'
    _bgl = _out_prefix+'.bgl'


    ### Main Actions ###

    # *.markers
    command = ' '.join(["awk", '\'{print $2 " " $4 " " $5 " " $6}\'', _bim, ">", _markers])
    # print(command)
    os.system(command)

    # *.{dat,nopheno.ped}
    command = "{plink} --recode --bfile {bfile} --keep-allele-order --alleleACGT --out {out}" \
                .format(plink=plink, bfile=_bfile, out=_out_prefix)
    # print(command)
    Bash_RUN_PLINK(command, _out_prefix)

    _ped = _out_prefix+'.ped'
    _map = _out_prefix+'.map'

    command = ' '.join(["awk", '\'{print "M " $2}\'', _map, ">", _out_prefix+'.dat'])
    # print(command)
    os.system(command)

    _dat = _out_prefix+'.dat'

    command = ' '.join(["cut -d ' ' -f1-5,7-", _ped, ">", _out_prefix+'.nopheno.ped'])
    # print(command)
    os.system(command)

    _ped_nopheno = _out_prefix+'.nopheno.ped'



    # Main Conversion (linkage2beagle.jar)
    command = "{java} -Xmx{java_mem} -jar {l2b} {data} {ped_nopheno}" \
                .format(java=_java, java_mem=_java_mem, l2b=_linkage2beagle, data=_dat, ped_nopheno=_ped_nopheno)
    # print(command)
    Bash_BEAGLE(command, _bgl)



    if not _f_save_intermediates:
        os.remove(_ped)
        os.remove(_map)

        os.remove(_ped_nopheno)
        os.remove(_dat)


    return _bgl, _markers



if __name__ == '__main__':

    """
    PLINK2BEAGLE.py
    : Conversion from PLINK genotype to BEAGLE(v3.x) file format.
    
    """

    bfile = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset"
    out_prefix = "/home/wansonchoi/sf_VirtualBox_Share/HATK/tests/20220314_BEAGLE/wtccc_filtered_58C_RA.hatk.300+300.imgt3470.header.subset.tt"
    plink = "/home/wansonchoi/miniconda3/envs/HATK/bin/plink"
    linkage2beagle = "/media/sf_VirtualBox_Share/HATK/dependency/linkage2beagle.jar"

    PLINK2BEAGLE(bfile, out_prefix, plink, linkage2beagle)