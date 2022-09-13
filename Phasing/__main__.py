#-*- coding: utf-8 -*-

import os, sys, re
from os.path import basename, dirname, join, exists
from shutil import which
import argparse, textwrap

from Phasing.BEAGLE_Phase import Phasing_wrapper, VCF2BEAGLE
# from Phasing.src.BEAGLE import BEAGLE
from src.PLINK2BEAGLE import PLINK2BEAGLE

std_MAIN_PROCESS_NAME = "\n[Phasing]: "
std_ERROR_MAIN_PROCESS_NAME = "\n[Phasing::ERROR]: "
std_WARNING_MAIN_PROCESS_NAME = "\n[Phasing::WARNING]: "

class HATK_Phasing(object):

    def __init__(self, _out, _bfile=None, _vcf=None, _bgl=None, _markers=None,
                 _linkage2beagle="dependency/linkage2beagle.jar", _beagle2vcf="dependency/beagle2vcf.jar",
                 _vcf2beagle="dependency/vcf2beagle.jar", _beagle=which("beagle"), _plink=which("plink"),
                 _java=which("java"), _java_mem='1G', _nthreads=1, _f_save_intermediates=False, _f_gzip=True):

        """
        Input must be in beagle file format anyway.
            1. PLINK -> Beagle
            2. VCF -> Beagle
            3. Beagle

        """

        self.out_prefix = _out
        self.out_dir = dir(_out)

        self.bgl = None
        self.markers = None
        self.bgl_phased = None # main result.

        self.input_file_format = getInputFileFormat(_bfile, _vcf, _bgl, _markers)

        ### Setting BEAGLE input.
        if self.input_file_format == "PLINK":
            self.bgl, self.markers = PLINK2BEAGLE(_bfile, self.out_prefix +".P2B",
                                                  _linkage2beagle, _plink, _java, _java_mem, _f_save_intermediates)
        elif self.input_file_format== "VCF":
            # use 'VCF2BEAGLE'
            pass
        else:
            self.bgl = _bgl
            self.markers = _markers

        ### Phasing (Main)
        self.bgl_phased = \
            Phasing_wrapper(self.bgl, self.markers, self.out_prefix, _beagle2vcf, _vcf2beagle, _beagle,
                            _java, _java_mem, _nthreads, _f_save_intermediates, _f_gzip)


        print(self.markers)
        print(self.bgl_phased)


def getInputFileFormat(_bfile, _vcf, _beagle, _markers):

    if bool(_bfile):
        return "PLINK"
    elif bool(_vcf):
        return "VCF"
    else:
        return "BEAGLE"



if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='Phasing',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        [Phasing]

    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--bfile", help="\nThe prefix for PLINK genotype(SNP) data file(.bed/.bim/.fam)\n\n")

    parser.add_argument("--bgl", help="\nBeagle(v3.x.x) file.\n\n")
    parser.add_argument("--markers", help="\nThe '*.markers' file along with Beagle(v3.x.x) file.\n\n")

    # parser.add_argument("--vcf", help="\nThe VCF file.\n\n")

    parser.add_argument("--out", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')

    # java
    parser.add_argument("--java-mem", "-mem",
                        help="\nHeap memory size for each BEAGLE implementation (default: 2g).\n\n", default="2g")

    parser.add_argument("--nthreads", "-nth",
                        help="\nThe number of threads for each BEAGLE implementation (default: 1).\n\n", default=1,
                        type=int)


    ##### <for Test> #####

    # 2019. 01. 10
    # args = parser.parse_args(["--variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
    #                           "--out", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    HATK_Phasing(args.out, args.bfile, None, args.bgl, args.markers,
                 _f_save_intermediates=args.save_intermediates, _java_mem=args.java_mem, _nthreads=args.nthreads)
