#-*- coding: utf-8 -*-

import os, sys, re
import shutil
from os.path import basename, dirname, join, exists
from shutil import which
import argparse, textwrap

from Phasing.BEAGLE_Phase import Phasing_wrapper, VCF2BEAGLE
# from Phasing.src.BEAGLE import BEAGLE
from src.PLINK2BEAGLE import PLINK2BEAGLE
from src.util import Exists, findExec

std_MAIN = "\n[Phasing]: "
std_ERROR = "\n[Phasing::ERROR]: "
std_WARNING = "\n[Phasing::WARNING]: "

class HATK_Phasing(object):

    """
    An action class for HLA fine-mapping analysis.
    """


    ### External Software ###
    plink = findExec("plink", std_ERROR+"'plink' command can't be found. Please install it.")

    beagle = findExec("beagle", std_ERROR+"'beagle' command can't be found. Please install it.")
    java = findExec("java", std_ERROR+"'java' command can't be found. Please install it.")

    beagle2vcf = findExec("beagle2vcf.jar",
        std_ERROR+"'beagle2vcf.jar' module can't be found.\n"
                  "Please (1) Download that jar file from 'https://faculty.washington.edu/browning/beagle_utilities/beagle2vcf.jar' AND \n"
                  "(2) Locate it in the 'dependency/' folder.")
    vcf2beagle = findExec("vcf2beagle.jar",
        std_ERROR+"'vcf2beagle.jar' module can't be found.\n"
                  "Please (1) Download that jar file from 'https://faculty.washington.edu/browning/beagle_utilities/vcf2beagle.jar' AND \n"
                  "(2) Locate it in the 'dependency/' folder.")
    linkage2beagle = findExec("linkage2beagle.jar",
        std_ERROR+"'linkage2beagle.jar' module can't be found.\n"
                  "Please (1) Download that jar file from 'https://faculty.washington.edu/browning/beagle_utilities/linkage2beagle.jar' AND \n"
                  "(2) Locate it in the 'dependency/' folder.")


    def __init__(self, _out, _bfile=None, _vcf=None, _bgl=None, _markers=None,
                 _java_mem='1G', _nthreads=1, _f_save_intermediates=False, _f_gzip=True, _f_CookHLA_REF=False):

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
                                                  self.linkage2beagle, self.plink, self.java, _java_mem,
                                                  _f_save_intermediates, _f_CookHLA_REF)
        elif self.input_file_format== "VCF":
            # use 'VCF2BEAGLE'
            pass
        else:
            self.bgl = _bgl
            self.markers = _markers

        ### Phasing (Main)
        self.bgl_phased = \
            Phasing_wrapper(self.bgl, self.markers, self.out_prefix, self.beagle2vcf, self.vcf2beagle, self.beagle,
                            self.java, _java_mem, _nthreads, _f_save_intermediates, _f_gzip)

        self.markers = shutil.copy(self.markers, self.markers.rstrip(".P2B.markers") + ".markers")

        # print(self.markers)
        # print(self.bgl_phased)

    def remove_P2B_input(self):
        if Exists(self.out_prefix + ".P2B.bgl"): os.remove(self.out_prefix + ".P2B.bgl")
        if Exists(self.out_prefix + ".P2B.markers"): os.remove(self.out_prefix + ".P2B.markers")


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
