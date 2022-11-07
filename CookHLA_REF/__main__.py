#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join
import argparse, textwrap

from src.HATK_Error import HATK_InputPreparation_Error, RaiseError
from src.PLINK_Bash import Bash_RUN_PLINK
from bMarkerGenerator.bMarker import bMarker
from Phasing.__main__ import HATK_Phasing
from CookHLA_REF.isBIMwithOldLabel import isBIMwithOldLabel, updateNameANDAlleles
from CookHLA_REF.exclude_AA_SNPS_INS import exclude_AA_SNPS_INS
from src.util import findExec, Exists

std_MAIN = "\n[CookHLA_REF]: "
std_ERROR = "\n[CookHLA_REF::ERROR]: "
std_WARNING = "\n[CookHLA_REF::WARNING]: "


class HATK_CookHLA_REF(object):

    plink = findExec("plink", std_ERROR+"'plink' command can't be found. Please install PLINK.")

    def __init__(self, _out_prefix, _bmarker=None, _bgl_phased=None,
                 _markers=None, _java_mem='1G', _nthreads=1, _f_only_SNP_HLA=True, _f_save_intermediates=False):

        self.out_dir = dirname(_out_prefix)
        self.out_prefix = _out_prefix
        self.out_prefix2 = join(self.out_dir, basename(self.out_prefix))

        self.bMARKER = bMarker(_bmarker) if _bmarker else None

        self.bgl = None
        self.markers = None

        self.bgl_phased = _bgl_phased if _bgl_phased else None # Debug

        self.HATK_Phasing = None

        # conditions
        self.f_only_SNP_HLA = _f_only_SNP_HLA
        self.to_exclude = None

        self.f_isBIMwithOldLabel = False
        self.update_name = None
        self.update_alleles = None


        """
        bMarker: *.{bed,bim,fam,FRQ.frq}
        bgl_phased: *.{bgl.phased,markers}
        
        2 components, 6 files are required.
        """

        # Preparing PLINK command for Converting/Subsetting bMarker
        """
        - Does all bim markers have old/new(HATK) label?
        - Exclude AA + SNPS + INS ?
        """

        command = [
            self.plink, "--make-bed",
            "--bfile {}".format(self.bMARKER.file_prefix),
            "--out {}".format(self.out_prefix2),
            "--allow-no-sex",
            "--keep-allele-order"
        ]

        self.f_isBIMwithOldLabel = isBIMwithOldLabel(self.bMARKER.bim)

        if not self.f_isBIMwithOldLabel:
            print(std_MAIN.lstrip("\n") +
                  "Marker label and Allele code will be converted to the old version. "
                  "(e.g. 'HLA_A*01:01 -> 'HLA_A_0101' / (A1='p', A2='a') -> ('P', 'A'))")
            self.update_name , self.update_alleles = \
                updateNameANDAlleles(self.bMARKER.bim, self.out_prefix2 + '.update_name',
                                     self.out_prefix2 + '.update_allele')
            command.append("--update-name {}".format(self.update_name))
            command.append("--update-alleles {}".format(self.update_alleles))

        if _f_only_SNP_HLA:
            print(std_MAIN.lstrip("\n") +
                  "Only SNP and HLA markers will be left. "
                  "(e.g. 'rs7762289' and 'HLA_A_0101')")
            self.to_exclude = exclude_AA_SNPS_INS(self.bMARKER.bim, self.out_prefix2+".to_exclude")
            command.append("--exclude {}".format(self.to_exclude))


        ## PLINK bash execution. (MAIN 1)
        bmarker_filtered = \
            Bash_RUN_PLINK(' '.join(command), self.out_prefix2, _f_save_log=True, _f_save_intermediates=_f_save_intermediates)

        if not _f_save_intermediates:
            if Exists(self.to_exclude): os.remove(self.to_exclude)
            if Exists(self.update_name): os.remove(self.update_name)
            if Exists(self.update_alleles): os.remove(self.update_alleles)

        self.bMARKER = bMarker(bmarker_filtered) # (***) The bMarker to be used for reference panel.

        ## FRQ generation (MAIN 1-1)
        if not Exists(self.bMARKER.FRQ):
            self.bMARKER.genFRQ()

        ## Phasing with the given bMarker. (MAIN 2)
        self.HATK_Phasing = \
            HATK_Phasing(_out_prefix, _bfile=self.bMARKER.file_prefix,
                         _java_mem=_java_mem, _nthreads=_nthreads, _f_save_intermediates=_f_save_intermediates,
                         _f_CookHLA_REF=True)

        self.bgl_phased = self.HATK_Phasing.bgl_phased
        self.markers = self.HATK_Phasing.markers

        if not _f_save_intermediates:
            self.HATK_Phasing.remove_P2B_input()

        """
        Intentionally, Beagle files('*.bgl' and '*.markers' files will be generated from bMarker PLINK file.
        (To deal with A1, A2 allele position problem in '*.markers' file.
        """

        # print(self.__repr__())


    def __repr__(self):

        str_bMARKER = \
            "- bMarker:\n{}\n".format(self.bMARKER)

        str_bgl_phased = \
            "- Phased beagle file: {}\n".format(self.bgl_phased)

        str_markers = \
            "- Beagle markers file: {}\n".format(self.markers)

        str_Phasing = \
            "- Phasing information:\n{}\n".format(self.HATK_Phasing)

        str_summary = \
            ''.join([
                str_bMARKER, str_bgl_phased, str_markers
            ]).rstrip("\n")
        return str_summary







if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='CookHLA_REF',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        [CookHLA_REF]
        
        - To generate a custom reference panel for CookHLA.
        - (1) HLA type and (2) Genotype data are required
        - HLA type should be from 'Sequencing-Based Typing' or '30x-WGS based inferred ones'.
            (i.e. not imputed ones.)
        - Reference panel generated with imputed HLA type could cause poor accuracy of HLA imputation. 

    #################################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--bmarker", help="\nBinary marker genotype file generated by the 'bMarkerGenerator'. (.bed/.bim/.fam)\n\n",
                        required=True)

    parser.add_argument("--bgl-phased", help="Phased beagle file(ex. \"*.bgl.phased\"). (DEBUG)\n\n")
    parser.add_argument("--markers", help="Beagle markers file(ex. \"*.markers\"). (DEBUG)\n\n")

    parser.add_argument("--out", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')

    # java
    parser.add_argument("--java-mem", "-mem",
                        help="\nHeap memory size for each BEAGLE implementation (default: 2g).\n\n", default="2g")

    parser.add_argument("--nthreads", "-nth",
                        help="\nThe number of threads for each BEAGLE implementation (default: 1).\n\n", default=1,
                        type=int)


    ##### <for Test> #####

    # args = parser.parse_args(["--variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
    #                           "--out", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    ref = HATK_CookHLA_REF(args.out, _bmarker=args.bmarker, _bgl_phased=args.bgl_phased, _markers=args.markers,
                           _java_mem=args.java_mem, _nthreads=args.nthreads, _f_save_intermediates=args.save_intermediates)
    print(ref)