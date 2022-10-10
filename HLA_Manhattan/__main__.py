#-*- coding: utf-8 -*-
import os, sys, re
from os.path import basename, dirname, join, isdir
from shutil import which
import numpy as np
import argparse, textwrap

from HLA_Manhattan.manhattan import Manhattan, Manhattan_OM
from src.util import Exists, checkFile, findExec
from src.HATK_Error import HATK_InputPreparation_Error, RaiseError

std_MAIN = "\n[Manhattan]: "
std_ERROR = "\n[Manhattan::ERROR]: "
std_WARNING = "\n[Manhattan::WARNING]: "


class HATK_Manhattan(object):

    # External software
    Rscript = findExec("Rscript", std_ERROR+"'Rscript' command can't be found. Please install R.")

    def __init__(self, _assoc:list, _out_prefix, _hg,
                 _point_color="#778899", _top_color="#FF0000", _point_size="15", _yaxis_unit="10",
                 _f_save_intermediates=False):

        """
        - given assocs are all same type, i.e. all PLINK assoc or all Omnibus assoc
        - All given assocs exist
        - How many assocs?
        """

        ### Main Variables ###
        self.assoc = [checkFile(item) for item in _assoc]
        self.N_assoc = len(self.assoc)
        self.out_prefix = _out_prefix
        self.out_dir = dirname(self.out_prefix)
        self.hg = _hg

        self.point_color = _point_color
        self.top_color = _top_color
        self.point_size = _point_size
        self.yaxis_unit = _yaxis_unit

        self.assoc_type = None
        self.ManhattanPlot = [] # Main result

        ### Main Actions ###
        # Type check of assocs. (is PLINK assoc or Omnibus?)
        func1 = np.vectorize(lambda x: x.endswith('.assoc.logistic') or x.endswith('.assoc.linear'))
        isPLINK = np.all(func1(self.assoc))
        func2 = np.vectorize(lambda x: x.endswith('.omnibus'))
        isOMNIBUS = np.all(func2(self.assoc))

        if isPLINK:
            self.assoc_type = "PLINK"
        elif isOMNIBUS:
            self.assoc_type = "OMNIBUS"
        else:
            raise HATK_InputPreparation_Error(
                std_ERROR + "Given assoociation test results are not homogenous.\n" +
                "Please check whether given assoc results are all from either PLINK('*.assoc.{logistic,linear}') or Omnibus Test('*.omnibus')."
            )


        # print(self.__repr__())
        if self.assoc_type == "PLINK":
            self.ManhattanPlot = \
                Manhattan(self.assoc, self.out_prefix, self.hg,
                          _pointcol=self.point_color, _topcol=self.top_color,
                          _pointsize=self.point_size, _yaxis_unit=self.yaxis_unit,
                          _Rscript=self.Rscript, _f_save_intermediates=_f_save_intermediates)
        elif self.assoc_type == "OMNIBUS":
            self.ManhattanPlot = \
                Manhattan_OM(self.assoc, self.out_prefix, self.hg, _pointcol=self.point_color, _topcol=self.top_color,
                             _pointsize=self.point_size, _yaxis_unit=self.yaxis_unit, _Rscript=self.Rscript,
                             _f_save_intermediates=_f_save_intermediates)

        # print(self.__repr__())


    def __repr__(self):
        str_assoc = \
            "- Association test result files: {}\n".format(self.assoc)
        str_assoc_type = \
            "- Assoc test type: {}\n".format(self.assoc_type)
        str_N_assoc = \
            "- The number of given assoc results: {}\n".format(self.N_assoc)
        str_outprefix = \
            "- Output prefix: {}\n".format(self.out_prefix)
        str_hg = \
            "- Human Genome version: {}\n".format(self.hg)

        str_point_color = \
            "- Point color: {}\n".format(self.point_color)
        str_top_color = \
            "- Top color: {}\n".format(self.top_color)
        str_point_size = \
            "- Point size: {}\n".format(self.point_size)
        str_yaxis_unit = \
            "- Y-axis unit: {}\n".format(self.yaxis_unit)

        str_manhattan = \
            "- Manhattan plot(s): {}\n".format(
                self.ManhattanPlot if self.assoc_type == "PLINK" else \
                ''.join("\n   - {}".format(item) for item in self.ManhattanPlot)
            )


        str_summary = ''.join([str_assoc, str_assoc_type, str_N_assoc, str_outprefix, str_hg,
                               str_point_color, str_top_color, str_point_size, str_yaxis_unit,
                               str_manhattan])
        return str_summary



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Manhattan',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        Manhattan Plot.


    #################################################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(ex. *.assoc.logistic).\n\n",
                        nargs='+', required=True)
    # parser.add_argument("--plot-label", "-pl", help="\nPlot Label\n\n", default="")
    parser.add_argument("--out", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("--hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)

    parser.add_argument("--point-color", "-pc", help="\nPoint color(ex. \"#778899\").\n\n", default="#778899")
    parser.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n",
                        default="#FF0000")

    parser.add_argument("--point-size", "-ps", help="\nGeneral point size (default: 15).\n\n", default="15")
    parser.add_argument("--yaxis-unit", "-yau", help="\nY-axis value(-log10(x)) unit (default : 10).\n\n",
                        default="10")

    # parser.add_argument("--HLA", help="\nHLA gene to plot heatmap.\n\n",
    #                     choices=['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'], nargs='+')

    ### for Testing

    # args = parser.parse_args(['-ar', '/Users/wansun/Git_Projects/HATK/MyOmnibusTest2/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.RA.AA_DRB1_96.omnibus',
    #                           '-hg', '18',
    #                           '-o', '/Users/wansun/Git_Projects/HATK/MyOmnibusTest2/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.RA.AA_DRB1_96.omnibus',
    #                           "--HLA", 'A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1'])

    ### for Publish
    args = parser.parse_args()
    print(args)

    # Manhattan(args.assoc_result, args.out, args.hg, _pointcol=args.point_color, _topcol=args.top_color,
    #           _pointsize=args.point_size, _yaxis_unit=args.yaxis_unit, _HLA=args.HLA)

    HATK_Manhattan(args.assoc_result, args.out, args.hg, _point_color=args.point_color, _top_color=args.top_color,
                   _point_size=args.point_size, _yaxis_unit=args.yaxis_unit)