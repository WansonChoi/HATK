# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
from shutil import which
from mpmath import log10


p_PLINK = which("plink")
p_RSCRIPT = which("Rscript")



def manhattan(_bfile,
              _lr, _out,
              _pointcol = "\#778899", _topcol = "\#FF0000",
              _min_pos = "29.60E6", _max_pos = "33.2E6",
              _gb = "data/HLA_Analysis/Plotting/manhattanplot/known_genes_build36/known_genes_build36_soumya.txt",
              _p_Rscript = which("Rscript"), _p_plink = which("plink"),
              _Rsrc = "src/HLA_Analysis/Plotting/manhattanplot/manhattan_HLA_HATK.R"):

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_MAIN_PROCESS_ERROR = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting Manhattan Plotting.")


    _out = _out if not _out.endswith('/') else _out.rstrip('/')

    INTERMEDIATE_PATH = os.path.dirname(_out)

    if not os.path.exists(INTERMEDIATE_PATH):
        os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))


    ##### < Argument Checking. > #####

    if _pointcol.startswith('#'):
        _pointcol = "\\" + _pointcol




    ##### < Loading .assoc.logistic file. > #####

    print(std_MAIN_PROCESS_NAME + "Loading \".assoc.logistic\".\n")
    logistic_result = pd.read_table(_lr, sep='\s+', engine='python', header=0, usecols=["SNP", "P"])

    print(logistic_result.head())

    MARKER_set = logistic_result.iloc[:, 0].tolist()
    # print(std_MAIN_PROCESS_NAME + "Marker Labels are \n{0}".format(MARKER_set))


    logistic_result = logistic_result.sort_values("P")

    TOP_LABEL = logistic_result.iat[0, 0]
    TOP_LABEL_value = logistic_result.iat[0, 1]
    print(std_MAIN_PROCESS_NAME + "Top signal is \"{0}(P-value : {1})\"".format(TOP_LABEL, TOP_LABEL_value))


    temp = -log10(TOP_LABEL_value)

    if temp <= 10:
        _yaxis = 10
    else:

        if temp % 5 == 0:
            _yaxis = temp
        else:
            _yaxis = (int(temp/5) + 1)*5


    print(std_MAIN_PROCESS_NAME + "maximum y-axis is {0}".format(_yaxis))



    ##### < Executing Rscript. > #####

    command = [_p_Rscript, _Rsrc,
               _lr, _out,
               _pointcol, _topcol,
               _min_pos, _max_pos, TOP_LABEL, str(_yaxis),
               _gb]

    command = ' '.join(command)
    print(command)

    os.system(command)






    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################
           
        manhattan.py
        
        Conducting manhattan plotting.
        
        
        
    #################################################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--bfile", help="\nInput Data file.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("--logistic-result", "-lr", help="\nOutput from logistic regression(\".assoc.logstic\").\n\n", required=True)
    parser.add_argument("--gene-build", "-gb", help="\nGene Build file.\n\n", required=True)

    parser.add_argument("--point-color", "-pc", help="\nPoint color(ex. \"#778899\").\n\n")
    parser.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n")

    parser.add_argument("--Rsrc", help="\nTop signal point color(ex. \"#FF0000\").\n\n",
                        default="src/HLA_Analysis/Plotting/manhattanplot/manhattan_HLA_HATK.R")


    # (2018. 8. 10.)
    # Arguments for "hg", "min.pos" and "max.pos" should be added later.


    ### for Testing

    # # AA
    # args = parser.parse_args(["--bfile", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/manhattanplot/data_Rev_merged.AA.CODED",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/MANHATTAN_test/data_Rev_merged.AA.CODED",
    #                           "-lr", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/manhattanplot/data_Rev_merged.AA.CODED.assoc.logistic",
    #                           "-gb", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/manhattanplot/known_genes_build36/known_genes_build36_soumya_chr6.txt",
    #                           ])

    # ## Original Sample
    # args = parser.parse_args(["--bfile", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/sample/Merged/merged",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/MANHATTAN_test2/merged",
    #                           "-lr", "/Users/wansun/Git_Projects/UC-CD-HLA/UC-CD-HLA/manhattan_plot/Figure1/Assoc_result_Han/01-association/All_CD.assoc.logistic",
    #                           "-gb", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/data/HLA_Analysis/Plotting/manhattanplot/known_genes_build36/known_genes_build36_soumya_chr6.txt",
    #                           "-pc", "#E0FFFF",
    #                           "--Rsrc", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/src/HLA_Analysis/Plotting/manhattanplot/manhattan_HLA_HATK.R"
    #                           ])


    ### for Publish
    args = parser.parse_args()
    # print(args)

    manhattan(args.bfile,
              args.logistic_result, args.out,
              _gb=args.gene_build,
              _pointcol=args.point_color,
              _Rsrc=args.Rsrc,
              )