# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
from shutil import which
from mpmath import log10


p_PLINK = which("plink")
p_RSCRIPT = which("Rscript")

def MakeLDMatrix(_bfile, _out, _p_plink):

    command = [_p_plink, "--noweb", "--r2", "square", "--bfile", _bfile, "--out", _out]
    command = ' '.join(command)

    print(command)
    os.system(command)


    return (_out + ".ld")


def manhattan(_bfile,
              _lr, _out,
              _pointcol = "\#FF0000",
              _min_pos = "29.60E6", _max_pos = "33.2E6",
              _gb = "data/HLA_Analysis/Plotting/manhattanplot/known_genes_build36/known_genes_build36_soumya.txt",
              _p_Rscript = which("Rscript"), _p_plink = which("plink"),
              _Rsrc = "src/HLA_Analysis/Plotting/manhattanplot/manhattan_HLA_HATK.R"):

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    std_MAIN_PROCESS_ERROR = "\n[%s::ERROR]: " % (os.path.basename(__file__))

    print(std_MAIN_PROCESS_NAME + "Conducting Manhattna Plotting.")


    _out = _out if not _out.endswith('/') else _out.rstrip('/')

    INTERMEDIATE_PATH = os.path.dirname(_out)

    if not os.path.exists(INTERMEDIATE_PATH):
        os.system(' '.join(["mkdir", "-p", INTERMEDIATE_PATH]))


    ##### < Argument Checking. > #####

    if _pointcol.startswith('#'):
        _pointcol = "\\" + _pointcol


    ##### < Making LD matrix. > #####

    if not os.path.isfile(_out+".ld"):

        print(std_MAIN_PROCESS_NAME + "Making LD Matrix.")
        f_LDMatrix =  MakeLDMatrix(_bfile, _out, _p_plink)

    else:

        print(std_MAIN_PROCESS_NAME + "LD Matrix {0} is found. Skipping making it.".format(_out+".ld"))
        f_LDMatrix = _out+".ld"




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



    ##### < Extracting R^2 data. > #####

    LD_Matrix = pd.read_table(f_LDMatrix, sep='\t| ', engine='python', header=None, names=MARKER_set).loc[:, TOP_LABEL]

    LD_Matrix.index = pd.Index(MARKER_set)

    print(std_MAIN_PROCESS_NAME + "Loading LD Matrix data.\n")
    print(LD_Matrix.head())
    print(LD_Matrix.tail())

    LD_Matrix.dropna().to_csv(_out+".r2totophit", sep='\t', header=False, index=True, na_rep="NA")


    # So far, Necessary input files have been prepared.




    ##### < Executing Rscript. > #####

    command = [_p_Rscript, _Rsrc,
               _lr, _out,
               _pointcol, _min_pos, _max_pos, TOP_LABEL, str(_yaxis),
               _gb, _out+".r2totophit"]

    command = ' '.join(command)
    print(command)

    os.system(command)



    ##### < Removing .ld and .r2totophit files. > #####




    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################
           
        manhattan.py
        
        Explanation.
        
        
        
    #################################################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--bfile", help="\nInput Data file.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("--logistic-result", "-lr", help="\nOutput from logistic regression(\".assoc.logstic\").\n\n", required=True)

    parser.add_argument("--gene-build", "-gb", help="\nGene Build file.\n\n", required=True)




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
    #                           ])


    ### for Publish
    args = parser.parse_args()
    # print(args)

    manhattan(args.bfile,
              args.logistic_result, args.out,
              _Rsrc="/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/src/HLA_Analysis/Plotting/manhattanplot/manhattan_HLA_HATK.R",
              _gb=args.gene_build)