# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
from shutil import which
from mpmath import log10


p_RSCRIPT = which("Rscript")


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def HATK_manhattan(_results_assoc, _plot_label, _out, _hg, _pointcol="#778899", _topcol="#FF0000", _min_pos="29.60E6", _max_pos="33.2E6",
                   _p_Rscript=which("Rscript")):


    if not isinstance(_results_assoc, list):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--logistic-result(-lr)\" wansn't given as list. Please check it again.\n")
    else:

        for item in _results_assoc:
            if not os.path.exists(item):
                print(std_ERROR_MAIN_PROCESS_NAME + "The given association result file({0}) doesn't exist. Please check it again.\n".format(item))
                sys.exit()


    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--out\" wasn't given. Please check it again.\n")
        sys.exit()

    if not bool(_hg):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hg"))
        sys.exit()

    if not bool(_plot_label):
        _plot_label = " "


    return manhattan(_results_assoc, _plot_label, _out, _hg, _pointcol="#778899", _topcol="#FF0000", _min_pos="29.60E6",
                     _max_pos="33.2E6", _p_Rscript=which("Rscript"))




def manhattan(_assoc_result, _out, _hg, _plot_label = "",
              _pointcol="#778899", _topcol="#FF0000", _min_pos="29.60E6", _max_pos="33.2E6",
              _pontsize="15", _yaxis_unit = "10", _p_Rscript=which("Rscript")):



    ########## < Core Variables > ##########

    # Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)

    # Point color
    _pointcol = re.escape(_pointcol if _pointcol.startswith('#') else ("#"+_pointcol))
    _topcol = re.escape(_topcol if _topcol.startswith('#') else ("#"+_topcol))


    # Paths
    p_src = "./src"
    p_data = "./data"



    ########## < Dependency and Argument Checking > ##########

    if not isinstance(_assoc_result, list):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--logistic-result(-lr)\" wansn't given as list. Please check it again.\n")
        sys.exit()

    # Logistic Regression
    print(std_MAIN_PROCESS_NAME + "Loading \".assoc.logistic\".\n")
    print("{0} logistic regression results are given.\n".format(len(_assoc_result)))


    l_TOP_LABEL = []
    l_yaxis = []

    for i in range(0, len(_assoc_result)):

        print("\n[{0}] : {1}\n".format(i, _assoc_result[i]))

        t_lr = pd.read_table(_assoc_result[i], sep='\s+', header=0, usecols=["SNP", "P"]).dropna().sort_values("P")

        # MARKER_set = t_lr.iloc[:, 0].tolist()
        # print(std_MAIN_PROCESS_NAME + "Marker Labels are \n{0}".format(MARKER_set))

        print(t_lr.head())

        TOP_LABEL = t_lr.iat[0, 0]
        TOP_LABEL_value = t_lr.iat[0, 1]

        print("\nTop signal is \"{0}(P-value : {1})\"".format(TOP_LABEL, TOP_LABEL_value))


        temp = -log10(TOP_LABEL_value)

        if temp <= 10:
            _yaxis = 10
        else:

            if temp % 5 == 0:
                _yaxis = temp
            else:
                _yaxis = (int(temp/5) + 1)*5


        print("maximum y-axis is {0}".format(_yaxis))

        l_TOP_LABEL.append(TOP_LABEL)
        l_yaxis.append(str(_yaxis))



    export_ar = ','.join(_assoc_result)
    export_TOP_LABEL = ','.join(l_TOP_LABEL)
    export_yaxis = ','.join(l_yaxis)



    # hg (Human Genome)
    _knownGene = os.path.join(p_data, 'known_genes/known_genes_soumya_chr6.hg{0}.txt'.format(_hg))


    print("\n")


    ########## < Plotting Manhattan > ##########

    print(std_MAIN_PROCESS_NAME + "Plotting Manhattan.\n")

    command = [_p_Rscript, os.path.join(p_src, "manhattan_HLA_HATK.R"),
               export_ar, _out,
               _pointcol, _pontsize, _topcol,
               _min_pos, _max_pos, export_TOP_LABEL, export_yaxis, _yaxis_unit,
               _knownGene, p_src]

    command = ' '.join(command)
    print(command)

    os.system(command)



    return _out + ".pdf"





if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################
           
        manhattan.py
        
        Conducting manhattan plotting.
        
        
        (Example 1 : hg18 / imgt3320 / single logistic regression result.)
        
        $ python3 manhattan.py \
            -lr example/Marker_Panel.hg18.imgt3320.MERGED.assoc.logistic \
            -o Marker_Panel.hg18.imgt3320.MERGED.assoc.logistic \
            -kg known_genes/known_genes_soumya_chr6.hg18.txt \
            -pc E0FFFF
        
        
        (Example 2 : hg18 / imgt3320 / multiple logistic regression result.)

        $ python3 manhattan.py \
            -lr \
            example/Marker_Panel.hg18.imgt3320.AA.CODED.assoc.logistic \
            example/Marker_Panel.hg18.imgt3320.HLA.assoc.logistic \
            example/Marker_Panel.hg18.imgt3320.SNPS.CODED.assoc.logistic \
            -o Marker_Panel.hg18.imgt3320.EACH.assoc.logistic \
            -kg known_genes/known_genes_soumya_chr6.hg18.txt \
            --point-color 33FF33        
        
    #################################################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--assoc-result", "-ar", help="\nAssociation test result file(ex. *.assoc.logistic).\n\n", nargs='+', required=True)
    parser.add_argument("--plot-label", "-pl", help="\nPlot Label\n\n", default="")
    parser.add_argument("--out", "-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("-hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"], metavar="hg", required=True)

    parser.add_argument("--point-color", "-pc", help="\nPoint color(ex. \"#778899\").\n\n", default="#778899")
    parser.add_argument("--top-color", "-tc", help="\nTop signal point color(ex. \"#FF0000\").\n\n", default="#FF0000")

    parser.add_argument("--point-size", "-ps", help="\nGeneral point size (default: 15).\n\n", default="15")
    parser.add_argument("--yaxis-unit", "-yau", help="\nY axis value(-log10(x)) unit (default : 10).\n\n", default="10")



    ### for Testing




    ### for Publish
    args = parser.parse_args()
    print(args)

    manhattan(args.assoc_result, args.out, args.hg, _plot_label=args.plot_label, _pointcol=args.point_color, _topcol=args.top_color,
              _pontsize=args.point_size, _yaxis_unit=args.yaxis_unit)