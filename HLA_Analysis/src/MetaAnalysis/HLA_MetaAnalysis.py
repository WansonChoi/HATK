# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
import numpy as np

from math import log
from scipy.stats import norm

########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


ComplementBase = {
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C'
}


def HATK_MetaAnalysis():


    return 0



def HLA_MetaAnalysis(_study1_lr, _study2_lr, _out=None, _study1_m=None, _study2_m=None):



    ########## < Checking files > ##########


    ### Logistic regression result.

    # study1
    df_study1_lr = isLogisticResult(_study1_lr)

    if df_study1_lr.shape[0] == 0:
        sys.exit()

    # print("df_study1_lr:\n{}\n".format(df_study1_lr))

    # study2
    df_study2_lr = isLogisticResult(_study2_lr)

    if df_study2_lr.shape[0] == 0:
        sys.exit()

    # print("df_study2_lr:\n{}\n".format(df_study2_lr))


    df_Main = df_study1_lr.merge(df_study2_lr, left_on=['SNP'], right_on=['SNP'])[['SNP', 'BETA_x', 'SE_x', 'BETA_y', 'SE_y']]
    # print("df_Main:\n{}\n".format(df_Main))




    ### (Optional) Marker information file (*.bim or *.markers)

    Flag_Markers = False

    if bool(_study1_m) and bool(_study2_m):

        # study1
        df_study1_m = isMarkerFile(_study1_m)

        if df_study1_m.shape[0] == 0:
            sys.exit()

        # print("df_study1_m:\n{}\n".format(df_study1_m))

        # study2
        df_study2_m = isMarkerFile(_study2_m)

        if df_study2_m.shape[0] == 0:
            sys.exit()

        # print("df_study2_m:\n{}\n".format(df_study2_m))


        if (df_study1_m.shape[0] > 0 and df_study2_m.shape[0] > 0):

            ### Applying Flipping
            f = Flip(df_Main, df_study1_m, df_study2_m)

            if type(f) == int and f == -1:
                print(std_ERROR_MAIN_PROCESS_NAME + "Please check the inner contents of marker information files.\n"
                                                    "('{file1}', '{file2}')".format(file1=_study1_m, file2=_study2_m))
                sys.exit()

            df_Main = f ## Flip information is applied here.

            # print("df_Main(Flipped) : \n{}\n".format(df_Main))

            Flag_Markers = True




    """
    # Inverse-Variance Method
    
    wk=1/sek**2
    wj=1/sej**2
    ivw=(betak*wk+betaj*wj)/(wk+wj)
    ivw.se=sqrt(1/(wk+wj))
    meta.z=ivw/ivw.se
    meta.p=2*pnorm(-abs(meta.z)) 
    """

    # df_Main.to_csv(_out+'.ongoing.txt', sep='\t', header=True, index=False)


    np.seterr(divide='ignore')


    w_study1 = 1/(np.square(df_Main['SE_x'].to_numpy()))
    w_study2 = 1/(np.square(df_Main['SE_y'].to_numpy()))

    beta_study1 = df_Main['BETA_x'].to_numpy()
    beta_study2 = df_Main['BETA_y'].to_numpy()

    ivw = (beta_study1*w_study1 + beta_study2*w_study2)/(w_study1 + w_study2)
    ivw_se = np.sqrt(1/(w_study1+w_study2))
    # print("ivw_se:{}".format(ivw_se))

    meta_z = np.divide(ivw,ivw_se)
    # print("meta_z : {}".format(meta_z))
    meta_p = 2*norm.cdf(-abs(meta_z))

    df_Main['BETA'] = ivw
    df_Main['SE'] = ivw_se
    df_Main['MetaP'] = meta_p
    df_Main = df_Main[['SNP', 'BETA', 'SE', 'MetaP']]

    # print(df_Main)

    if bool(_out):
        df_Main.to_csv(_out+'.meta', sep='\t', header=True, index=False)
        return _out+'.meta'
    else:
        return df_Main




def isLogisticResult(_lr):


    df_NULL = pd.DataFrame([])


    if os.path.exists(_lr):

        if _lr.endswith('.assoc.logistic'):

            # PLINK Logistic Regression Result.
            df_RETURN = pd.read_csv(_lr, sep='\s+', header=0, usecols=['SNP', 'OR', 'SE', 'STAT', 'P'])

            df_RETURN['BETA'] = df_RETURN['OR'].map(lambda x : log(x))

            return df_RETURN

        else:
            # Check Header manually.

            f_lr = open(_lr, 'r')

            l_Header = f_lr.readline()
            l_Header = re.split(r'\s+', l_Header.rstrip('\n'))

            # print(l_Header)


            if 'SNP' not in l_Header:
                print(std_ERROR_MAIN_PROCESS_NAME + "Column 'SNP' can't be found in the header('{}').\n"
                                                    "Please check it again.".format(l_Header))
                return df_NULL

            l_having = []

            # Columns 'BETA' and 'SE' are needed.

            if not ('BETA' in l_Header or 'OR' in l_Header):
                # At least either 'BETA' or 'OR' column is needed.
                print(std_ERROR_MAIN_PROCESS_NAME + "Neither 'BETA' or 'OR'(Odds Ratio) column can be found in the header('{}').\n"
                                                    "Please check it again.".format(l_Header))
                return df_NULL

            else:

                if 'BETA' in l_Header:
                    l_having.append('BETA')
                elif 'OR' in l_Header:
                    l_having.append('OR')


            if not ('SE' in l_Header or 'P' in l_Header):
                print(std_ERROR_MAIN_PROCESS_NAME + "Neither 'SE' or 'P' column can be found in the header('{}').\n"
                                                    "Please check it again.".format(l_Header))
                return df_NULL

            else:
                if 'SE' in l_Header:
                    l_having.append('SE')
                elif 'P' in l_Header:
                    l_having.append('P')

            f_lr.close()


            df_RETURN = pd.read_csv(_lr, sep='\s+', header=0, usecols=['SNP']+l_having)

            if not 'BETA' in l_Header and 'OR' in l_Header:
                # Make 'BETA' based on 'OR'
                df_RETURN['BETA'] = df_RETURN['OR'].map(lambda x : log(x))

            if not 'SE' in l_Header and 'P' in l_Header:
                # Make 'SE' based on 'P'
                df_RETURN['SE'] = (df_RETURN['BETA'] / df_RETURN['P'].map(lambda x : norm.ppf(x/2))).map(lambda x : abs(x))

            return df_RETURN


    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Given Logistic Regression result file('{}') doesn't exist.\n"
                                            "Please check it again.".format(_lr))
        return df_NULL




def isMarkerFile(_m):

    df_NULL = pd.DataFrame([])


    if os.path.exists(_m):

        if _m.endswith('.bim'):
            return pd.read_csv(_m, sep='\s+', header=None, names=['Chr', 'SNP', 'GD', 'BP', 'al1', 'al2'], usecols=['SNP', 'BP', 'al1', 'al2'])
        elif _m.endswith('.markers'):
            return pd.read_csv(_m, sep='\s+', header=None, names=['SNP', 'BP', 'al1', 'al2'])
        else:
            # Check Header manually.
            f_m = open(_m, 'r')

            l_Header = f_m.readline()
            l_Header = re.split(r'\s+', l_Header.rstrip('\n'))

            print(l_Header)

            """
            Temporarily put off.
            """

            return df_NULL




    else:

        print(std_ERROR_MAIN_PROCESS_NAME + "Given marker information file('{}') doesn't exist.\n"
                                            "Please check it again.".format(_m))
        return df_NULL




def Flip(_df_Main, _df_study1_m, _df_study2_m):

    # print(std_MAIN_PROCESS_NAME + "Flipping.")


    df_markers = _df_study1_m.merge(_df_study2_m, left_on='SNP', right_on='SNP')
    # print("df_markers:\n{}\n".format(df_markers))


    df_Main2 = _df_Main.merge(df_markers, left_on='SNP', right_on='SNP', how='left')
    # print("df_Main2:\n{}\n".format(df_Main2))


    f_toComplement = df_Main2[['al1_x', 'al2_x', 'al1_y', 'al2_y']] \
                        .apply(set, axis=1) \
                        .map(lambda x : len(x) == 4)

    # print("f_toComplement:\n{}\n".format(f_toComplement))

    df_Main2['toComplement'] = f_toComplement


    count = 0

    l_toFlip = []

    for item in df_Main2[['SNP', 'al1_x', 'al2_x', 'al1_y', 'al2_y', 'toComplement']].itertuples():

        [idx, SNP, al1_x, al2_x, al1_y, al2_y, toComplement] = item

        if toComplement:
            al1_y = ComplementBase[al1_y]
            al2_y = ComplementBase[al2_y]


        ### Main Classification for Flipping.
        if al1_x == al1_y and al2_x == al2_y:
            # Nothing to Flip
            l_toFlip.append(False)
        elif al1_x == al2_y and al2_x == al1_y:
            l_toFlip.append(True)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Wrong allele characters.\n"
                                                "{SNP} : ({al1_x}, {al2_x}, {al1_y}, {al2_y})\n"
                                                "".format(al1_x=al1_x, al2_x=al2_x, al1_y=al1_y, al2_y=al2_y, SNP=SNP))
            return -1



        count += 1
        # if count > 5: break

    # if len(l_toFlip) != df_Main2.shape[0]:
    #     print(std_ERROR_MAIN_PROCESS_NAME + "Wrong output of `l_toFlip`.")
    #     return -1

    sr_toFlip = pd.Series(l_toFlip)
    df_Main2['FLIP'] = sr_toFlip

    df_Main2['BETA_y'] = df_Main2['BETA_y'] * sr_toFlip.map(lambda x : -1 if x else 1)
    # df_Main2.drop('toComplement', inplace=True, axis=1)


    return df_Main2




if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        HLA_MetaAnalysis.py
        
        MetaAnalysis with Inverse-variance weighting method.


    #################################################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--s1-logistic-result", "-s1lr", help="\nLogistic Regression Result file of Study 1(ex. *.assoc.logistic).\n\n", required=True)
    parser.add_argument("--s1-markers", "-s1m", help="\nMarker information which used in Study 1 Logistic Regression(ex. *.bim, *.markers).\n\n")

    parser.add_argument("--s2-logistic-result", "-s2lr", help="\nLogistic Regression Result file of Study 2(ex. *.assoc.logistic).\n\n", required=True)
    parser.add_argument("--s2-markers", "-s2m", help="\nMarker information which used in Study 2 Logistic Regression(ex. *.bim, *.markers).\n\n")

    parser.add_argument("--out", "-o", help="\nOuput file prefix\n\n", required=True)



    ### for Testing

    # Linux
    # args = parser.parse_args(['-s1lr', '03-meta-of-binary/All_CD.assoc.logistic',
    #                           '-s1m', '03-meta-of-binary/All_CD.bim',
    #                           '-s2lr', '03-meta-of-binary/RIKEN_CD_HLA-imputation-binary.logistic.txt',
    #                           '-s2m', '03-meta-of-binary/RIKEN_CD_HLA-imputation-binary.logistic.markers',
    #                           '-o', '03-meta-of-binary/meta.out.python.txt'])

    # OS X
    # args = parser.parse_args(['-s1lr', '/Users/wansun/Git_Projects/HATK/tests/03-meta-of-binary/All_CD.assoc.logistic',
    #                           '-s1m', '/Users/wansun/Git_Projects/HATK/tests/03-meta-of-binary/All_CD.bim',
    #                           '-s2lr', '/Users/wansun/Git_Projects/HATK/tests/03-meta-of-binary/RIKEN_CD_HLA-imputation-binary.logistic.txt',
    #                           '-s2m', '/Users/wansun/Git_Projects/HATK/tests/03-meta-of-binary/RIKEN_CD_HLA-imputation-binary.logistic.markers',
    #                           '-o', '/Users/wansun/Git_Projects/HATK/tests/03-meta-of-binary/HATK_meta'])


    ### for Publish
    args = parser.parse_args()
    print(args)

    HLA_MetaAnalysis(args.s1_logistic_result, args.s2_logistic_result, args.out, args.s1_markers, args.s2_markers)