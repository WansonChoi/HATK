# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def NfieldDictionary(_AA_dict, _SNP_dict, _maptable, _OUTPUT_FORMAT, _out=None):

    if not (0 < _OUTPUT_FORMAT < 4):
        print("Bizzare output format({})".format(_OUTPUT_FORMAT))
        sys.exit()


    if _OUTPUT_FORMAT == 1:
        p_Field = re.compile(r'^\w+\*\d{2,3}')
    elif _OUTPUT_FORMAT == 2:
        p_Field = re.compile(r'^\w+\*\d{2,3}(:\d{2,3}){0,1}[A-Z]?')
    elif _OUTPUT_FORMAT == 3:
        p_Field = re.compile(r'^\w+\*\d{2,3}(:\d{2,3}){0,2}[A-Z]?')
    else:
        p_Field = None


    df_AA = pd.read_csv(_AA_dict, '\s+', header=None, dtype=str, names=['Allele', 'Seq'])
    # print("_AA_dict :\n{}\n".format(df_AA.head()))

    df_SNP = pd.read_csv(_SNP_dict, '\s+', header=None, dtype=str, names=['Allele', 'Seq'])
    # print("_SNP_dict :\n{}\n".format(df_SNP.head()))


    d_df_maptable = {HLA_names[i]: pd.read_csv(_maptable[HLA_names[i]], '\s+', header=[0,1,2]) for i in range(len(HLA_names))}

    # for k, v in d_df_maptable.items():
    #     print("HLA : {}\n{}\n".format(k, v.head()))



    ##### < (1) Trimming AA dictionary > #####

    sr_Allele_AA = df_AA.loc[:, 'Allele']
    # print("sr_Allele_AA :\n{}\n".format(sr_Allele_AA.head()))

    if sr_Allele_AA.str.match(p_Field).all():

        df_AA['Allele_cutted'] = sr_Allele_AA.map(lambda x : p_Field.match(x).group())
        df_AA_Nfield = df_AA.loc[:, ['Allele_cutted', 'Seq']].drop_duplicates(['Allele_cutted'])
        # print("df_AA_Nfield :\n{}\n".format(df_AA_Nfield.head(70)))

        if _out:
            df_AA_Nfield.to_csv(_out+'.AA.{}field.txt'.format(_OUTPUT_FORMAT), sep='\t', header=False, index=False)
        else:
            df_AA_Nfield.to_csv(_AA_dict, sep='\t', header=False, index=False) # Override.

    else:
        print("Wrong allele names in Amino acid dictionary({})".format(_AA_dict))
        return -1



    ##### < (2) Trimming SNP dictionary > #####

    sr_Allele_SNP = df_SNP.loc[:, 'Allele']

    if sr_Allele_SNP.str.match(p_Field).all():

        df_SNP['Allele_cutted'] = sr_Allele_SNP.map(lambda x : p_Field.match(x).group())
        df_SNP_Nfield = df_SNP.loc[:, ['Allele_cutted', 'Seq']].drop_duplicates(['Allele_cutted'])
        # print("df_SNP_Nfield :\n{}\n".format(df_SNP_Nfield.head(70)))

        if _out:
            df_SNP_Nfield.to_csv(_out+'.SNP.{}field.txt'.format(_OUTPUT_FORMAT), sep='\t', header=False, index=False)
        else:
            df_SNP_Nfield.to_csv(_SNP_dict, sep='\t', header=False, index=False) # Override.

    else:
        print("Wrong allele names in DNA sequence(SNP) dictionary({})".format(_AA_dict))
        return -1



    ##### < (3) Trimming Maptable > #####

    _maptable_new = {}

    # for i in range(1):
    for i in range(len(HLA_names)):

        # print("HLA : {}\n{}".format(HLA_names[i], d_df_maptable[HLA_names[i]].head()))

        sr_Allele = d_df_maptable[HLA_names[i]].iloc[:, 0]
        # print(sr_Allele.head())

        if sr_Allele.str.match(p_Field).all():

            sr_Allele_cutted = sr_Allele.map(lambda x : p_Field.match(x).group())
            f_dup = sr_Allele_cutted.duplicated()
            # print("sr_Allele_cutted :\n{}".format(sr_Allele_cutted))
            df_maptable = pd.concat([sr_Allele_cutted, d_df_maptable[HLA_names[i]].iloc[:, 1:]], axis=1).loc[~f_dup, :]
            # print(df_maptable.head(70))

            if _out:
                df_maptable.to_csv(_out+'.maptable_{}.{}field.txt'.format(HLA_names[i], _OUTPUT_FORMAT), sep='\t',
                                   header=True, index=False)
                _maptable_new[HLA_names[i]] = _out+'.maptable_{}.{}field.txt'.format(HLA_names[i], _OUTPUT_FORMAT)
            else:
                df_maptable.to_csv(_maptable[HLA_names[i]], sep='\t', header=True, index=False)



    if _out:
        return [_out+'.AA.{}field.txt'.format(_OUTPUT_FORMAT), _out+'.SNP.{}field.txt'.format(_OUTPUT_FORMAT), _maptable_new]
    else:
        return [_AA_dict, _SNP_dict, _maptable]






if __name__ == '__main__':

    """

    < NfieldDictionary.py >


    - 

    """

    [_AA_dict, _SNP_dict, _maptable, _OUTPUT_FORMAT, _out] = sys.argv[1:]

    ## Ubuntu
    # _AA_dict = '/home/wanson/Git_Projects/HATK/tests/IMGT3370_test/HLA_DICTIONARY_AA.hg18.imgt3370.txt'
    # _SNP_dict = '/home/wanson/Git_Projects/HATK/tests/IMGT3370_test/HLA_DICTIONARY_SNPS.hg18.imgt3370.txt'
    # _maptable = {HLA_names[i]: '/home/wanson/Git_Projects/HATK/tests/IMGT3370_test/HLA_MAPTABLE_{}.hg18.imgt3370.txt'.format(HLA_names[i]) for i in range(len(HLA_names))}
    # _OUTPUT_FORMAT =3

    ## OS X
    # _AA_dict = '/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190924/HLA_DICTIONARY_AA.hg18.imgt3370.txt'
    # _SNP_dict = '/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190924/HLA_DICTIONARY_SNPS.hg18.imgt3370.txt'
    # _maptable = {HLA_names[i]: '/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190924/HLA_MAPTABLE_{}.hg18.imgt3370.txt'.format(HLA_names[i]) for i in range(len(HLA_names))}
    # _out = '/Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190924/NfieldTest'
    # _OUTPUT_FORMAT =3



    NfieldDictionary(_AA_dict, _SNP_dict, _maptable, _OUTPUT_FORMAT, _out=_out)