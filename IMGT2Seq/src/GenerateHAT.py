# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd

p_GenePrefix = re.compile(r'^(\w+\*)')


def GenerateHAT(_Nomenclature_2009, _Allelelist, _G_group, _P_group, _imgt, _out):


    _out = _out+'.imgt{}'.format(_imgt)
    OUTPUT_dir = os.path.dirname(_out)


    ### (1) Reading 'Nomenclature_2009.txt' file.
    __ALLELE_2009__ = pd.read_csv(_Nomenclature_2009, sep='\s+', header=None, dtype=str, skiprows=[0,1],
                                  names=['OLD_Allele', 'Allele'], na_values='None').dropna()
    # print("__ALLELE_2009__ :\n{}\n".format(__ALLELE_2009__.head(30)))


    ### (2) Reading 'Allelelist.txt'
    __ALLELE_LIST__ = pd.read_csv(_Allelelist, sep=',', header=0, comment='#', dtype=str)
    # print("__ALLELE_LIST__ :\n{}\n".format(__ALLELE_LIST__.head()))


    ### (3) Reading 'wmda/hla_nom_g.txt'
    __ALLELE_Ggroup__ = ExpandItems(pd.read_csv(_G_group, sep=';', dtype=str, header=None, comment='#').fillna('0'), 'Ggroup')
    # print("__ALLELE_Ggroup__ :\n{}\n".format(__ALLELE_Ggroup__.head()))

    ### (4) Reading 'wmda/hla_nom_p.txt'
    __ALLELE_Pgroup__ = ExpandItems(pd.read_csv(_P_group, sep=';', dtype=str, header=None, comment='#').fillna('0'), 'Pgroup')
    # print("__ALLELE_Pgroup__ :\n{}\n".format(__ALLELE_Pgroup__.head()))




    ## Merging above three information.

    df_merge0 = pd.merge( __ALLELE_LIST__, __ALLELE_2009__, left_on='Allele', right_on='Allele', how='left').fillna('0')
    # df_merge0.to_csv('/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/allele+old.txt', sep='\t', header=True, index=True)
    # print("df_merge0:\n{}\n".format(df_merge0))


    df_merge1 = pd.merge(df_merge0, __ALLELE_Ggroup__, how='left', right_on='Allele', left_on='Allele').fillna('0')
    # print("df_merge1:\n{}\n".format(df_merge1))
    # df_merge1.to_csv('/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/allele+old+group.txt', sep='\t', header=True, index=True)


    df_merge2 = pd.merge(df_merge1, __ALLELE_Pgroup__, how='left', right_on='Allele', left_on='Allele').fillna('0')
    # print("df_merge2:\n{}\n".format(df_merge2))
    # df_merge2.to_csv('/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/allele+old+group+pgroup.txt', sep='\t', header=True, index=True)


    # Extracting HLA gene name string.
    sr_HLA = df_merge2.loc[:, 'Allele'].str.extract(r'^(\w+)\*', expand=False)
    sr_HLA.name = 'HLA'
    # print(sr_HLA)


    __HAT__ = pd.concat([sr_HLA, df_merge2], axis=1)
    __HAT__ = __HAT__.applymap(lambda x: p_GenePrefix.sub('', x))
    __HAT__.columns = ['HLA', 'AlleleID', 'STANDARD', 'OLD', 'Ggroup', 'Pgroup']

    # print(__HAT__.head())
    __HAT__.to_csv(_out+'.hat', sep='\t', header=True, index=False)

    return _out+'.hat'






def ExpandItems(_df, _type):


    l_OUTPUT = []

    count = 0

    for row in _df.itertuples():

        chunk = re.split(string=row[2], pattern='/')

        if len(chunk) > 1:
            ### Expading items to rows
            l_OUTPUT.extend([[row[1]+eachAllele, row[3]] for eachAllele in chunk])

        else:
            l_OUTPUT.append([row[1]+row[2], row[3]])

        count += 1
        # if count > 5 : break

    df_OUTPUT = pd.DataFrame(l_OUTPUT)
    df_OUTPUT.columns = ['Allele', _type]
    # print(df_OUTPUT.head())

    return df_OUTPUT




if __name__ == '__main__':

    """
    
    < GenerateHAT.py >
    
    
    - Generate *.hat file to be used in 'NomenCleaner.py'.
    
    - (Necessary files in the IPD-IMGT/HLA database folder)
        1. Allelelist.txt
        2. wmda/hla_nom_g.txt
        3. wmda/hla_nom_p.txt
    
    """

    [_Nomenclature_2009, _Allelelist, _G_group, _P_group, _imgt, _out] = sys.argv[1:]

    # os x
    # _Nomenclature_2009 = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/Nomenclature_2009.txt'
    # _Allelelist = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/Allelelist.txt'
    # _G_group = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/wmda/hla_nom_g.txt'
    # _P_group = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/wmda/hla_nom_p.txt'
    # _imgt = '3320'
    # _out = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/HLA_ALLELE_TABLE'



    GenerateHAT(_Nomenclature_2009, _Allelelist, _G_group, _P_group, _imgt, _out)