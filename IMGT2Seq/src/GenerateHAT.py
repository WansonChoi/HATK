# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

p_GenePrefix = re.compile(r'^(\w+\*)')


def GenerateHAT(_Nomenclature_2009, _Allelelist, _G_group, _P_group, _imgt, _out,
                _f_Answer=True):


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


    if _f_Answer:

        ### Generating Answer file to test NomenCleaner.py
        # print("### Generating Answer file to test NomenCleaner.py")

        # Only to 9 HLAs in `HLA_names`.
        __HAT__.set_index('HLA', inplace=True, drop=False)
        d__HAT__ = {HLA_names[i] : __HAT__.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}

        # for k, v in d__HAT__.items():
        #     print("HLA : {}\n{}\n".format(k,v.head()))


        l__HAT_4field__ = []
        l__HAT_1field__ = []
        l__HAT_2field__ = []
        l__HAT_3field__ = []

        for i in range(0, len(HLA_names)):
        # for i in range(0, 1):

            t_sr_Allele = d__HAT__[HLA_names[i]].loc[:, 'Allele']
            # print("t_sr_Allele : \n{}".format(t_sr_Allele.head()))


            ### Preparing columns of 1,2,3-field

            # 1-field
            t_sr_Allele_1field = t_sr_Allele.str.extract(r'^(\d{2,3})').fillna('0')
            # print("sr_Allele_1field : \n{}".format(t_sr_Allele_1field.head(20)))

            # 2-field
            t_sr_Allele_2field = t_sr_Allele.str.extract(r'^(\d{2,3}:\d{2,3}[A-Z]?)').fillna('0')
            # print("sr_Allele_2field : \n{}".format(t_sr_Allele_2field.head(20)))

            # 3-field
            t_sr_Allele_3field = t_sr_Allele.str.extract(r'^(\d{2,3}:\d{2,3}:\d{2,3}[A-Z]?)').fillna('0')
            # print("sr_Allele_3field : \n{}".format(t_sr_Allele_3field.head(20)))



            df_Alleles = pd.concat([t_sr_Allele_1field, t_sr_Allele_2field, t_sr_Allele_3field, t_sr_Allele], axis=1)
            df_Alleles.columns = ['1field', '2field', '3field', 'Allele']

            df_Alleles = pd.concat([d__HAT__[HLA_names[i]].iloc[:, :2], df_Alleles, d__HAT__[HLA_names[i]].iloc[:, 3:]], axis=1)
            # print("df_Allele : \n{}".format(df_Alleles.head(30)))


            ### Gathering 1st item.

            # 1-field
            df_1field_1stItems = Gather1stRows(df_Alleles, '1field')
            # print("df_1field_1stItems : \n{}".format(df_1field_1stItems.head()))

            # 2-field
            df_2field_1stItems = Gather1stRows(df_Alleles, '2field')
            # print("df_2field_1stItems : \n{}".format(df_2field_1stItems.head()))

            # 3-field
            df_3field_1stItems = Gather1stRows(df_Alleles, '3field')
            # print("df_3field_1stItems : \n{}".format(df_3field_1stItems.head()))



            l__HAT_4field__.append(df_Alleles)
            l__HAT_1field__.append(df_1field_1stItems)
            l__HAT_2field__.append(df_2field_1stItems)
            l__HAT_3field__.append(df_3field_1stItems)


        # Alleles
        __HAT_4field__ = pd.concat(l__HAT_4field__, axis=0)
        __HAT_4field__.to_csv(_out+'.4field.answer.hat', sep='\t', header=True, index=False)

        # 1field
        __HAT_1field__ = pd.concat(l__HAT_1field__, axis=0)
        __HAT_1field__.to_csv(_out+'.1field.answer.hat', sep='\t', header=True, index=False)

        # 2field
        __HAT_2field__ = pd.concat(l__HAT_2field__, axis=0)
        __HAT_2field__.to_csv(_out+'.2field.answer.hat', sep='\t', header=True, index=False)

        # 3field
        __HAT_3field__ = pd.concat(l__HAT_3field__, axis=0)
        __HAT_3field__.to_csv(_out+'.3field.answer.hat', sep='\t', header=True, index=False)



    return __HAT__






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



def Gather1stRows(_df, _column):


    sr_Allele = _df.loc[:, _column]

    l_1stItems = []

    for item in sr_Allele.unique():

        flag_item = (sr_Allele == item)

        if flag_item.any():

            df_temp = _df.loc[sr_Allele == item, :].iloc[0, :].to_frame().transpose()
            # print("df_temp : \n{}".format(df_temp))
            # print(type(df_temp))

            l_1stItems.append(df_temp)

    df_1stItem = pd.concat(l_1stItems, axis=0)
    # print("df_1stItem : \n{}".format(df_1stItem))

    return df_1stItem









if __name__ == '__main__':

    """
    
    < GenerateHAT.py >
    
    
    - Generate *.hat file to be used in 'NomenCleaner.py'.
    
    - (Necessary files in the IPD-IMGT/HLA database folder)
        1. Allelelist.txt
        2. wmda/hla_nom_g.txt
        3. wmda/hla_nom_p.txt
    
    """

    # [_Allelelist, _hg, _out] = sys.argv[1:]

    # os x
    _Nomenclature_2009 = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/Nomenclature_2009.txt'
    _Allelelist = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/Allelelist.txt'
    _G_group = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/wmda/hla_nom_g.txt'
    _P_group = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/IMGTHLA/IMGTHLA3320/wmda/hla_nom_p.txt'
    _imgt = '3320'
    _out = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/HLA_ALLELE_TABLE'



    GenerateHAT(_Nomenclature_2009, _Allelelist, _G_group, _P_group, _imgt, _out, _f_Answer=False)
    # GenerateHAT(_Allelelist, _G_group, _P_group, _imgt, _out, _f_Answer=True)