# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

def XgroupDictionary(_hat, _AA_dict, _SNP_dict, _maptable, _out, _OUTPUT_FORMAT=5,
                     HLA_names = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')):

    Xgroup = None

    if _OUTPUT_FORMAT == 5:  # Ggroup

        Xgroup = 'Ggroup'

    elif _OUTPUT_FORMAT == 6:  # Pgroup

        Xgroup = 'Pgroup'

    """

    1. Get 'df_1st_item'.
    2. Modify AA dict.
    3. Modify SNPS dict.
    4. Modify Maptable.
    5. Output File write.

    """

    ### < 1. Get 'df_1st_item'. > ###
    (df_1st_item, df_1st_item2) = get1stItem(_hat, Xgroup)
    # print("df_1st_item2:\n{}\n".format(df_1st_item2))

    ### < 2. Modify AA dict. > ###
    df_AA_dict_Xgroup = Subset_Dict(_AA_dict, df_1st_item2)
    # print("df_AA_dict_Xgroup:\n{}\n".format(df_AA_dict_Xgroup))

    ### < 3. Modify SNPS dict. > ###
    df_SNPS_dict_Xgroup = Subset_Dict(_SNP_dict, df_1st_item2)
    # print("df_SNPS_dict_Xgroup:\n{}\n".format(df_SNPS_dict_Xgroup))

    ### < 4. Modify Maptable. > ###
    dict_Ggroup_MapTable_HLA = {
        hla: \
            Subset_Maptable(_maptable[hla], df_1st_item2) \
        for hla in HLA_names
    }



    ### < 5. Output File write. > ###
    if _out:
        # Output prefix specified.
        df_AA_dict_Xgroup.to_csv(_out + '.HLA_DICTIONARY.AA.{}.txt'.format(Xgroup), sep='\t', header=False, index=False)
        df_SNPS_dict_Xgroup.to_csv(_out + '.HLA_DICTIONARY.SNPS.{}.txt'.format(Xgroup), sep='\t', header=False,
                                   index=False)

        for hla in HLA_names:
            dict_Ggroup_MapTable_HLA[hla][0].to_csv(_out + '.HLA_MAPTABLE.{}.{}.txt'.format(hla, Xgroup), sep='\t',
                                                    header=True, index=True)

    else:
        # Override.
        df_AA_dict_Xgroup.to_csv(_AA_dict, sep='\t', header=False, index=False)
        df_SNPS_dict_Xgroup.to_csv(_SNP_dict, sep='\t', header=False, index=False)

        for hla in HLA_names:
            dict_Ggroup_MapTable_HLA[hla][0].to_csv(_maptable[hla], sep='\t', header=True, index=True)


    return 0





def get1stItem(_hat, Xgroup, HLA_names = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')):

    hat = pd.read_csv(_hat, sep='\s+', header=0, dtype=str)
    #     print("hat:\n{}\n".format(hat))

    """
    have to do the gropuby with HLA first,
    because there are same G-group labels in different HLAs.

    (ex.)
    A*01:01:01G
    DPB1*01:01:01G

    """

    hat_HLA = hat.groupby('HLA')
    dict_hat_HLA = dict(list(hat_HLA))
    #     print(dict_hat_HLA['A'])

    ### Main iteration ###

    l_df_HLA = []

    #     for hla in ['A']:
    for hla in HLA_names:

        df_HLA = dict_hat_HLA[hla]
        #         print("df_HLA:\n{}\n".format(df_HLA))

        l_HLA_Ggroup = []

        #     df_Ggroup = df_HLA.groupby('Ggroup')
        for _Ggroup, _df in df_HLA.groupby(Xgroup):

            if _Ggroup == '0': continue

            l_HLA_Ggroup.append(_df.drop_duplicates(Xgroup))  # Gathering the 1st item.

        df_temp = pd.concat(l_HLA_Ggroup)
        l_df_HLA.append(df_temp)

    df_1st_item = pd.concat(l_df_HLA)  ### return 1
    #     print("df_1st_item:\n{}\n".format(df_1st_item))

    df_4field = df_1st_item[['HLA', 'STANDARD']].apply(lambda x: '*'.join(x), axis=1)
    #     print(df_4field)
    df_Xgroup = df_1st_item[['HLA', Xgroup]].apply(lambda x: '*'.join(x), axis=1)
    #     print(df_Xgroup)

    df_1st_item2 = pd.concat([df_4field, df_Xgroup], axis=1)
    df_1st_item2.columns = ['4field', Xgroup]

    return (df_1st_item, df_1st_item2)



def Subset_Dict(_dict, df_1st_item2):

    df_dict = pd.read_csv(_dict, sep='\s+', header=None, dtype=str, names=['Allele', 'Seq'])
    #     print("df_dict:\n{}\n".format(df_dict))

    df_merge0 = df_1st_item2.merge(df_dict, left_on='4field', right_on='Allele')
    #     print("df_merge0:\n{}\n".format(df_merge0))

    df_dict_Xgruop = df_merge0.iloc[:, [1, 3]]
    #     print("df_dict_Xgruop:\n{}\n".format(df_dict_Xgruop))

    return df_dict_Xgruop



def Subset_Maptable(_maptable_HLA, df_1st_item2):

    df_maptable_HLA = pd.read_csv(_maptable_HLA, sep='\s+', header=[0, 1, 2], dtype=str, index_col=0)
    #     print("df_maptable_HLA:\n{}\n".format(df_maptable_HLA))

    #     print("df_1st_item2:\n{}\n".format(df_1st_item2))

    sr_index = df_maptable_HLA.index.to_series(name='Allele')
    #     print(sr_index)

    ### (Optional) Subset table.
    df_subset_table = df_maptable_HLA.index.to_frame(name='Allele').merge(df_1st_item2, left_on='Allele',
                                                                          right_on='4field', how='left')
    #     print("df_subset_table:\n{}\n".format(df_subset_table))

    """
    (1) 
    df_merge0 = df_maptable_HLA.index.to_frame(name='Allele').merge(df_1st_item2, left_on='Allele', right_on='4field')
    print("df_merge0:\n{}\n".format(df_merge0))
    print(df_merge0.shape) # 52

    (2) f_sub = sr_index.isin(df_1st_item2['4field']) # 52

    I checked both give the same result.
    """

    ### Subset

    f_sub = sr_index.isin(df_1st_item2['4field'])
    df_maptable_HLA = df_maptable_HLA[f_sub]
    #     print("df_maptable_HLA (Subsetted):\n{}\n".format(df_maptable_HLA))

    ### New index as Xgroup

    dict_Xgroup = df_1st_item2.set_index('4field').iloc[:, 0].to_dict()
    #     print(dict_Xgroup)

    sr_index_new = df_maptable_HLA.index.to_series().map(dict_Xgroup)  # Transform as Xgroup label.
    #     print("sr_index_new:\n{}\n".format(sr_index_new))

    df_maptable_HLA.index = sr_index_new
    #     print("df_maptable_HLA (Xgroup):\n{}\n".format(df_maptable_HLA))

    return (df_maptable_HLA, df_subset_table)