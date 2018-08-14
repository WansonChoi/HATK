# -*- coding: utf-8 -*-


import pandas as pd
import argparse, textwrap




def MakeMap(_inputfile, _TYPE, _HLA, _OUTPUT = "output",
            _return_as_dataframe=False):

    """
    """

    ########## < Core Varialbes > ##########

    HLA_names = ("A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1")
    isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}

    # Final output which is exactly map file contents.
    df_MAPFILE = pd.DataFrame()
    df_forRETURN = pd.DataFrame()

    # (2018.2.18)
    forMap = pd.DataFrame()
    map_entities = []




    forMap = pd.read_table(_inputfile, sep='\t', header=None)

    POS = forMap.iloc[:, 0].tolist()

    Chr = ['6' for i in range(0, len(POS))]
    GM = ['0' for i in range(0, len(POS))]


    # print("Making MapFile for AA")
    # print(forMap.head())

    for i in range(0, len(forMap)):

        """
        col 0 : genomic_position
        col 1 : relative_position
        col 2 : codon value
        col 3 : AA value
        """

        genomic_position = forMap.iat[i, 0]
        relative_position = forMap.iat[i, 1]

        if genomic_position != 'i' and relative_position != 'i':

            map_entities.append('_'.join([("AA" if _TYPE == "AA" else "SNP"), _HLA, str(relative_position), str(genomic_position)]))

        else:
            # Indel part.
            # It is assumed that indel won't appear in 'i == 0' case.
            POS[i] = str(int((int(forMap.iat[i - 1, 0]) + int(forMap.iat[i + 1, 0])) / 2))
            map_entities.append('_'.join(["INDEL", _HLA, forMap.iat[i - 1, 1] + 'x' + forMap.iat[i + 1, 1], str(int((int(forMap.iat[i - 1, 0]) + int(forMap.iat[i + 1, 0])) / 2))]))



    # print(Chr)
    # print(len(Chr))
    # print(map_entities)
    # print(len(map_entities))
    # print(GM)
    # print(len(GM))
    # print(POS)
    # print(len(POS))


    df_MAPFILE = pd.DataFrame.from_dict({"Chr" : Chr,
                                         "Label" : map_entities,
                                         "GM" : GM,
                                         "POS" : POS},
                                        orient='columns')

    df_MAPFILE = df_MAPFILE[["Chr", "Label", "GM", "POS"]]


    if _return_as_dataframe:

        # Just return the DataFrame which has the result of Mapfile.
        return df_MAPFILE

    else:

        if _TYPE == "AA":
            df_MAPFILE.to_csv(_OUTPUT+'.AA.map', sep='\t', header=False, index=False)
        elif _TYPE == "SNPS":
            df_MAPFILE.to_csv(_OUTPUT + '.SNPS.map', sep='\t', header=False, index=False)


        return 0







if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        MakeMap.py


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')
    parser.add_argument("-i", help="Input file name\n\n", required=True, metavar='INPUT')
    parser.add_argument("-o", help="Output file name prefix\n\n", required=True, metavar='OUTPUT')
    parser.add_argument("--type", help="Sequence type to deal with(Amino Acids[AA] or SNPs[SNPS]\n\n", required=True, choices=["AA", "SNPS"], metavar='TYPE')
    parser.add_argument("--HLA", help="HLA type\n\n", required=True, choices=["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"], metavar='HLA')



    ##### < for Test > #####

    # args = parser.parse_args(["-i", "./IMGT_A.AA.forMAP.txt", "-o", "IMGT_A", "--type", "AA", "--HLA", "A"])
    # args = parser.parse_args(["-i", "./IMGT_C.AA.forMAP.txt", "-o", "IMGT_C", "--type", "AA", "--HLA", "C"])
    # args = parser.parse_args(["-i", "./IMGT_B.AA.forMAP.txt", "-o", "IMGT_B", "--type", "AA", "--HLA", "B"])
    # args = parser.parse_args(["-i", "./IMGT_DRB.AA.forMAP.txt", "-o", "IMGT_DRB", "--type", "AA", "--HLA", "DRB1"])
    # args = parser.parse_args(["-i", "./IMGT_DQA.AA.forMAP.txt", "-o", "IMGT_DQA", "--type", "AA", "--HLA", "DQA1"])
    # args = parser.parse_args(["-i", "./IMGT_DQB.AA.forMAP.txt", "-o", "IMGT_DQB", "--type", "AA", "--HLA", "DQB1"])
    # args = parser.parse_args(["-i", "./IMGT_DPA.AA.forMAP.txt", "-o", "IMGT_DPA", "--type", "AA", "--HLA", "DPA1"])
    # args = parser.parse_args(["-i", "./IMGT_DPB.AA.forMAP.txt", "-o", "IMGT_DPB", "--type", "AA", "--HLA", "DPB1"])

    # args = parser.parse_args(["-i", "./IMGT_A.SNPS.forMAP.txt", "-o", "IMGT_A", "--type", "SNPS", "--HLA", "A"])
    # args = parser.parse_args(["-i", "./IMGT_C.SNPS.forMAP.txt", "-o", "IMGT_C", "--type", "SNPS", "--HLA", "C"])
    # args = parser.parse_args(["-i", "./IMGT_B.SNPS.forMAP.txt", "-o", "IMGT_B", "--type", "SNPS", "--HLA", "B"])
    # args = parser.parse_args(["-i", "./IMGT_DRB.SNPS.forMAP.txt", "-o", "IMGT_DRB", "--type", "SNPS", "--HLA", "DRB1"])
    # args = parser.parse_args(["-i", "./IMGT_DQA.SNPS.forMAP.txt", "-o", "IMGT_DQA", "--type", "SNPS", "--HLA", "DQA1"])
    # args = parser.parse_args(["-i", "./IMGT_DQB.SNPS.forMAP.txt", "-o", "IMGT_DQB", "--type", "SNPS", "--HLA", "DQB1"])
    # args = parser.parse_args(["-i", "./IMGT_DPA.SNPS.forMAP.txt", "-o", "IMGT_DPA", "--type", "SNPS", "--HLA", "DPA1"])
    # args = parser.parse_args(["-i", "./IMGT_DPB.SNPS.forMAP.txt", "-o", "IMGT_DPB", "--type", "SNPS", "--HLA", "DPB1"])


    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    MakeMap(args.i, args.o, args.type, args.HLA)