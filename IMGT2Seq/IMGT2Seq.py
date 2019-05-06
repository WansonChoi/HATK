# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd
import argparse, textwrap
from glob import glob
import multiprocessing as mp


########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


class IMGT_Dictionary(object):

    def __init__(self, _imgt, _hg, _out, *args, **kwargs):

        """

        Either (1) generate a new Dictionary files or (2) find existing dictionaries.

        """

        ########## < Assigning arguments > ##########

        self.imgt = _imgt
        self.hg = _hg
        self.out = _out

        # Main Outputs
        self.__dict_AA__ = None
        self.__dict_SNPS__ = None
        self.__IAT__ = None
        self.__d_MapTable__ = None

        # Flags to check
        self.f__dict_AA__ = False
        self.f__dict_SNPS__ = False
        self.f__IAT__ = False
        self.f__d_MapTable__ = False

        # Optional arguments
        self.no_indel = kwargs["_no_indel"]
        self.multiprocess = kwargs["_multiprocess"]
        self.save_intermediates = kwargs["_save_intermediates"]
        self.imgt_dir = kwargs["_imgt_dir"]

        # Summary string
        self.summary_string = ""


        dir_path = os.path.dirname(self.out)
        version_label = "hg{}.imgt{}".format(self.hg, self.imgt)




        ### Checking results
        # print(std_MAIN_PROCESS_NAME + "Arguments for HLA Dictionary")
        # for k, v in self.__dict__.items():
        #     print("{} : {}".format(k, v))




        ########## < Checking existing results > ##########

        # (1) dict_AA
        dict_AA_txt = os.path.join(dir_path, "HLA_DICTIONARY_AA.{}.txt".format(version_label))
        dict_AA_map = os.path.join(dir_path, "HLA_DICTIONARY_AA.{}.map".format(version_label))

        if os.path.exists(dict_AA_txt) and os.path.exists(dict_AA_map):
            self.__dict_AA__ = os.path.join(dir_path, "HLA_DICTIONARY_AA.{}".format(version_label))
            self.f__dict_AA__ = True

        # (2) dict_SNPS
        dict_SNPS_txt = os.path.join(dir_path, "HLA_DICTIONARY_SNPS.{}.txt".format(version_label))
        dict_SNPS_map = os.path.join(dir_path, "HLA_DICTIONARY_SNPS.{}.map".format(version_label))

        if os.path.exists(dict_SNPS_txt) and os.path.exists(dict_SNPS_map):
            self.__dict_SNPS__ = os.path.join(dir_path, "HLA_DICTIONARY_SNPS.{}".format(version_label))
            self.f__dict_SNPS__ = True

        # (3) iat
        iat = os.path.join(dir_path, "HLA_INTEGRATED_ALLELE_TABLE.{}.iat".format(version_label))

        if os.path.exists(iat):
            self.__IAT__ = iat
            self.f__IAT__ = True

        # (4) maptables
        self.__d_MapTable__ = {hla_name: None for hla_name in HLA_names}

        f_temp = True

        for hla_name in HLA_names:
            t_maptable = os.path.join(dir_path, "HLA_MAPTABLE_{}.{}.txt".format(hla_name, version_label))
            if os.path.exists(t_maptable):
                self.__d_MapTable__[hla_name] = t_maptable
            else:
                f_temp = False

        self.f__d_MapTable__ = f_temp




        # if (self.f__dict_AA__ and self.f__dict_SNPS__ and self.f__IAT__ and self.f__d_MapTable__):
        #
        #     ### All of previously generated dictionary files are available.
        #
        #     self.summary_string = \
        #         ''.join([self.summary_string,
        #                  "< IMGT2Sequence(Using existing results) >\n" \
        #                  "- HLA Amino Acids : {}\n" \
        #                  "- HLA SNPs : {}\n" \
        #                  "- Integrated Allele Table : {}\n" \
        #                  "- Maptables for heatmap : \n" \
        #                  "   A   : {A}\n" \
        #                  "   B   : {B}\n" \
        #                  "   C   : {C}\n" \
        #                  "   DPA1: {DPA1}\n" \
        #                  "   DPB1: {DPB1}\n" \
        #                  "   DQA1: {DQA1}\n" \
        #                  "   DQB1: {DQB1}\n" \
        #                  "   DRB1: {DRB1}\n".format(self.__dict_AA__, self.__dict_SNPS__, self.__IAT__, **self.__d_MapTable__)
        #                  ])
        #
        # else:
        #
        #     ### Not all of previously generated dictionary files are available.
        #
        #     self.__dict_AA__, self.__dict_SNPS__, self.__IAT__, self.__d_MapTable__ = \
        #         IMGT2Sequence(self.imgt, self.hg, self.out,
        #                       _no_Indel=self.no_indel, _MultiP=self.multiprocess,
        #                       _save_intermediates=self.save_intermediates, _imgt_dir=self.imgt_dir)
        #
        #
        #     self.f__dict_AA__ = True
        #     self.f__dict_SNPS__ = True
        #     self.f__IAT__ = True
        #     self.f__d_MapTable__ = True
        #
        #
        #     self.summary_string = \
        #         ''.join([self.summary_string,
        #                  "< IMGT2Sequence(Newly generated.) >\n" \
        #                  "- HLA Amino Acids : {}\n" \
        #                  "- HLA SNPs : {}\n" \
        #                  "- Integrated Allele Table : {}\n" \
        #                  "- Maptables for heatmap : \n" \
        #                  "   A   : {A}\n" \
        #                  "   B   : {B}\n" \
        #                  "   C   : {C}\n" \
        #                  "   DPA1: {DPA1}\n" \
        #                  "   DPB1: {DPB1}\n" \
        #                  "   DQA1: {DQA1}\n" \
        #                  "   DQB1: {DQB1}\n" \
        #                  "   DRB1: {DRB1}\n".format(self.__dict_AA__, self.__dict_SNPS__, self.__IAT__, **self.__d_MapTable__)
        #                  ])


        # temporary hardcoding ---
        self.__dict_AA__, self.__dict_SNPS__, self.__IAT__, self.__d_MapTable__ = \
            IMGT2Seq(self.imgt, self.hg, self.out, _no_Indel=self.no_indel, _MultiP=self.multiprocess,
                     _save_intermediates=self.save_intermediates, _imgt_dir=self.imgt_dir)


        self.f__dict_AA__ = True
        self.f__dict_SNPS__ = True
        self.f__IAT__ = True
        self.f__d_MapTable__ = True


        self.summary_string = \
            ''.join([self.summary_string,
                     "< IMGT2Sequence(Newly generated.) >\n" \
                     "- HLA Amino Acids : {}\n" \
                     "- HLA SNPs : {}\n" \
                     "- Integrated Allele Table : {}\n" \
                     "- Maptables for heatmap : \n" \
                     "   A   : {A}\n" \
                     "   B   : {B}\n" \
                     "   C   : {C}\n" \
                     "   DPA1: {DPA1}\n" \
                     "   DPB1: {DPB1}\n" \
                     "   DQA1: {DQA1}\n" \
                     "   DQB1: {DQB1}\n" \
                     "   DRB1: {DRB1}\n".format(self.__dict_AA__, self.__dict_SNPS__, self.__IAT__, **self.__d_MapTable__)
                     ])
        # --Hardcoding#


    def __bool__(self):
        return (self.f__dict_AA__ and self.f__dict_SNPS__ and self.f__IAT__ and self.f__d_MapTable__)


    def __str__(self):
        return self.summary_string


    def getDictionaries(self):
        return [self.__dict_AA__, self.__dict_SNPS__, self.__IAT__, self.__d_MapTable__]







def IMGT2Seq(_imgt, _hg, _out, _no_Indel=False, _MultiP=False, _save_intermediates=False, _imgt_dir=None,
             _no_prime = True, _p_data="./data", _p_src="./src"):

    """
    """

    ########## < Core Variables > ##########


    # (2018. 8. 28.) This order must be kept.
    raw_HLA_names = ["A", "B", "C", "DPA", "DPB", "DQA", "DQB", "DRB"]
    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


    ### Dictionaries for Raw files.
    TARGET_prot_files = {}
    TARGET_nuc_files = {}
    TARGET_gen_files = {}


    ### Variables for Paths.
    p_data = _p_data
    p_src = _p_src


    ### raw MapTables
    d_MapTables = {}

    ### OUTPUT prefix

    # Preparing intermediate paths.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)):
        INTERMEDIATE_PATH = os.path.dirname(_out)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)


    _OUTPUT_AA_RETURN = os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_AA.hg{0}.imgt{1}'.format(_hg, _imgt))
    _OUTPUT_SNPS_RETURN = os.path.join(INTERMEDIATE_PATH, 'HLA_DICTIONARY_SNPS.hg{0}.imgt{1}'.format(_hg, _imgt))
    _OUTPUT_IAT_RETURN = os.path.join(INTERMEDIATE_PATH, 'HLA_INTEGRATED_ALLELE_TABLE.hg{0}.imgt{1}'.format(_hg, _imgt))


    ########## < Dependency Checking > ##########

    """
    1. Necessary Python source files.
    2. Necessary Static Data files
    3. Checking each files of "*_nuc.txt", "*_prot.txt", "*_gen.txt".
    """


    ##### < Necessary Python source files > #####

    # ProcessIMGT.py
    if os.path.exists(os.path.join(p_src, "ProcessIMGT.py")):
        from src.ProcessIMGT import ProcessIMGT
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "\"src/ProcessIMGT.py\" doesn't exist!")
        sys.exit()

    # ClassifyGroups.py
    if os.path.exists(os.path.join(p_src, "ClassifyGroups.py")):
        from src.ClassifyGroups import ClassifyGroups
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "\"src/ClassifyGroups.py\" doesn't exist!")
        sys.exit()



    ##### < Necessary Static Data files > #####

    IMGTHLA_directory = os.path.join(p_data, "IMGTHLA{0}".format(_imgt)) if not bool(_imgt_dir) else _imgt_dir
    # If `_imgt_dir` is not given, then use default IMGTHLA data in `p_data`.

    # print(IMGTHLA_directory)

    if not os.path.isdir(IMGTHLA_directory):
        print(std_ERROR_MAIN_PROCESS_NAME + "The given IMGT-HLA directory \"{0}\" can't be found!".format(IMGTHLA_directory))
        sys.exit()

    # (2018. 10. 26.) Loading `HLA_INTEGRATED_TABLE` file is moved to "ProcessIMGT.py".



    ##### < "*_nuc.txt", "*_prot.txt", "*_gen.txt" > #####

    files_in_alignments = glob(IMGTHLA_directory + '/alignments/*')

    for i in range(0, len(HLA_names)):

        ### <*_prot.txt, *_nuc.txt>
        found_prot = list(filter(re.compile("/" + raw_HLA_names[i] + "_prot.txt").search, files_in_alignments))
        found_nuc = list(filter(re.compile("/" + raw_HLA_names[i] + "_nuc.txt").search, files_in_alignments))

        if len(found_prot) == 1:
            # When found exactly.
            TARGET_prot_files[HLA_names[i]] = found_prot.pop()
            TARGET_nuc_files[HLA_names[i]] = found_nuc.pop()

        else:
            # If not, then try same job to the names with number "1".
            found_prot = list(filter(re.compile("/" + HLA_names[i] + "_prot.txt").search, files_in_alignments))
            found_nuc = list(filter(re.compile("/" + HLA_names[i] + "_nuc.txt").search, files_in_alignments))

            if len(found_prot) == 1:
                # When found exactly in 2nd trial.
                TARGET_prot_files[HLA_names[i]] = found_prot.pop()
                TARGET_nuc_files[HLA_names[i]] = found_nuc.pop()
            else:
                # Not found finally.
                print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}_prot.txt\" or \"{0}_nuc.txt\" file can't be found.".format(HLA_names[i]))
                sys.exit()

        ### <*_gen.txt>

        found_gen = list(filter(re.compile("/" + raw_HLA_names[i] + "_gen.txt").search, files_in_alignments))

        if len(found_gen) == 1:
            # When found exactly.
            TARGET_gen_files[HLA_names[i]] = found_gen.pop()
        else:
            # If not, then try same job to the names with number "1".
            found_gen = list(filter(re.compile("/" + HLA_names[i] + "_gen.txt").search, files_in_alignments))

            if len(found_gen) == 1:
                # When found exactly in 2nd trial.
                TARGET_gen_files[HLA_names[i]] = found_gen.pop()
            else:
                # Not found finally.
                print(std_ERROR_MAIN_PROCESS_NAME + "\"{0}_gen.txt\" file can't be found.".format(HLA_names[i]))
                sys.exit()

    # print("\nTarget prot and gen files : \n")
    # print(TARGET_prot_files)
    # print(TARGET_gen_files)


    ##### < "Allelelist.txt", "hla_nom_g.txt", "hla_nom_p.txt" > #####

    # (2018. 10. 12.)

    t_allelelist = os.path.join(IMGTHLA_directory, "Allelelist.txt")
    t_p_Group = os.path.join(IMGTHLA_directory, "wmda/hla_nom_g.txt")
    t_p_Proup = os.path.join(IMGTHLA_directory, "wmda/hla_nom_p.txt")

    if not os.path.exists(t_allelelist):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"Allelelist.txt\" file can't be found.\n")
        sys.exit()

    if not os.path.exists(t_p_Group):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"hla_nom_g.txt\" file can't be found.\n")
        sys.exit()

    if not os.path.exists(t_p_Proup):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"hla_nom_p.txt\" file can't be found.\n")
        sys.exit()




    ########## < Control Flags > ##########

    _1_MAKING_DICTIONARY = 1
    _2_CLASSIFY_GROUP = 1
    CLEAN_UP = 0


    if _1_MAKING_DICTIONARY:

        ########## < 1. Making HLA Dictionary. > ##########

        # print(std_MAIN_PROCESS_NAME + "[1] Making HLA dictionary file.")

        l_df_Seqs_AA = []
        l_df_forMAP_AA = []

        l_df_Seqs_SNPS = []
        l_df_forMAP_SNPS = []



        if _MultiP > 1:

            pool = mp.Pool(processes=_MultiP)

            dict_Pool = {HLA_names[i]: pool.apply_async(ProcessIMGT, (_out, HLA_names[i], _hg, _imgt,
                                                                      TARGET_nuc_files[HLA_names[i]], TARGET_gen_files[HLA_names[i]], TARGET_prot_files[HLA_names[i]],
                                                                      p_data, _no_Indel, _save_intermediates))
                         for i in range(0, len(HLA_names))}

            pool.close()
            pool.join()


        for i in range(0, len(HLA_names)):

            if _MultiP > 1:

                t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE = dict_Pool[HLA_names[i]].get()

            else:
                t_df_Seqs_SNPS, t_df_Seqs_AA, t_df_forMAP_SNPS, t_df_forMAP_AA, t_MAPTABLE \
                    = ProcessIMGT(_out, HLA_names[i], _hg, _imgt,
                                  TARGET_nuc_files[HLA_names[i]], TARGET_gen_files[HLA_names[i]], TARGET_prot_files[HLA_names[i]],
                                  p_data, _no_Indel=_no_Indel, _save_intermediates=_save_intermediates)



            # Dictionary AA.
            l_df_Seqs_AA.append(t_df_Seqs_AA)
            l_df_forMAP_AA.append(MakeMap(HLA_names[i], "AA", t_df_forMAP_AA))

            # Dictionary SNPS.
            l_df_Seqs_SNPS.append(t_df_Seqs_SNPS)
            l_df_forMAP_SNPS.append(MakeMap(HLA_names[i], "SNPS", t_df_forMAP_SNPS))

            # raw MapTable
            d_MapTables[HLA_names[i]] = t_MAPTABLE


        HLA_DICTIONARY_AA = pd.concat(l_df_Seqs_AA, axis=0)
        HLA_DICTIONARY_AA_map = pd.concat(l_df_forMAP_AA, axis=0)

        HLA_DICTIONARY_SNPS = pd.concat(l_df_Seqs_SNPS, axis=0)
        HLA_DICTIONARY_SNPS_map = pd.concat(l_df_forMAP_SNPS, axis=0)


        # Exporting AA dictionary.
        HLA_DICTIONARY_AA.to_csv(_OUTPUT_AA_RETURN + ".txt", sep='\t', header=False, index=True)

        HLA_DICTIONARY_AA_map.to_csv(_OUTPUT_AA_RETURN + ".map", sep='\t', header=False, index=False)

        # Exporting SNPS dictionary.
        HLA_DICTIONARY_SNPS.to_csv(_OUTPUT_SNPS_RETURN+".txt", sep='\t', header=False, index=True)

        HLA_DICTIONARY_SNPS_map.to_csv(_OUTPUT_SNPS_RETURN+".map", sep='\t', header=False, index=False)


    if _2_CLASSIFY_GROUP:


        ########## < 2. Making *.iat file. > ##########

        # print(std_MAIN_PROCESS_NAME + "[2] Making \"*.iat\" file.")

        t_dict_AA = _OUTPUT_AA_RETURN + ".txt"
        t_dict_SNPS = _OUTPUT_SNPS_RETURN+ ".txt"

        _OUTPUT_IAT_RETURN = ClassifyGroups(t_allelelist, t_p_Group, t_p_Proup, t_dict_AA, t_dict_SNPS, _OUTPUT_IAT_RETURN)




    if CLEAN_UP:

        ########## < 5. Removing unnecessary files. > ##########

        print(std_MAIN_PROCESS_NAME + "[5] Removing unnecessary files")

    # print("\n")
    # print(std_MAIN_PROCESS_NAME + "Making HLA Dictionary is done.\n")


    return [_OUTPUT_AA_RETURN, _OUTPUT_SNPS_RETURN, _OUTPUT_IAT_RETURN, d_MapTables]






#################### < Core Functions > ####################

def MakeMap(_hla, _type, _df_forMAP, _no_prime=True):

    if not (_type == "AA" or _type == "SNPS"):
        return -1

    # if _type == "AA" and (_df_forMAP.shape[1] != 3):
    #     return -1
    #
    # if _type == "SNPS" and (_df_forMAP.shape[1] != 6):
    #     return -1

    """
    (1) Chr
    (2) Label
    (3) GD
    (4) Genomic Positions
    """


    if not _no_prime:

        _df_forMAP = _df_forMAP.astype(str)

        sr_Chr = pd.Series(["6" for i in range(0, _df_forMAP.shape[0])])
        l_Label = []
        sr_GD = pd.Series(["0" for i in range(0, _df_forMAP.shape[0])])
        sr_GenPos = _df_forMAP.iloc[:, 1]

        # Processing `l_Label`.

        p = re.compile('-?\d+x-?\d+')

        for i in range(0, _df_forMAP.shape[0]):

            t_rel_pos = _df_forMAP.iat[i, 0]
            t_gen_pos = _df_forMAP.iat[i, 1]
            t_type = _df_forMAP.iat[i, 2]

            main_label = '_'.join([("INDEL" if bool(p.match(t_rel_pos)) else "AA" if _type == "AA" else "SNPS"),
                                   _hla, t_rel_pos, t_gen_pos, t_type])

            if _type == "SNPS":

                additional_label = ["AA"]

                t_AA_rel_pos = _df_forMAP.iat[i, 3]
                t_AA_gen_pos = _df_forMAP.iat[i, 4]
                #             t_AA_type = _df_forMAP.iat[i, 5]

                if t_AA_rel_pos != "nan" and t_AA_rel_pos != "NaN":
                    additional_label.append(t_AA_rel_pos)

                if t_AA_gen_pos != "nan" and t_AA_gen_pos != "NaN":
                    additional_label.append(t_AA_gen_pos)


                if len(additional_label) > 1:
                    additional_label = '_'.join(additional_label)
                    main_label = '_'.join([main_label, additional_label])

            l_Label.append(main_label)

        df_MAP = pd.concat([sr_Chr, pd.Series(l_Label), sr_GD, sr_GenPos], axis=1)




    else:

        sr_Chr = pd.Series(("6" for i in range(0, _df_forMAP.shape[0])))
        l_Label = []
        sr_GD = pd.Series(("0" for i in range(0, _df_forMAP.shape[0])))
        sr_GenPos = _df_forMAP.iloc[:, 1]


        p = re.compile('-?\d+x-?\d+') # The pattern for "INDEL"

        for row in _df_forMAP.astype(str).itertuples():

            """
            row[0] := index
            
            row[1] := rel_pos
            row[2] := gen_pos
            row[3] := type
            
            """

            # t = "INDEL" if bool(p.match(row[1])) else "AA" if _type == "AA" else "SNPS" # (2019. 03. 03.) SNP2HLA README.v2.txt수정하다가, Indel label 때문에 수정하기 전.

            if bool(p.match(row[1])):

                if _type == "AA":
                    t = "INDEL_AA"
                elif _type == "SNPS":
                    t = "INDEL_SNPS"

                main_label = '_'.join([t, _hla, *row[1:-1]])


            else:

                if _type == "AA":
                    t = "AA"
                elif _type == "SNPS":
                    t = "SNPS"

                main_label = '_'.join([t, _hla, *row[1:]])



            l_Label.append(main_label)

        l_Label = pd.Series(l_Label)

        # Output map file.
        df_MAP = pd.concat([sr_Chr, l_Label, sr_GD, sr_GenPos], axis=1)




    return df_MAP




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        IMGT2Sequence.py



    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
    parser.add_argument("-hg", help="\nHuman Genome version(18, 19 or 38)\n\n", choices=["18", "19", "38"],
                        metavar="hg", required=True)
    parser.add_argument("-o", help="\nOutput File Prefix\n\n", metavar="OUTPUT", required=True)

    parser.add_argument("-imgt", help="\nIMGT-HLA data version(ex. 370, 3300)\n\n", metavar="IMGT_Version",
                        required=True)

    parser.add_argument("--no-indel", help="\nExcluding indel in HLA sequence outputs.\n\n", action='store_true')
    parser.add_argument("--multiprocess", help="\nSetting off parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8], nargs='?', default=1, const=8)
    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')
    parser.add_argument("--imgt-dir", help="\nIn case User just want to specify the directory of IMGT data folder.\n\n")




    ##### < for Test > #####

    # (2018. 9. 17.)

    # Dictionary / hg18 / imgt370
    # args = parser.parse_args(["-hg", "18", "-o", "MAKEDICTIONARY_v2_hg18_imgt370/makedictionary.hg18.imgt370", "-imgt", "370"])

    # Dictionary / hg18 / imgt3320
    # args = parser.parse_args(["-hg", "18", "-o", "MAKEDICTIONARY_v2_hg18_imgt3320/makedictionary.hg18.imgt3320", "-imgt", "3320"])

    # Dictionary / hg19 / imgt3320
    # args = parser.parse_args(["-hg", "19", "-o", "MAKEDICTIONARY_v2_hg19_imgt3320/makedictionary.hg19.imgt3320", "-imgt", "3320"])




    ##### < for Publish > #####

    args = parser.parse_args()
    print(args)


    ##### < Main function Execution. > #####
    IMGT2Seq(_imgt=args.imgt, _hg=args.hg, _out=args.o, _no_Indel=args.no_indel, _MultiP=args.multiprocess,
             _save_intermediates=args.save_intermediates, _imgt_dir=args.imgt_dir)
