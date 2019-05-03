# -*- coding: utf-8 -*-

import os, sys, re
import logging
import numpy as np
import argparse, textwrap
import src.HLA_Analysis.HLA_Analysis_modules as hla_m



########## < Core Global Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


# (2018. 12. 19.)
class Study(object):

    def __init__(self, *args, **kwargs):

        """
        (Data)
        - HLA allele data (*.rhped, *.hped, *.chped)
        - Platform info. (AXIOM, etc.)
        - Normal SNP data (*.ped/*.map or *.bed/*.bim/*.fam)
        - Phenotype data (*.phe)
        - Phenotype name
        - Covariate data (*.cov or *.covar)
        - Covariate name
        - Condition data (*.condvars or just string)
        - Reference allele (*.ref or *.refallele)
        - Human Genome version
        - Output Prefix

        ex)
        mystudy = Study(platform=args.platform, rhped=args.rhped, hped=args.hped, chped=args.chped,
                        pheno=args.pheno, pheno_name=args.pheno_name, ...)

        * --input := the common prefix to point out next files
            (1) *.hped or *.chped
            (2) normal_SNPs
            (3) phenotype file(*.phe)
            (4) covariate file(*.covar)
            (5) condition (*.condvars)
            (6) reference allele (*.refallele)

        """


        ########## < [0] Assigning arguments > ##########

        # Dictionry Info.(It will be replaced by the instance of class "HLA_Dictionary" in IMGT2Sequences.
        self.dict_AA = kwargs["dict_AA"]
        self.dict_SNPS = kwargs["dict_SNPS"]
        self.iat = kwargs["iat"]
        self.d_MapTable = kwargs["d_MapTable"]

        # HLA ped
        self.platform = kwargs["platform"]
        self.rhped = kwargs["rhped"] # supposed to be 'None' or 'list'
        self.hped = None
        self.chped = None

        # Traditional genome SNP variants
        self.variants = None

        # Plink
        self.phenotype = None
        self.pheno_name = kwargs["pheno_name"]

        self.covar = None
        self.covar_name = kwargs["covar_name"]

        self.condition = None
        self.condition_list = None

        self.refallele = None

        # verison and etc.
        self.hg = kwargs["hg"]
        self.out = kwargs["out"]

        self.summary_string = \
        "< Summary info of HLA Study >\n" \
        "- Used IAT file : {}\n" \
        "- Human Genome : hg{}\n".format(self.iat, self.hg)


        # "--input" arguments to point out next 6 arguments.
        if kwargs["input"]:
            self.hped = kwargs["input"]+".hped" if os.path.exists(kwargs["input"]+".hped") else None
            self.chped = kwargs["input"]+".chped" if os.path.exists(kwargs["input"]+".chped") else None

            # `normal_SNPs` => just as prefix.
            self.variants = kwargs["input"] if os.path.exists(kwargs["input"] + ".bed") and os.path.exists(kwargs["input"] + ".bim") and os.path.exists(kwargs["input"] + ".fam") else None
            self.phenotype = kwargs["input"]+".phe" if os.path.exists(kwargs["input"]+".phe") else kwargs["input"]+".pheno" if os.path.exists(kwargs["input"]+".pheno") else None
            self.covar = kwargs["input"]+".covar" if os.path.exists(kwargs["input"]+".covar") else kwargs["input"]+".cov" if os.path.exists(kwargs["input"]+".cov") else None
            # `condition` => assuming it is a file with suffix('.condvars')
            self.condition_list = kwargs["input"] + ".condvars" if os.path.isfile(kwargs["input"] + ".condvars") else None
            self.refallele = kwargs["input"]+".refallele" if os.path.isfile(kwargs["input"]+".refallele") else kwargs["input"]+".ref" if os.path.exists(kwargs["input"]+".ref") else None


        # Override the results of "--input" arguments when the arguments for those 6 stuffs are given explicitly.
        if kwargs["hped"]:
            self.hped = kwargs["hped"]
        if kwargs["chped"]:
            self.chped = kwargs["chped"]

        if kwargs["variants"]:
            self.variants = kwargs["variants"]
        if kwargs["pheno"]:
            self.phenotype = kwargs["pheno"]
        if kwargs["covar"]:
            self.covar = kwargs["covar"]
        if kwargs["covar_name"]:
            self.covar_name = kwargs["covar_name"]
        if kwargs["condition"]:
            self.condition = kwargs["condition"]
        if kwargs["condition_list"]:
            self.condition_list = kwargs["condition_list"]
        if kwargs["refallele"]:
            self.refallele = kwargs["refallele"]


        """
        [0] Assigning Arguments
        [1] Preparing & Checking necessary input files
            [1-1] cleaned HLA ped file
            [1-2] Traditional genomewide SNPs
            [1-3] Phenotype info. and Phenotype name
            [1-4] Covariate file and Covariate names to use
            [1-5] Condition or Condition_list
        
        [2] b:MarkerGenerator
        [3] Association Test
            [3-1] Logistic Regression
            [3-2] Omnibus Test
            [3-3] Binary Test
        
        [4] Plotting
            [4-1] Manhattan Plot
            [4-2] Heatmap(for each HLA genes) 
        """



        ########## < [1] Preparing & Checking necessary input files > ##########

        ### [1-1] HLA ped (rhped, hped, or chped).

        self.__HLA_type__ = None

        if self.platform and self.rhped:
            from src.HLA2HPED.HLA2HPED import HATK_HLA2HPED
            from src.NomenCleaner.NomenCleaner import HATK_NomenCleaner
            self.__HLA_type__ = HATK_NomenCleaner(HATK_HLA2HPED(self.rhped, self.out, self.platform), 1, self.iat,
                                                  self.out, 4)

            self.summary_string = ''.join([self.summary_string, "- Platform : {}\n".format(self.platform)])

        elif self.chped:
            self.__HLA_type__ = self.chped

        elif self.hped:
            from src.NomenCleaner.NomenCleaner import HATK_NomenCleaner
            self.__HLA_type__ = HATK_NomenCleaner(self.hped, 1, self.iat, self.out, 4)

        else:
            # self.logger.error("Please check the arguments for HLA type data.\n")
            print(std_ERROR_MAIN_PROCESS_NAME + "No available HLA ped file! Please check (1) {}, (2) {}, (3) {} options again.".format("--input", "-chped", "-hped"))
            sys.exit()

        self.summary_string = ''.join([self.summary_string, "- (Cleaned) HLA ped: {}\n".format(self.__HLA_type__)])



        ### [1-2] Variants (normal SNP variants)

        if self.variants:

            self.summary_string = ''.join([self.summary_string,
                                           "- Variants(Normal SNPs) : {}\n".format(self.variants)])

        else:

            self.summary_string = ''.join([self.summary_string,
                                           "- Variants(Normal SNPs) : Not_given\n"])


        ### [1-3] Phenotye (*.phe and "--pheno-name")

        if self.phenotype:

            if not self.pheno_name:

                # only one phenotype columns is given.
                f = open(self.phenotype, 'r')
                columns = re.split(r'\s+', f.readline().rstrip('\n'))
                f.close()

                if columns[0] == "FID" and columns[1] == "IID":
                    if len(columns) == 3:
                        self.pheno_name = columns[2]
                        self.summary_string = ''.join([self.summary_string,
                                                       "- Phenotype file : {}\n"
                                                       "- Phenotype name : {}\n".format(self.phenotype, self.pheno_name)])

                    elif len(columns) > 3:
                        print(std_ERROR_MAIN_PROCESS_NAME + "The phenotype unspecified.\n"
                                                            "Please specify the phenotype to use.(\"--pheno-name\")\n")
                        sys.exit()

                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "The header line must start with \"FID\" and \"IID\".\n"
                                                        "Please check \"http://zzz.bwh.harvard.edu/plink/data.shtml#pheno\"\n")
                    sys.exit()

        else:

            c_phe = np.genfromtxt(self.__HLA_type__, usecols=[5], dtype=str)

            f_phe = (c_phe == '-9').all() # (2019. 01. 10.) Temporary hard-coding.
            # f_phe = (c_phe == '-9').all() or (c_phe == '0').all()

            if f_phe:
                print(std_ERROR_MAIN_PROCESS_NAME +
                      "The phenotype information must be given in the 6th column of \"*.chped\" file or as \"*.phe\" file. "
                      "Please check the phenotype information of your study.")
                sys.exit()
            else:
                self.summary_string = ''.join([self.summary_string,
                                               "- Phenotype name : {}\n"
                                               "- Phenotype file : {}\n".format("Not_given", "Included in the 6th column of \"*.chped\" file.")])


        ### [1-4] Covariate (*.cov(ar) and "--covar-name")

        if self.covar:
            # When *.cov(ar) file is given.

            if not self.covar_name:
                # Just handle this case as error.
                print(std_ERROR_MAIN_PROCESS_NAME + "The covariate file is given(\"{}\" : {}), but the name of covariate wasn't given(\"{}\").".format("--covar", self.covar, "--covar-name"))
                sys.exit()
            else:
                # Ready to use covariate file with covar_name.
                self.summary_string = ''.join([self.summary_string,
                                               "- Covariate name : {}\n"
                                               "- Covariate file : {}\n".format(self.covar_name, self.covar)])


        ### [1-5] Condition or Condition List

        if bool(self.condition) ^ bool(self.condition_list):

            # --condition-list
            if self.condition_list and isinstance(self.condition_list, str) and os.path.isfile(self.condition_list):

                self.summary_string = ''.join([self.summary_string,
                                               "- Condition List(as file) : {}\n".format(self.condition_list)])

            # --condition
            elif self.condition and isinstance(self.condition, str):

                self.summary_string = ''.join([self.summary_string,
                                               "- Condition : {}\n".format(self.condition)])

            else:
                print(std_WARNING_MAIN_PROCESS_NAME + "Something wrong with \"{}\" or \"{}\" arguments. Just skip them enforcibly.".format("--condition", "--condition-list"))
                self.summary_string = ''.join([self.summary_string,
                                               "- (Error) Condition or Condition_list was not applied. \n"])



        else:

            if self.condition and self.condition_list:
                print(std_WARNING_MAIN_PROCESS_NAME + "Both \"{}\" and \"{}\" arguments are given. \"{}\" argument is prioritized to be used enforcibly.".format("--condition", "--condition-list", "--condition"))





        ########## < [2] b:MarkerGenerator > ##########

        from src.b_MarkerGenerator.b_MarkerGenerator import b_MarkerGenerator

        # Generating b:Markers
        self.__b_MARKER_PANELS__ = b_MarkerGenerator(self.__HLA_type__, self.out, self.hg, self.dict_AA, self.dict_SNPS, _variants=self.variants)

        # byproducts of b:MarkerGenerators
        self.__b_Marker_AA__ = self.__b_MARKER_PANELS__+".AA.CODED"
        self.__b_Marker_HLA__ = self.__b_MARKER_PANELS__+".HLA"

        self.refallele = self.__b_MARKER_PANELS__+".refallele"
        self.refallele_AA = hla_m.MakeDefaultReferenceAllele(self.__b_Marker_AA__)
        self.refallele_HLA = hla_m.MakeDefaultReferenceAllele(self.__b_Marker_HLA__)


        # print("__b_Marker_Panel__ : {}\n"
        #       "__b_Marker_AA__ : {}\n"
        #       "__b_Marker_HLA__ : {}\n"
        #       "__refallele__ : {}".format(self.__b_MARKER_PANELS__, self.__b_Marker_AA__, self.__b_Marker_HLA__, self.refallele))

        self.summary_string = ''.join([self.summary_string,
                                       "- Generated Marker Panel : \n"
                                       "   1. __b_Marker__ : {}\n"
                                       "   2. __b_Marker_AA__ : {}\n"
                                       "   3. __b_Marker_HLA__ : {}\n".format(self.__b_MARKER_PANELS__, self.__b_Marker_AA__, self.__b_Marker_HLA__)])

        ########## < [3] Association Test > ##########

        ### [3-1] Logistic Regression
        self.assoc_logistic = hla_m.Logistic_Regression(_bfile=self.__b_MARKER_PANELS__, _out=self.out,
                                                        _covar=self.covar, _covar_names=self.covar_name,
                                                        _phe=self.phenotype, _phe_name=self.pheno_name,
                                                        _condition=self.condition, _condition_list=self.condition_list,
                                                        _ref_allele=self.refallele)

        self.assoc_logistic_AA = hla_m.Logistic_Regression(_bfile=self.__b_Marker_AA__, _out=self.__b_Marker_AA__,
                                                        _covar=self.covar, _covar_names=self.covar_name,
                                                        _phe=self.phenotype, _phe_name=self.pheno_name,
                                                        _condition=self.condition, _condition_list=self.condition_list,
                                                        _ref_allele=self.refallele)

        self.assoc_logistic_HLA = hla_m.Logistic_Regression(_bfile=self.__b_Marker_HLA__, _out=self.__b_Marker_HLA__,
                                                        _covar=self.covar, _covar_names=self.covar_name,
                                                        _phe=self.phenotype, _phe_name=self.pheno_name,
                                                        _condition=self.condition, _condition_list=self.condition_list,
                                                        _ref_allele=self.refallele)

        # print("- Association Results(logistic regression) : \n"
        #       "   1. Merged(whole) : {}\n"
        #       "   2. Amino Acids : {}\n"
        #       "   3. HLA alleles : {}".format(self.assoc_logistic, self.assoc_logistic_AA, self.assoc_logistic_HLA))

        self.summary_string = ''.join([self.summary_string,
                                       "- Association Results(logistic regression) : \n"
                                       "   1. Merged(whole) : {}\n"
                                       "   2. Amino Acids : {}\n"
                                       "   3. HLA alleles : {}\n".format(self.assoc_logistic, self.assoc_logistic_AA, self.assoc_logistic_HLA)])


        ### [3-2] Omnibus Test

        ### [3-3] Binary Test



        ########## < [4] Plotting > ##########

        from src.HLA_Analysis.Plotting.manhattanplot.manhattan import manhattan, HATK_manhattan
        from src.HLA_Analysis.Plotting.heatmap.heatmap import HEATMAP, HATK_HEATMAP


        ### [4-1] Manhattan
        # => Only `self.assoc_logistic`

        self.plot_Manhattan = HATK_manhattan([self.assoc_logistic], _plot_label="Cancer Test", _out=self.out, _hg=self.hg)

        self.summary_string = ''.join([self.summary_string,
                                       "- Manhattan plot : {}\n".format(self.plot_Manhattan)])


        ### [4-2] Heatmap
        # => Using `assoc_logistic_AA` and `self.assoc_logistic_HLA`
        self.Dict_HeatmapPlot_byHLA = {HLA_names[i]: -1 for i in range(0, len(HLA_names))} # dictionary initialization

        for i in range(0, len(HLA_names)):
            self.Dict_HeatmapPlot_byHLA[HLA_names[i]] = HATK_HEATMAP(HLA_names[i], self.out+".heatmap.{}".format(HLA_names[i]), self.d_MapTable[HLA_names[i]],
                                                                     self.assoc_logistic_AA, self.assoc_logistic_HLA)

        # print("- Heatmap plots by HLA : \n"
        #       "   A   : {A}\n"
        #       "   B   : {B}\n"
        #       "   C   : {C}\n"
        #       "   DPA1: {DPA1}\n"
        #       "   DPB1: {DPB1}\n"
        #       "   DQA1: {DQA1}\n"
        #       "   DQB1: {DQB1}\n"
        #       "   DRB1: {DRB1}\n".format(**self.Dict_HeatmapPlot_byHLA))

        self.summary_string = ''.join([self.summary_string,
                                       "- Heatmap plots by HLA : \n"
                                       "   A   : {A}\n"
                                       "   B   : {B}\n"
                                       "   C   : {C}\n"
                                       "   DPA1: {DPA1}\n"
                                       "   DPB1: {DPB1}\n"
                                       "   DQA1: {DQA1}\n"
                                       "   DQB1: {DQB1}\n"
                                       "   DRB1: {DRB1}\n".format(**self.Dict_HeatmapPlot_byHLA)])



    def __str__(self):
        return self.summary_string




def HATK_ASSOC1_Logistic_Regression(_input, _out, _covar=None, _covar_names=None, _phe=None, _phe_name=None,
                                    _condition=None, _condition_list=None, _ref_allele=None, _b_marker_panels=None):


    ##### < Core Variable > #####

    bFILE = None
    COVAR = None
    PHE = None
    REFERENCE_ALLELE = None


    # bed, bim, fam (Main panels for association tests)
    if bool(_b_marker_panels) and \
            (os.path.exists(_b_marker_panels+".bed") and os.path.exists(_b_marker_panels+".bim") and os.path.exists(_b_marker_panels+".fam")):

        bFILE = _b_marker_panels
        REFERENCE_ALLELE = (_b_marker_panels + ".refallele") if os.path.exists(_b_marker_panels + ".refallele") else hla_m.MakeDefaultReferenceAllele(_b_marker_panels)
        # _out = _b_marker_panels

    else:

        if bool(_input) and \
           (os.path.exists(_input+".bed") and os.path.exists(_input+".bim") and os.path.exists(_input+".fam")):

            bFILE = _input
            REFERENCE_ALLELE = (_input + ".refallele") if os.path.exists(_input + ".refallele") else hla_m.MakeDefaultReferenceAllele(_input)

        else:

            print(std_ERROR_MAIN_PROCESS_NAME + "The output panel of \"b:MarkerGenerator\" can't be found. Please check them again.\n")
            sys.exit()


    # _input (*.covar, covar-name, *.phe, etc...)
    if bool(_input):

        COVAR = (_input + ".covar") if os.path.exists(_input + ".covar") else (_input + ".cov") if os.path.exists(_input + ".cov") else None
        PHE = (_input + ".phe") if os.path.exists(_input + ".phe") else (_input + ".pheno") if os.path.exists(_input + ".pheno") else None
        # REFERENCE_ALLELE = (_input + ".refallele") if os.path.exists(_input + ".refallele") else hla_m.MakeDefaultReferenceAllele(_input)

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--input(-i)"))
        sys.exit()


    # _covar
    if bool(_covar):

        # When the argument "--covar" is given. (with or without "--covar-name" argument)
        if os.path.exists(_covar):
            COVAR = (_covar)
        elif os.path.exists(_covar):
            COVAR = (_covar)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The covariate file \"{0}\" given by the argument \"--covar\" doesn't exist. Please check that again.\n")
            sys.exit()

    # _pheno
    if bool(_phe):

        # When the argument "--pheno" is given. (with or without "--pheno-name" argument)
        if os.path.exists(_phe):
            PHE = (_phe)
        elif os.path.exists(_phe):
            PHE = (_phe)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The phenotype file \"{0}\" given by the argument \"--pheno\" doesn't exist. Please check that again.\n")
            sys.exit()


    # _covar_names
    if not bool(COVAR):

        print(std_MAIN_PROCESS_NAME + "No Covariate file has given.\n")

        if bool(_covar_names):
            print(std_WARNING_MAIN_PROCESS_NAME + "Covaraite file han't given but the argument \"--covar-names\" has given. "
                                                  "The argument \"--covar-names\" will be ignored.\n")
            _covar_names = None
    else:
        print(std_MAIN_PROCESS_NAME + "Given covariate file is \"{0}\"".format(COVAR))


    # _pheno_name
    if not bool(PHE):

        print(std_MAIN_PROCESS_NAME + "No Phenotype file has given.\n")

        if bool(_phe_name):
            print(std_WARNING_MAIN_PROCESS_NAME + "Phenotype file han't given but the argument \"--pheno-name\" has given. "
                                                  "The argument \"--pheno-name\" will be ignored.\n")
            _phe_name = None
    else:
        print(std_MAIN_PROCESS_NAME + "Given phenotype file is \"{0}\"".format(PHE))


    # _ref_allele
    if bool(_ref_allele):

        # When the argument "--reference-allele" is given.
        if os.path.exists(_ref_allele):
            REFERENCE_ALLELE = _ref_allele
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The reference allele file \"{0}\" given by the argument \"--reference-allele\" doesn't exist. Please check that again.\n")
            sys.exit()

    if not bool(REFERENCE_ALLELE):
        print(std_MAIN_PROCESS_NAME + "No Reference Allele file has given.\n")
    else:
        print(std_MAIN_PROCESS_NAME + "Given Reference Allele file is \"{0}\"".format(REFERENCE_ALLELE))


    # _codition and _condition_list
    if bool(_condition) and bool(_condition_list):
        print(std_ERROR_MAIN_PROCESS_NAME + "Either \"{0}\" or \"{1}\" should be given.\n".format("--condition", "--conditino-list"))
        sys.exit()

    elif not bool(_condition) and bool(_condition_list):

        # When _condition_list is given.
        if not os.path.isfile(_condition_list):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given file of the argument \"{0}\" doesn't exist. Please check it again.\n".format(_condition_list))
            sys.exit()


    # _out
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out(-o)"))
        sys.exit()
    else:
        # Intermediate path.
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)



    print(std_MAIN_PROCESS_NAME + "Conducting Logistic Regression for \"{0}\" (Plink v1.07).\n".format(bFILE))

    return hla_m.Logistic_Regression(_bfile=bFILE, _out=_out, _covar=COVAR, _covar_names=_covar_names, _phe=PHE,
                                     _phe_name=_phe_name, _condition=_condition, _condition_list=_condition_list,
                                     _ref_allele=REFERENCE_ALLELE)






def HATK_ASSOC2_Omnibus_Test(_input, _out, _phased, _phe, _phe_name, _covar, _covar_names, _condition="NA", _condition_list="NA"):

    """
    Association Test limited to AA markers.
    """

    ##### < Core Variable > #####

    PHE = None
    COVAR = None
    PHASED = None
    BGL = None


    # _input
    if bool(_input):
        COVAR = (_input + ".covar") if os.path.exists(_input + ".covar") else (_input + ".cov") if os.path.exists(_input + ".cov") else None
        PHE = (_input + ".phe") if os.path.exists(_input + ".phe") else (_input + ".pheno") if os.path.exists(_input + ".pheno") else None

        PHASED = (_input + ".aa") if os.path.exists(_input + ".aa") else None
        BGL = (_input + ".bgl.phased") if os.path.exists(_input + ".bgl.phased") else None

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--input(-i)"))
        sys.exit()



    # _covar
    if bool(_covar):

        # When the argument "--covar" is given. (with or without "--covar-name" argument)
        # Change `COVAR` from the one given by "_input" argument to another one given by "_covar" argument.
        if os.path.exists(_covar):
            COVAR = (_covar)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The covariate file \"{0}\" given by the argument \"--covar\" doesn't exist. Please check that again.\n")
            sys.exit()



    # _pheno
    if bool(_phe):

        # When the argument "--pheno" is given. (with or without "--pheno-name" argument)
        if os.path.exists(_phe):
            PHE = (_phe)
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "The phenotype file \"{0}\" given by the argument \"--pheno\" doesn't exist. Please check that again.\n")
            sys.exit()

        # If pheno file isn't given, then Omnibus Test can't be implemented and will be aborted.
        # Checking pheno file will be done below with checking pheno_name



    # _covar_names
    if not bool(COVAR):

        print(std_MAIN_PROCESS_NAME + "No Covariate file has given.\n")

        if bool(_covar_names):
            print(std_WARNING_MAIN_PROCESS_NAME + "Covaraite file han't given but the argument \"--covar-names\" has given. "
                                                  "The argument \"--covar-names\" will be ignored.\n")
            _covar_names = "NA"
    else:
        if not bool(_covar_names):
            # covar file given, but covar_name not given.
            _covar_names = "NA"
        else:
            # covar file and covar_name both given. => Just let it be passed to the function.
            print(std_MAIN_PROCESS_NAME + "Given column name(s) of covariate file is \"{0}\"".format(_covar_names))



    # _pheno_name
    if not bool(PHE):

        print(std_ERROR_MAIN_PROCESS_NAME + "Phenotype information must be given.\n")
        sys.exit()

        # 사실 여기서 ped파일 phenotype까지 확인해다가 ped파일에 phenotype정보 준비된거면 얘를 따로 파일로 빼다가 쓰고,
        # 아니면 여기서 abort시키는 코드를 추가해야함.

        # 우선 지금은 그냥 phenotype 파일 형태로 안주어지면 무조건 abort시키는 형태로.

        # (2018. 11. 2.) Retrieve this part later.
        # # On the othre hand, Pheno_name must be given.
        # print(std_MAIN_PROCESS_NAME + "No Phenotype file has given.\n")
        #
        # if bool(_phe_name):
        #     print(std_WARNING_MAIN_PROCESS_NAME + "Phenotype file han't given but the argument \"--pheno-name\" has given. "
        #                                           "The argument \"--pheno-name\" will be ignored.\n")
        #     _phe_name = "NA"

    else:
        if not bool(_phe_name):
            # phenotype file given, but phenotype_name not given.
            print(std_ERROR_MAIN_PROCESS_NAME + "Phenotype file has given({0}), but Phenotype name(\"--pheno-name\") hasn't given".format(_phe))
            sys.exit()
        else:
            # phenotype file and phenotype_name both given. => Just let it be passed to the function.
            print(std_MAIN_PROCESS_NAME + "Given column name(s) of phenotype file is \"{0}\"".format(_phe_name))



    # _codition and _condition_list
    if bool(_condition) and bool(_condition_list):

        print(std_ERROR_MAIN_PROCESS_NAME + "Either \"{0}\" or \"{1}\" should be given.\n".format("--condition", "--conditino-list"))
        sys.exit()

    elif not bool(_condition) and bool(_condition_list):

        """
        OmnibusTest script works only with condition_list not with condition as file.
        This code block is useless but some function to transform the contents of condition file to condition_list would be great.
        """

        # _condition_list(given as file)
        if not os.path.isfile(_condition_list):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given file of the argument \"{0}\" doesn't exist. Please check it again.\n".format(_condition_list))
            sys.exit()

    elif bool(_condition) and not bool(_condition_list):

        # _condition(given as csv string) => Just let it be passed to the function.
        pass



    # _out
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out(-o)"))
        sys.exit()
    else:
        # Intermediate path.
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)



    # _phased
    if bool(_phased):

        if os.path.exists(_phased):
            # Whatever *.aa file was given by the argument "--intput" or not, it will be overridden by "--phased" argument.
            PHASED = _phased
        else:

            # (1) *.aa파일이 어떤 이유에서든 존재하진 않으나 *.bgl.phased 파일로 *.aa파일은 만들 수 있는 경우.
            # (2) *.bgl.phased파일도 --input argument로 안주어지는 경우.


            if not bool(PHASED):

                # *.aa file can't be acquired by either "--input" or "--phased".

                if bool(BGL):
                    PHASED = hla_m.GetPhasedAlleles(_input, _out, )
                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "None of phased file(*.aa) or beagle phased file(*.bgl.phased) given. Please check them again.\n")
                    sys.exit()

            else:
                # *.aa file is prepared by "--input" argument. Just let it pass to the function.
                pass

    """
    ### (2018. 10. 11.) 여기 잠깐 보류해놓음. 나중에 각 implementation wrapper 각각 도입하고 여기좀 다시 돌아올것.
    ### (2018. 11. 1.) 여기 또 아직 해결 안 된 문제가 뭐였냐면 covar, pheno파일의 특정 컬럼을 covar_name, pheno_name으로 지시 못하는거.
    
    # (2018. 11. 2.)
    아직도 해결 안된 문제 : 
    (1) phenotype file & pheno-name형태로 주어지지 않고 그냥 ped에 phenotype 정보가 주어지는 경우,
    (2) Condition이 파일로 주어진 경우, csv string으로 바꾸는 기능.
    """

    return hla_m.Omnibus_Test(_input, _out, PHASED, PHE, _phe_name, COVAR, _covar_names, _condition)





def HATK_META_ANALYSIS(_out, _rassoc):


    ##### < Argument Checking > #####

    # _out
    if not bool(_out):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--out\" wasn't given. Please check it again.\n")
        sys.exit()
    else:
        # Intermediate path.
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)


    # _rassoc
    if not (bool(_rassoc) or isinstance(_rassoc, list)):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--results-assoc\" has something wrong. Please check it again.\n")
        sys.exit()


    # (2018. 8. 6.) As far as input association result files have extension defined by plink(ex. .assoc, .fisher, .assoc.logstic, etc.), Meta-analysis can be conducted.
    # If the result of Omnibus test also could be used in meta-analysis maybe...


    return hla_m.Meta_Analysis(_out, _rassoc)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###############################################################################       
        
        HLA_Analysis.py
        
        The modules prepared in "HLA_Analysis_modules.py" will be used in this 
        script. Checking and Preprocessing arguments will be conducted too.
        
        made by Wanson Choi.
        
    ###############################################################################
                                             '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-i", help="\nInput Data file.\n\n", required=True)
    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)


    # for Publish
    args = parser.parse_args()
    # print(args)

    # for Testing


    # HLA_Analysis()