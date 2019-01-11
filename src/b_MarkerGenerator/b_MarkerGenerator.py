# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
# from platform import platform


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


# def HATK_b_MarkerGenerator(_CHPED, _OUT, _hg, _dictionary_AA, _dictionary_SNPS, _plain_SNP_DATA=None):
def HATK_b_MarkerGenerator(_b_MarkerGenerator):

    def wrapper_function(*args, **kwargs):

        if __name__ == "__main__":

            # ##### Additional Argument processing
            #
            # t_dict_AA = ""
            # t_dict_SNPS = ""
            #
            # if (args.dict_AA != "Not_given" and args.dict_SNPS != "Not_given"):
            #
            #     # When all HLA DICTIONARY information is given properly,
            #
            #     t_dict_AA = args.dict_AA
            #     t_dict_SNPS = args.dict_SNPS
            #
            #
            # elif (args.dict_AA == "Not_given" and args.dict_SNPS == "Not_given"):
            #
            #     # No values are given to HLA DICTIONARY related options.
            #
            #     # Abort
            #     print("\n[Error]: None of HLA DICTIONARY files are given. Please check them all again.")
            #     print('{"-dict-AA", "-dict-SNPS"}\n')
            #     sys.exit()
            #
            #
            #
            # else:
            #     # Abort
            #     print("\n[Error]: Not all of HLA DICTIONARY files are given. Please check them all again.")
            #     print('{"-dict-AA", "-dict-SNPS"}\n')
            #     sys.exit()
            pass

        else:

            print(args)
            print(kwargs)

            _CHPED = args[0]
            _OUT = args[1]
            _HG = args[2]
            _dict_AA = args[3]
            _dict_SNPS = args[4]
            _variants = kwargs["_variants"]


            if not bool(_CHPED):
                print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hped"))
                sys.exit()
            else:
                # Checking whether it went through "NomenCleaner.py".
                if not _CHPED.endswith(".chped"):
                    print(std_ERROR_MAIN_PROCESS_NAME + 'Given ped file should be processed by "NomenCleaner.py"\n')
                    sys.exit()


            if not bool(_HG):
                print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("-hg"))
                sys.exit()


            if not (bool(_dict_AA) and bool(_dict_SNPS)):
                print(std_ERROR_MAIN_PROCESS_NAME + 'One of "-dict-AA" or "-dict-SNPS" options(or Both) are not given. Please check them again.\n')
                sys.exit()

        return _b_MarkerGenerator(*args, **kwargs)

    return wrapper_function



@HATK_b_MarkerGenerator
def b_MarkerGenerator(_CHPED, _OUT, _hg, _dictionary_AA, _dictionary_SNPS, _variants=None,
                      _p_src=None, _p_dependency=None):


    ########## < Core Variables > ##########

    ### Major Path Variables

    # [1] src (with "_p_src")
    p_src_b_MarkerGenerator = "src/b_MarkerGenerator" if not bool(_p_src) else _p_src
    p_dependency = "dependency" if not bool(_p_dependency) else _p_dependency


    # [2] dependency (with "_p_dependency")
    _p_plink = os.path.join(p_dependency, "plink")
    # _p_plink = os.path.join(p_dependency, "plink_mac" if not bool(re.search(pattern="Linux", string=platform())) else "plink_linux")
    _p_beagle = os.path.join(p_dependency, "beagle.jar")
    _p_linkage2beagle = os.path.join(p_dependency, "linkage2beagle.jar")


    ### Dictionary Files

    _dictionary_AA_seq = _dictionary_AA + ".txt" # From now on, official extension of HLA sequence information dictionary is ".txt". (2018. 9. 25.)
    _dictionary_AA_map = _dictionary_AA + ".map"

    _dictionary_SNPS_seq = _dictionary_SNPS + ".txt"
    _dictionary_SNPS_map = _dictionary_SNPS + ".map"


    ### Intermediate path.
    OUTPUT = _OUT if not _OUT.endswith('/') else _OUT.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)


    ### Flag for plain SNP markers.
    f_plain_SNP = bool(_variants)


    ########## < Checking Dependencies > ##########

    ### Other Software.

    if not os.path.exists(_p_plink):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(os.path.dirname(_p_plink)))
        sys.exit()

    if not os.path.exists(_p_beagle):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'Beagle 3' (http://faculty.washington.edu/browning/beagle/beagle.html#download) in '{0}'\n".format(os.path.dirname(_p_beagle)))
        sys.exit()

    if not os.path.exists(_p_linkage2beagle):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'linkage2beagle.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) in '{0}'\n".format(os.path.dirname(_p_linkage2beagle)))
        sys.exit()


    ### Dictionary Information for HLA sequence

    if not os.path.exists(_dictionary_AA_map):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_AA.map' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_AA_map)))
        sys.exit()

    if not os.path.exists(_dictionary_AA_seq):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_AA.txt' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_AA_seq)))
        sys.exit()

    if not os.path.exists(_dictionary_SNPS_map):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_SNPS.map' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_SNPS_map)))
        sys.exit()

    if not os.path.exists(_dictionary_SNPS_seq):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_SNPS.txt' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_SNPS_seq)))
        sys.exit()


    ### Source Code Scripts

    # New version with Python.

    if not os.path.exists(os.path.join(p_src_b_MarkerGenerator, "HLAtoSequences.py")):
        print(std_MAIN_PROCESS_NAME + "Error. 'HLAtoSequences.py' not found in '{0}'".format(p_src_b_MarkerGenerator))
        sys.exit()
    else:
        from src.b_MarkerGenerator.HLAtoSequences import HLAtoSequences

    if not os.path.exists(os.path.join(p_src_b_MarkerGenerator, "encodeVariants.py")):
        print(std_MAIN_PROCESS_NAME + "Error. 'encodeVariants.py' not found in '{0}'".format(p_src_b_MarkerGenerator))
        sys.exit()
    else:
        from src.b_MarkerGenerator.encodeVariants import encodeVariants

    if not os.path.exists(os.path.join(p_src_b_MarkerGenerator, "encodeHLA.py")):
        print(std_MAIN_PROCESS_NAME + "Error. 'encodeHLA.py' not found in '{0}'".format(p_src_b_MarkerGenerator))
        sys.exit()
    else:
        from src.b_MarkerGenerator.encodeHLA import encodeHLA




    ########## < Core Variables 2 > ##########

    # Input 1 : HLA type data
    HLA_DATA = _CHPED

    # Input 2 : Plain SNP data
    if f_plain_SNP:

        SNP_DATA = _variants
        SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, os.path.basename(_variants))


    # Output prefix
    AA_CODED = OUTPUT + '.AA.CODED'
    HLA_CODED = OUTPUT + ".HLA"
    SNPS_CODED = OUTPUT + '.SNPS.CODED'

    plink = ' '.join([_p_plink, "--noweb", "--silent"])
    beagle = ' '.join(["java", "-Xmx4000m", "-jar", _p_beagle])
    linkage2beagle = ' '.join(["java", "-Xmx2000m", "-jar", _p_linkage2beagle])




    ########## <Flags for Code Block> ##########

    ENCODE_AA = 1
    ENCODE_HLA = 1
    ENCODE_SNPS = 1

    EXTRACT_FOUNDERS = 1
    MERGE = 1
    QC = 1

    PREPARE = 0
    PHASE = 0
    CLEANUP = 0

    # (2019. 01. 10.) Last three code block won't be implemented



    ########## <Making Reference Panel> ##########

    print(std_MAIN_PROCESS_NAME + "Making Reference Panel for \"{0}\"\n".format(OUTPUT))
    index = 1

    if ENCODE_AA:

        '''
        echo "[$i] Generating amino acid sequences from HLA types.";  @ i++
        ./HLAtoSequences.pl $HLA_DATA HLA_DICTIONARY_AA.txt AA > $OUTPUT.AA.ped
        cp HLA_DICTIONARY_AA_hg19.map $OUTPUT.AA.map # hg19
        # cp HLA_DICTIONARY_AA.map $OUTPUT.AA.map

        echo "[$i] Encoding amino acids positions." ;  @ i++
        ./encodeVariants.pl $OUTPUT.AA.ped $OUTPUT.AA.map $OUTPUT.AA.CODED

        plink --file $OUTPUT.AA.CODED --missing-genotype 0 --make-bed --out $OUTPUT.AA.TMP
        awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.AA.TMP.bim | grep -v INS | cut -f2 > to_remove
        plink --bfile $OUTPUT.AA.TMP --exclude to_remove --make-bed --out $OUTPUT.AA.CODED

        # rm $OUTPUT.AA.TMP*; rm to_remove
        # rm $OUTPUT.AA.???
        '''

        print("[{}] Generating Amino acid(AA)sequences from HLA types.".format(index))

        HLAtoSequences(HLA_DATA, _dictionary_AA_seq, "AA", _out=OUTPUT)

        os.system(' '.join(["cp", _dictionary_AA_map, OUTPUT + '.AA.map']))

        index += 1


        print("[{}] Encoding Amino acids positions.".format(index))

        encodeVariants(OUTPUT + '.AA.ped', OUTPUT + '.AA.map', OUTPUT + '.AA.CODED')  # previously "enCODED".

        # command for checking output from encodeVariant.py(.pl)
        command = ' '.join(
            [plink, "--file", OUTPUT + '.AA.CODED', "--missing-genotype 0", "--make-bed", "--out", OUTPUT + '.AA.TMP'])
        # print(command)
        os.system(command)

        # MakeReference
        # command = ' '.join(
        #     ["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.AA.TMP.bim', "|",
        #      "grep -v {0}".format("INDEL"), "|", "cut -f2", ">",
        #      os.path.join(INTERMEDIATE_PATH, "to_remove")])

        # HLA2MARKER
        command = ' '.join(
            ["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.AA.TMP.bim', "|",
             "cut -f2", ">",
             os.path.join(INTERMEDIATE_PATH, "to_remove")])

        """
        In previous framework originally created by Sherman Jia, Only insertions were dealt with as a marker "INS".
        In the new version of Framework, marker label is "INDEL".
        """

        # print(command)
        os.system(command)

        command = ' '.join(
            [plink, "--bfile", OUTPUT + '.AA.TMP',
             "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"),
             "--make-bed",
             "--out", OUTPUT + '.AA.CODED'])
        # print(command)
        os.system(command)

        rm_tlist = (OUTPUT + '.AA.TMP*', os.path.join(INTERMEDIATE_PATH, "to_remove"), OUTPUT + '.AA.???')

        for i in rm_tlist:
            # print(i)
            os.system("rm " + i)

        index += 1



    if ENCODE_HLA:

        print("[{}] Encoding HLA alleles.".format(index))

        encodeHLA(HLA_DATA, OUTPUT + ".HLA", _hg)

        command = ' '.join([plink, "--file", OUTPUT + '.HLA', "--make-bed", "--out", OUTPUT + '.HLA'])
        # print(command)
        os.system(command)

        index += 1



    if ENCODE_SNPS:

        print("[{}] Generating DNA(SNPS) sequences from HLA types.".format(index))

        HLAtoSequences(HLA_DATA, _dictionary_SNPS_seq, "SNPS", OUTPUT)

        command = ' '.join(["cp", _dictionary_SNPS_map, OUTPUT + '.SNPS.map'])
        # print(command)
        os.system(command)

        index += 1


        print("[{}] Encoding SNP positions.".format(index))

        encodeVariants(OUTPUT + '.SNPS.ped', OUTPUT + '.SNPS.map', OUTPUT + '.SNPS.CODED')

        command = ' '.join([plink, "--file", OUTPUT + '.SNPS.CODED', "--missing-genotype 0", "--make-bed", "--out",
                            OUTPUT + '.SNPS.TMP'])
        # print(command)
        os.system(command)


        # # MakeReference
        # command = ' '.join(
        #     ["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.SNPS.TMP.bim', "|",
        #      "grep -v {0}".format("INDEL"), "|", "cut -f2", ">",
        #      os.path.join(INTERMEDIATE_PATH, "to_remove")])

        # HLA2MARKER
        command = ' '.join(
            ["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.SNPS.TMP.bim', "|", "cut -f2", ">",
             os.path.join(INTERMEDIATE_PATH, "to_remove")])
        # print(command)
        os.system(command)

        command = ' '.join(
            [plink, "--bfile", OUTPUT + '.SNPS.TMP',
             "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"),
             "--make-bed",
             "--out", OUTPUT + '.SNPS.CODED'])
        # print(command)
        os.system(command)

        rm_tlist = (OUTPUT + '.SNPS.TMP*', os.path.join(INTERMEDIATE_PATH, 'to_remove'), OUTPUT + '.SNPS.???')

        for i in rm_tlist:
            # print(i)
            os.system("rm " + i)

        index += 1


    # Above 3 code blocks are implmented whether `_variants`(Normal genomewide SNPs) is given or not.
    # Next parts are implented sligtly differently depending on `_variants`.


    if _variants:

        ### `_variants` is given.
        # Almost same as original "MakeReference".

        if EXTRACT_FOUNDERS:

            print("[{}] Extracting founders.".format(index))

            """
            if ($EXTRACT_FOUNDERS) then
                echo "[$i] Extracting founders."; @ i++
                plink --bfile $SNP_DATA --filter-founders --mind 0.3 --alleleACGT --make-bed --out $SNP_DATA.FOUNDERS
            
                # Initial QC on Reference SNP panel
                plink --bfile $SNP_DATA.FOUNDERS --hardy        --out $SNP_DATA.FOUNDERS.hardy  # 진짜 92명에 대해 position별로 HWE test한 결과
                plink --bfile $SNP_DATA.FOUNDERS --freq         --out $SNP_DATA.FOUNDERS.freq   # 실제 --freq 옵션이 allele frequency계산해주는 옵션임. 
                plink --bfile $SNP_DATA.FOUNDERS --missing      --out $SNP_DATA.FOUNDERS.missing
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe      | awk ' $9 < 0.000001 { print $2 }' | sort -u > remove.snps.hardy 
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.freq.frq       | awk ' $5 < 0.01 { print $2 } '             > remove.snps.freq
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.missing.lmiss  | awk ' $5 > 0.05 { print $2 } '              > remove.snps.missing
                cat remove.snps.*                                            | sort -u                                     > all.remove.snps
            
                plink --bfile $SNP_DATA.FOUNDERS --allow-no-sex --exclude all.remove.snps --make-bed --out $SNP_DATA.FOUNDERS.QC
            
                # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file
            
                plink --bfile $OUTPUT.HLA --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.HLA.FOUNDERS
                plink --bfile $OUTPUT.SNPS.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.SNPS.FOUNDERS
                plink --bfile $OUTPUT.AA.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.AA.FOUNDERS
            
                rm remove.snps.*
            endif
            """

            command = ' '.join([plink, "--bfile", SNP_DATA, "--filter-founders", "--mind 0.3", "--alleleACGT", "--make-bed", "--out", SNP_DATA2+'.FOUNDERS'])
            # print(command)
            os.system(command)


            # Initial QC on Reference SNP panel
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--hardy", "--out", SNP_DATA2+'.FOUNDERS.hardy'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--freq", "--out", SNP_DATA2+'.FOUNDERS.freq'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--missing", "--out", SNP_DATA2+'.FOUNDERS.missing'])
            # print(command)
            os.system(command)

            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.hardy.hwe', "|", "awk", "' $9 < 0.000001 { print $2 }'", "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.hardy")])
            # print(command)
            os.system(command)
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.freq.frq', "|", "awk", "' $5 < 0.01 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.freq")])
            # print(command)
            os.system(command)
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.missing.lmiss', "|", "awk", "' $5 > 0.05 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.missing")])
            # print(command)
            os.system(command)
            command = ' '.join(["cat", os.path.join(INTERMEDIATE_PATH, "remove.snps.*"), "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "all.remove.snps")])
            # print(command)
            os.system(command)


            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--allow-no-sex", "--exclude", os.path.join(INTERMEDIATE_PATH, "all.remove.snps"), "--make-bed", "--out", SNP_DATA2+'.FOUNDERS.QC'])
            # print(command)
            os.system(command)

            # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

            command = ' '.join([plink, "--bfile", OUTPUT+'.HLA', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.HLA.FOUNDERS'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.SNPS.FOUNDERS'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", OUTPUT+'.AA.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.AA.FOUNDERS'])
            # print(command)
            os.system(command)

            command = ' '.join(["rm", os.path.join(INTERMEDIATE_PATH, "remove.snps.*")])
            # print(command)
            os.system(command)


            index += 1


        if MERGE:

            print("[{}] Merging SNP, HLA, and amino acid datasets.".format(index))

            """
            echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
            echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
            echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
            echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
            plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
            rm $OUTPUT.HLA.???
            rm $OUTPUT.AA.CODED.???
            rm $OUTPUT.SNPS.CODED.???
            rm merge_list
    
            """

            TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")

            # plain_SNP_Data
            # original MakeReference.

            command = ' '.join(
                ["echo", OUTPUT + '.HLA.FOUNDERS.bed', OUTPUT + '.HLA.FOUNDERS.bim', OUTPUT + '.HLA.FOUNDERS.fam', ">",
                 TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(
                ["echo", OUTPUT + '.AA.FOUNDERS.bed', OUTPUT + '.AA.FOUNDERS.bim', OUTPUT + '.AA.FOUNDERS.fam', ">>",
                 TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(
                ["echo", OUTPUT + '.SNPS.FOUNDERS.bed', OUTPUT + '.SNPS.FOUNDERS.bim', OUTPUT + '.SNPS.FOUNDERS.fam',
                 ">>", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(
                [plink, "--bfile", SNP_DATA2 + '.FOUNDERS.QC', "--merge-list", TMP_merged_list, "--make-bed", "--out",
                 OUTPUT + '.MERGED.FOUNDERS'])
            # print(command)
            os.system(command)

            rm_tlist = (OUTPUT + '.HLA.???', OUTPUT + '.AA.CODED.???', OUTPUT + '.SNPS.CODED.???', TMP_merged_list)

            # (2019. 01. 11.) Decided to remove them after association test is done.
            # for i in rm_tlist:
            #     # print(i)
            #     os.system("rm " + i)

            index += 1


        if QC:

            print("[{}] Performing quality control.".format(index))

            """
            plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
            awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
            awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order
    
            # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
            plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT
    
            # Calculate allele frequencies
            plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
            rm $SNP_DATA.FOUNDERS.*
            rm $OUTPUT.MERGED.FOUNDERS.*
            rm $OUTPUT.*.FOUNDERS.???
            rm allele.order
            rm all.remove.snps
    
            """

            TMP_allele_order = os.path.join(INTERMEDIATE_PATH, "allele.order")
            TMP_all_remove_snps = os.path.join(INTERMEDIATE_PATH, "all.remove.snps")

            command = ' '.join(
                [plink, "--bfile", OUTPUT + '.MERGED.FOUNDERS', "--freq", "--out", OUTPUT + '.MERGED.FOUNDERS.FRQ'])
            # print(command)
            os.system(command)

            command = ' '.join(["awk", "'{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}'",
                                OUTPUT + '.MERGED.FOUNDERS.FRQ.frq', ">", TMP_all_remove_snps])
            # print(command)
            os.system(command)

            # command = ' '.join(["awk", '\'{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}\'', OUTPUT+'.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
            command = ' '.join(
                ["awk", '\'{if (NR > 1){if (($3 == "a" && $4 == "p") || ($4 == "a" && $3 == "p")){print $2 "\tp"}}}\'',
                 OUTPUT + '.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
            # print(command)
            os.system(command)

            # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
            command = ' '.join(
                [plink, "--bfile", OUTPUT + '.MERGED.FOUNDERS', "--reference-allele", TMP_allele_order, "--exclude",
                 TMP_all_remove_snps, "--geno 0.5", "--make-bed", "--out", OUTPUT])
            # print(command)
            os.system(command)

            # Calculate allele frequencies
            command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--freq", "--out", OUTPUT + '.FRQ'])
            # print(command)
            os.system(command)

            rm_tlist = (
            SNP_DATA2 + '.FOUNDERS.*', OUTPUT + '.MERGED.FOUNDERS.*', OUTPUT + '.*.FOUNDERS.???', TMP_allele_order,
            TMP_all_remove_snps)

            for i in rm_tlist:
                # print(i)
                os.system("rm " + i)

            index += 1



        if PREPARE:

            print("[{}] Preparing files for Beagle.".format(index))

            """
            [Source from Buhm Han.]
    
            awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
            plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
            awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
            cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
    
            echo "[$i] Converting to beagle format.";  @ i++
            linkage2beagle pedigree=$OUTPUT.nopheno.ped data=$OUTPUT.dat beagle=$OUTPUT.bgl standard=true > $OUTPUT.bgl.log
    
    
            [Source from Yang.]
    
            awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
            plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
            plink --bfile $OUTPUT --recode --transpose --out $OUTPUT
            # awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
            # cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
    
            echo "[$i] Converting to beagle format.";  @ i++
            beagle2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnout $OUTPUT.vcf
    
            I will make this code block based on source given by Yang. for now.
    
            """


            command = ' '.join(
                ["awk", '\'{print $2 " " $4 " " $5 " " $6}\'', OUTPUT + '.bim', ">", OUTPUT + '.markers'])
            print(command)
            os.system(command)

            command = ' '.join(
                [plink, "--bfile", OUTPUT, "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT])
            print(command)
            os.system(command)

            command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT + '.map', ">", OUTPUT + '.dat'])
            print(command)
            os.system(command)

            command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT + '.ped', ">", OUTPUT + '.nopheno.ped'])
            print(command)
            os.system(command)

            index += 1

            print("[{}] Converting to beagle format.".format(index))

            command = ' '.join([linkage2beagle, "pedigree=" + OUTPUT + '.nopheno.ped', "data=" + OUTPUT + '.dat',
                                "beagle=" + OUTPUT + '.bgl', "standard=true", ">", OUTPUT + '.bgl.log'])
            print(command)
            os.system(command)

            index += 1


        if PHASE:

            print("[{}] Phasing reference using Beagle (see progress in $OUTPUT.bgl.log).".format(index))

            # Put this part postponed. (2017.11.29. by B. Han.)
            # Anyway introduced phasing by beagle. (2018. 7. 16.)

            '''
            beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log
    
            '''

            command = ' '.join([beagle, "unphased=" + OUTPUT + '.bgl',
                                "nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000",
                                "log=" + OUTPUT + '.phasing', ">>", OUTPUT + '.bgl.log'])
            print(command)
            os.system(command)


            index += 1



        if CLEANUP:

            print("[{}] Removing unnecessary files.".format(index))

            '''
            rm $OUTPUT.nopheno.ped
            rm $OUTPUT.bgl.gprobs
            rm $OUTPUT.bgl.r2
            rm $OUTPUT.bgl
            rm $OUTPUT.ped
            rm $OUTPUT.map
            rm $OUTPUT.dat
            rm $OUTPUT.phasing.log
            '''

            rm_tlist = ('.nopheno.ped', '.bgl.gprobs', '.bgl.r2', '.bgl', '.ped', '.map', '.dat')

            for i in rm_tlist:
                print("rm " + OUTPUT + i)
                os.system("rm " + OUTPUT + i)

            index += 1


        print("[{}] Making reference panel for HLA-AA,SNPS,HLA and Normal variants(SNPs) is Done!".format(index))

        __return__ = OUTPUT



    else:

        ### No `_variants` is given.
        # Only for HLA markers.
        # "EXTRACT_FOUNDERS" code block is not included.

        print("(Just HLA type data.)")


        if MERGE:

            print("[{}] Merging SNP, HLA, and amino acid datasets.".format(index))

            """
            echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
            echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
            echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
            echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
            plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
            rm $OUTPUT.HLA.???
            rm $OUTPUT.AA.CODED.???
            rm $OUTPUT.SNPS.CODED.???
            rm merge_list
    
            """

            TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")


            command = ' '.join(
                ["echo", AA_CODED + '.bed', AA_CODED + '.bim', AA_CODED + '.fam', ">", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(
                ["echo", HLA_CODED + '.bed', HLA_CODED + '.bim', HLA_CODED + '.fam', ">>", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(
                ["echo", SNPS_CODED + '.bed', SNPS_CODED + '.bim', SNPS_CODED + '.fam', ">>", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join([plink, "--merge-list", TMP_merged_list, "--make-bed", "--out", OUTPUT + '.MERGED.TMP'])
            # print(command)
            os.system(command)


            # (2019. 01. 11.) Decided to remove them after association test is done.
            # rm_list = [AA_CODED, HLA_CODED, SNPS_CODED]
            # file_extensions = [".bed", ".bim", ".fam", ".log", ".map", ".ped"]
            #
            # for i in rm_list:
            #     for j in file_extensions:
            #         # print("rm " + i)
            #         os.system("rm " + i + j)

            index += 1


        if QC:

            """
            plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
            awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
            awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order
        
            # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
            plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT
        
            # Calculate allele frequencies
            plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
            rm $SNP_DATA.FOUNDERS.*
            rm $OUTPUT.MERGED.FOUNDERS.*
            rm $OUTPUT.*.FOUNDERS.???
            rm allele.order
            rm all.remove.snps
            
            """

            print("[{}] (only HLA) Generating \"*.reference\" file. (\"p\", \"a\")".format(index))


            TMP_allele_order = OUTPUT + ".refallele"

            command = ' '.join(
                ["awk", '\'{if (NR > 1){if (($5 == "a" && $6 == "p") || ($6 == "a" && $5 == "p")){print $2 "\tp"}}}\'',
                 OUTPUT + '.MERGED.TMP.bim', ">", TMP_allele_order])
            # print(command)
            os.system(command)

            command = ' '.join(
                [plink, "--make-bed", "--bfile", OUTPUT + '.MERGED.TMP', "--reference-allele", TMP_allele_order,
                 "--out", OUTPUT]) # (2019. 01. 10.) Final output as just output prefix(`OUTPUT`)
            # print(command)
            os.system(command)

            rm_tlist = [OUTPUT + '.MERGED.TMP.*', TMP_merged_list]

            for i in rm_tlist:
                # print(i)
                os.system("rm " + i)

            index += 1


        if PREPARE:

            print("[{}] Preparing files for Beagle.".format(index))

            """
            [Source from Buhm Han.]
            
            awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
            plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
            awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
            cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
        
            echo "[$i] Converting to beagle format.";  @ i++
            linkage2beagle pedigree=$OUTPUT.nopheno.ped data=$OUTPUT.dat beagle=$OUTPUT.bgl standard=true > $OUTPUT.bgl.log
            
            
            [Source from Yang.]
            
            awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
            plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
            plink --bfile $OUTPUT --recode --transpose --out $OUTPUT
            # awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
            # cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
        
            echo "[$i] Converting to beagle format.";  @ i++
            beagle2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnout $OUTPUT.vcf
            
            I will make this code block based on source given by Yang. for now.
            
            """

            # *.markers
            command = ' '.join(
                ["awk", '\'{print $2 " " $4 " " $5 " " $6}\'', OUTPUT + '.bim', ">", OUTPUT + '.markers'])
            print(command)
            os.system(command)

            command = ' '.join(
                [plink, "--bfile", OUTPUT, "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT])
            print(command)
            os.system(command)

            # *.dat
            command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT + '.map', ">", OUTPUT + '.dat'])
            print(command)
            os.system(command)

            # *.nopheno.ped
            command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT + '.ped', ">", OUTPUT + '.nopheno.ped'])
            print(command)
            os.system(command)

            print("[{}] Converting to beagle format.".format(index))

            command = ' '.join(
                [linkage2beagle, "pedigree=" + OUTPUT + '.nopheno.ped', "data=" + OUTPUT + '.dat',
                 "beagle=" + OUTPUT + '.bgl', "standard=true", ">", OUTPUT + '.bgl.log'])
            print(command)
            os.system(command)


            index += 1



        if PHASE:

            print("[{}] Phasing reference using Beagle (see progress in $OUTPUT.bgl.log).".format(index))

            # Put this part postponed. (2017.11.29. by B. Han.)
            # Anyway introduced phasing by beagle. (2018. 7. 16.)

            '''
            beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log
    
            '''

            command = ' '.join(
                [beagle, "unphased=" + OUTPUT + '.bgl',
                 "nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000",
                 "log=" + OUTPUT + '.phasing', ">>", OUTPUT + '.bgl.log'])
            print(command)
            os.system(command)


            index += 1


        if CLEANUP:

            print("[{}] Removing unnecessary files.".format(index))
            '''
            rm $OUTPUT.nopheno.ped
            rm $OUTPUT.bgl.gprobs
            rm $OUTPUT.bgl.r2
            rm $OUTPUT.bgl
            rm $OUTPUT.ped
            rm $OUTPUT.map
            rm $OUTPUT.dat
            rm $OUTPUT.phasing.log
            '''

            rm_tlist = ('.nopheno.ped', '.bgl.gprobs', '.bgl.r2', '.bgl', '.ped', '.map', '.dat')

            for i in rm_tlist:
                print("rm " + OUTPUT + ".MERGED" + i) # (2019. 01. 10.) ".MERGED" should be removed.
                os.system("rm " + OUTPUT + ".MERGED" + i)

            index += 1

        print("[{}] Making reference panel for HLA-AA,SNPS,HLA is Done!".format(index))

        # __return__ =


    return _OUT


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        b_MarkerGenerator.py
        
        Generating markers based on HLA sequence information dictionary(generated by "MakeDictionary").

        (ex.)
        : python3 b_MarkerGenerator.py  
            -chped ./data/b_MarkerGenerator.py/HAPMAP_CEU_HLA.4field.ped 
            -hg 18 
            -o ./Trial_HAPMAP_CEU
            -dict-AA ./data/b_MarkerGenerator.py/HLA_DICTIONARY_AA.hg18.imgt370
            -dict-SNPS ./data/b_MarkerGenerator.py/HLA_DICTIONARY_SNPS.hg18.imgt370 

        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-variants", help="\nInput variants data file(.bed/.bim/.fam)\n\n")
    parser.add_argument("-chped", help="\nHLA Type Data(.chped)\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
    parser.add_argument("-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("-dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n")
    parser.add_argument("-dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")



    ##### <for Test> #####

    # 2019. 01. 10
    # args = parser.parse_args(["-chped", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU_HLA.imgt370.4field.chped",
    #                           "-variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
    #                           "-hg", "18",
    #                           "-dict-AA", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_AA.hg18.imgt370",
    #                           "-dict-SNPS", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_SNPS.hg18.imgt370"
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)


    # Implementing Main Function.
    b_MarkerGenerator(_CHPED=args.chped, _OUT=args.o, _hg=args.hg, _dictionary_AA=args.dict_AA, _dictionary_SNPS=args.dict_SNPS)