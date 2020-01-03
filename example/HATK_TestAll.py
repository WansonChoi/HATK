# -*- coding: utf-8 -*-

import os, sys, re


def HATK_TestAll(_out=None, _control_flags=None):
    
    if bool(_out):

        if _out.endswith('/'):
            os.makedirs(_out, exist_ok=True)
        else:
            if bool(os.path.dirname(_out)):
                INTERMEDIATE_PATH = os.path.dirname(_out)
                os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
            else:
                os.makedirs(_out, exist_ok=True)

            _out = _out + '/'
    else:
        _out = ''

    ##### Control Flags

    if bool(_control_flags):

        [WHOLE,
        HLA2HPED,
        IMGT2SEQ,
        BMG,
        NOMENCLENAER,
        OMNIBUS,
        LOGISTIC,
        META,
        MANHATTAN,
        HEATMAP] = _control_flags

    else:

        # Manual Control
        WHOLE = 1
        HLA2HPED = 1
        IMGT2SEQ = 1
        BMG = 1
        NOMENCLENAER = 1
        OMNIBUS = 1
        LOGISTIC = 1
        META = 1
        MANHATTAN = 1
        HEATMAP = 1




    if WHOLE:

        ### 1. Whole Implementation
        print("\n\n### 1. Whole Implementation\n\n")

        command = \
        "python3 HATK.py \
            --variants example/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --hped example/wtccc_filtered_58C_RA.hatk.300+300.hped \
            --2field \
            --pheno example/wtccc_filtered_58C_RA.hatk.300+300.phe \
            --pheno-name RA \
            --out {OUT}MyHLAStudy/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --imgt 3320 \
            --hg 18 \
            --imgt-dir example/IMGTHLA3320 \
            --multiprocess 8 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyHLAStudy.log'))

        O = os.system(command)

    if HLA2HPED:

        ### 2. _0_HLA2HPED
        print("\n\n### 2. _0_HLA2HPED\n\n")

        # 1. xHLA
        command = \
        "python HATK.py \
            --hla2hped \
            --platform xHLA \
            --out {OUT}MyHLA2HPED_xHLA/MyxHLA \
            --rhped \
            example/HLA2HPED/xHLA/test.json \
            example/HLA2HPED/xHLA/test2.json \
            example/HLA2HPED/xHLA/test3.json > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyHLA2HPED_xHLA.log'))

        O = os.system(command)

        # 2. HIBAG
        command = \
        "python HATK.py \
            --hla2hped \
            --platform HIBAG \
            --out {OUT}MyHLA2HPED_HIBAG/MyHIBAG \
            --rhped \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-A.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-B.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-C.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DPA1.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DPB1.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DQA1.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DQB1.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DRB1.out > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyHLA2HPED_HIBAG.log'))

        O = os.system(command)


        command = \
        "python HATK.py \
            --hla2hped \
            --platform HIBAG \
            --out {OUT}MyHLA2HPED_HIBAG/MyHIBAG_withNA \
            --rhped \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-A.out \
            NA \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-C.out \
            NA \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DPB1.out \
            NA \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DQB1.out \
            example/HLA2HPED/HIBAG/HIBAG_TestResult.HLA-DRB1.out > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyHLA2HPED_HIBAG_withNA.log'))

        O = os.system(command)

        # 3. Axiom
        command = \
        "python HATK.py \
            --hla2hped \
            --platform AXIOM \
            --out {OUT}MyHLA2HPED_Axiom/MyAxiom \
            --rhped \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_A_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_B_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_C_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_DPA1_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_DPB1_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_DQA1_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_DQB1_Results.txt \
            example/HLA2HPED/Axiom/AxiomHLA_4dig_DRB1_Results.txt > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyHLA2HPED_Axiom.log'))

        O = os.system(command)

    if IMGT2SEQ:

        ### 3. _1_IMGT2Seq
        print("\n\n### 3. _1_IMGT2Seq\n\n")

        command = \
        "python3 HATK.py \
            --imgt2seq \
            --hg 18 \
            --imgt 3320 \
            --out {OUT}MyIMGT2Seq/ExamplePrefix.hg18.imgt3320 \
            --imgt-dir example/IMGTHLA3320 \
            --multiprocess 8 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyIMGT2Seq.log'))

        O = os.system(command)

    if BMG:

        ### 4. _2_bMarkerGenerator
        print("\n\n### 4. _2_bMarkerGenerator\n\n")

        # with '--variant'
        command = \
        "python3 HATK.py \
            --bmarkergenerator \
            --variants example/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --chped example/wtccc_filtered_58C_RA.hatk.300+300.imgt3320.2field.chped \
            --out {OUT}MybMarkerGenerator/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --hg 18 \
            --dict-AA example/RESULT_EXAMPLE/HLA_DICTIONARY_AA.hg18.imgt3320 \
            --dict-SNPS example/RESULT_EXAMPLE/HLA_DICTIONARY_SNPS.hg18.imgt3320 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MybMarkerGenerator.log'))

        O = os.system(command)

        # without '--variant'
        command = \
        "python3 HATK.py \
            --bmarkergenerator \
            --chped example/wtccc_filtered_58C_RA.hatk.300+300.imgt3320.2field.chped \
            --out {OUT}MybMarkerGenerator/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.novariants \
            --hg 18 \
            --dict-AA example/RESULT_EXAMPLE/HLA_DICTIONARY_AA.hg18.imgt3320 \
            --dict-SNPS example/RESULT_EXAMPLE/HLA_DICTIONARY_SNPS.hg18.imgt3320 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MybMarkerGenerator.novariants.log'))

        O = os.system(command)

    if NOMENCLENAER:

        ### 5. _3_NomenCleaner
        print("\n\n### 5. _3_NomenCleaner\n\n")

        # [1. Specifying the output nomenclature]
        command = \
        "python3 HATK.py \
            --nomencleaner \
            --hat example/RESULT_EXAMPLE/HLA_ALLELE_TABLE.imgt3320.hat \
            --hped example/wtccc_filtered_58C_RA.hatk.300+300.hped \
            --out {OUT}MyNomenCleaner/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.2field \
            --2field \
            --imgt 3320 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyNomenCleaner.2field.log'))
        O = os.system(command)

        # [2. Specifying NO output nomenclature]
        command = \
        "python3 HATK.py \
            --nomencleaner \
            --hat example/RESULT_EXAMPLE/HLA_ALLELE_TABLE.imgt3320.hat \
            --hped example/wtccc_filtered_58C_RA.hatk.300+300.hped \
            --out {OUT}MyNomenCleaner/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --imgt 3320 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyNomenCleaner.default.log'))
        O = os.system(command)

    if OMNIBUS:

        ### 6. _4-1_Omnibus Test
        print("\n\n### 6. _4-1_Omnibus Test\n\n")

        # Typical
        command = \
        "python3 HATK.py \
            --omnibus \
            --fam example/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.fam \
            --phased example/OmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.bgl.phased \
            --out {OUT}MyOmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.1 \
            --pheno example/wtccc_filtered_58C_RA.hatk.300+300.phe \
            --pheno-name RA > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyOmnibusTest.bgl_phased.log'))

        O = os.system(command)

        # with *.aa file.
        command = \
        "python3 HATK.py \
            --omnibus \
            --fam example/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.fam \
            --aa example/OmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.aa \
            --out {OUT}MyOmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.2 \
            --pheno example/wtccc_filtered_58C_RA.hatk.300+300.phe \
            --pheno-name RA > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyOmnibusTest.aa.log'))

        O = os.system(command)

        # with condition
        command = \
        "python3 HATK.py \
            --omnibus \
            --fam example/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.fam \
            --aa example/OmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.aa \
            --out {OUT}MyOmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --pheno example/wtccc_filtered_58C_RA.hatk.300+300.phe \
            --pheno-name RA \
            --condition AA_DRB1_96 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyOmnibusTest.condition.log'))

        O = os.system(command)

    if LOGISTIC:
        ### 7. _4-2_Logistic Regression
        print("\n\n### 7. _4-2_Logistic Regression\n\n")

        command = \
        "python HATK.py \
            --logistic \
            --input example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
            --pheno-name RA \
            --pheno example/wtccc_filtered_58C_RA.hatk.300+300.phe \
            --out {OUT}MyLogisticReg/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyLogisticReg.log'))

        O = os.system(command)

    if META:
        ### 8. _4-3_MetaAnalysis
        print("\n\n### 8. _4-3_MetaAnalysis\n\n")

        command = \
        "python HATK.py \
            --metaanalysis \
            --s1-logistic-result example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.assoc.logistic \
            --s2-logistic-result example/MetaAnalysis/RESULT_EXAMPLE_wtccc_filtered_NBS_CD.hatk.300+300.chr6.hg18.assoc.logistic \
            --out {OUT}MyMetaAnalysis/RESULT_EXAMPLE_wtccc_filtered.RAvsCD > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyMetaAnalysis.plain.log'))

        O = os.system(command)

        command = \
        "python HATK.py \
            --metaanalysis \
            --s1-logistic-result example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.assoc.logistic \
            --s1-bim example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.bim \
            --s2-logistic-result example/MetaAnalysis/RESULT_EXAMPLE_wtccc_filtered_NBS_CD.hatk.300+300.chr6.hg18.assoc.logistic \
            --s2-bim example/MetaAnalysis/RESULT_EXAMPLE_wtccc_filtered_NBS_CD.hatk.300+300.chr6.hg18.bim \
            --out {OUT}MyMetaAnalysis/RESULT_EXAMPLE_wtccc_filtered.RAvsCD2.flip > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyMetaAnalysis.flip.log'))

        O = os.system(command)

    if MANHATTAN:
        ### 9. _5-1_Manhattan
        print("\n\n### 9. _5-1_Manhattan\n\n")

        command = \
        "python HATK.py \
            --manhattan \
            --ar example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.assoc.logistic \
            --imgt 3320 \
            --hg 18 \
            --out {OUT}MyManhattan/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyManhattan.log'))

        O = os.system(command)

        
        command = \
        "python HATK.py \
            --manhattan \
            --ar example/OmnibusTest/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.RA.NA.omnibus \
            --imgt 3320 \
            --hg 18 \
            --out {OUT}MyManhattan/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.RA.NA.omnibus \
            --HLA A DRB1 DQA1 DQB1 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyManhattan.omnibus.log'))

        O = os.system(command)

    if HEATMAP:
        ### 10. _5-2_Heatmap
        print("\n\n### 10. _5-2_Heatmap\n\n")

        command = \
        "python HATK.py \
            --heatmap \
            --ar example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.assoc.logistic \
            --HLA A \
            --maptable example/RESULT_EXAMPLE/HLA_MAPTABLE_A.hg18.imgt3320.txt \
            --out {OUT}MyHeatmap/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 > {LOG}".format(OUT=_out, LOG=os.path.join(_out, 'MyHeatmap.log'))

        O = os.system(command)




if __name__ == '__main__':

    """
    
    Final module to test the BASH execution of each module of HATK.
    Usage examples in each README file will be implemented sequentially.

    1. Whole implementation
    2. _0_HLA2HPED
    3. _1_IMGT2Seq
    4. _2_bMarkerGenerator
    5. _3_NomenCleaner
    6. _4-1_OminbusTest
    7. _4-2_LogisticRegr
    8. _4-3_Meta
    9. _5-1_Manhattan
    10. _5-2_Heatmap
    
    """

    [_out, _control_flags] = sys.argv[1:]

    if _control_flags == 'All' or _control_flags == 'all':
        _control_flags = [1,1,1,1,1, 1,1,1,1,1]

    elif bool(re.match(pattern=r'^\d(,\d){9}$', string=_control_flags)):
        _control_flags = list(map(lambda x : int(x), _control_flags.split(',')))
        print(_control_flags)
    else:
        print("Error. Wrong _control_flags!")
        sys.exit()

    # print(_control_flags)

    # HATK_TestAll(_out)

    # _control_flag = [0,0,0,0,0,0,0,0,0,1]

    HATK_TestAll(_out=_out, _control_flags=_control_flags)
    # HATK_TestAll()
