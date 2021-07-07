#-*- coding: utf-8 -*-

import os, sys, re

from IMGT2Seq.__main__ import HATK_IMGT2Seq
from bMarkerGenerator.bMarkerGenerator import HATK_bMarkerGenertor
from HLA2HPED.HLA2HPED import HATK_HLA2HPED
from NomenCleaner.NomenCleaner import HATK_NomenCleaner
from HLA_Analysis.HLA_Analysis import HATK_LogisticRegression
from HLA_Heatmap.heatmap import HATK_Heatmap
from HLA_Manhattan.manhattan import HATK_Manhattan
from HLA_Analysis.HLA_Analysis import HATK_OmibusTest
from HLA_Analysis.src.MetaAnalysis.HLA_MetaAnalysis import HATK_MetaAnalysis



std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



class HLA_Study(object):


    def __init__(self, _args):

        """

        Main implementation body of HATK.
        HATK will work mainly by creating `HLA_Study` instance.

        """

        # print(_args)


        """
        Preprocessing major arguments
        (1) _args.input
        (2) _args.out
        (3) _args.imgt2sequence, ..., etc. (Module flags)
        """

        # (1) _args.input

        if _args.input:

            """
            * --input := the common prefix to point out next files
                (2) normal_SNPs (*.bed, *.bim, *.fam)
                (3) phenotype file(*.phe, *.pheno)
                (4) covariate file(*.covar, *.cov)
                (6) reference allele (*.refallele, *.ref)
            """

            # `normal_SNPs` => just as prefix.
            if not _args.variants:
                _args.variants = _args.input if os.path.exists(_args.input + ".bed") and os.path.exists(_args.input + ".bim") and os.path.exists(_args.input + ".fam") else None
            if not _args.pheno:
                _args.pheno = _args.input + ".phe" if os.path.exists(_args.input + ".phe") else _args.input + ".pheno" if os.path.exists(_args.input + ".pheno") else None
            if not _args.covar:
                _args.covar = _args.input + ".covar" if os.path.exists(_args.input + ".covar") else _args.input + ".cov" if os.path.exists(_args.input + ".cov") else None

            if not _args.reference_allele:
                _args.reference_allele = _args.input + ".refallele" if os.path.isfile(_args.input + ".refallele") else _args.input + ".ref" if os.path.exists(_args.input + ".ref") else None


            if not (_args.variants or _args.pheno or _args.covar or _args.reference_allele):
                print(std_ERROR_MAIN_PROCESS_NAME + "None of files related to '--variants', '--pheno', '--covar', '--reference-allele' arguments were overriden.\n"
                                                    "Please check '--input' argument again.")
                sys.exit()


            # *.aa and *.bgl.phased (related to Ominibus Test.)
            if not _args.aa and os.path.exists(_args.input+'.aa'):
                _args.aa = _args.input+'.aa'

            if not _args.phased and os.path.exists(_args.input+'.bgl.phased'):
                _args.phased = _args.input+'.bgl.phased'


        # (2) _args.out

        if not _args.out:
            print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not given. Please check it again.\n'.format("--out"))
            sys.exit()
        else:
            _args.out = _args.out if not _args.out.endswith('/') else _args.out.rstrip('/')
            if bool(os.path.dirname(_args.out)): os.makedirs(os.path.dirname(_args.out), exist_ok=True)




        # (3) Module Flags

        flag_MODULEs = [_args.imgt2seq,
                        _args.bmarkergenerator,
                        _args.hla2hped,
                        _args.nomencleaner,
                        _args.logistic, _args.omnibus, _args.metaanalysis,
                        _args.heatmap,
                        _args.manhattan]

        f_sum = sum(flag_MODULEs)
        # print("f_sum : {}".format(f_sum))


        if f_sum == 0:



            ########## < [0] IMGT2Seq > ##########

            myIMGT2Seq = HATK_IMGT2Seq(_args.imgt, _args.hg, _args.out,
                                       _args.oneF, _args.twoF, _args.threeF, _args.fourF, _args.Ggroup, _args.Pgroup,
                                       _no_indel=_args.no_indel, _multiprocess=_args.multiprocess,
                                       _save_intermediates=_args.save_intermediates,
                                       _imgt_dir=_args.imgt_dir)

            if myIMGT2Seq:
                print(std_MAIN_PROCESS_NAME + "IMGT2Seq result : \n{}".format(myIMGT2Seq))

                _args.dict_AA, _args.dict_SNPS, _args.hat, _args.maptable = myIMGT2Seq.getResult()

            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Failed to preprocess HLA sequence information.")
                print(std_MAIN_PROCESS_NAME + "IMGT2Seq result : \n{}".format(myIMGT2Seq))
                sys.exit()




            ########## < [1] Checking HLA Typing information > ##########

            if _args.chped:
                pass


            elif _args.hped:

                print(std_MAIN_PROCESS_NAME + "Given HPED file('{}') is to be processed by NomenCleaner.".format(_args.hped))

                myNomenCleaner = HATK_NomenCleaner(_args.hped, _args.hat, _args.imgt, _args.out,
                                                   __oneF=_args.oneF, __twoF=_args.twoF, __threeF=_args.threeF, __fourF=_args.fourF,
                                                   __Ggroup=_args.Ggroup, __Pgroup=_args.Pgroup,
                                                   __f_NoCaption=_args.NoCaption, __leave_NotFound=_args.leave_NotFound)

                _args.chped = myNomenCleaner.getResult()

                if _args.chped == -1:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Failed to process HPED file('{}') by NomenCleaner.".format(_args.hped))
                    sys.exit()


            elif _args.rhped:

                ### Step 1 : rhped to hped (HLA2HPED)
                print(std_MAIN_PROCESS_NAME + "Given raw HPED file(s)('{}') is/are to be processed by HLA2HPED.".format(_args.rhped))

                myHLA2HPED = HATK_HLA2HPED(_args.rhped, _args.out, _args.platform)

                _args.hped = myHLA2HPED.getResult()

                if _args.hped == -1:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Failed to process raw HPED file('{}') by HLA2HPED.")
                    sys.exit()


                ### Step 2 : hped to chped (NomenCleaner)
                print(std_MAIN_PROCESS_NAME + "Generated HPED file('{}') is to be processed by NomenCleaner.".format(_args.hped))

                myNomenCleaner = HATK_NomenCleaner(_args.hped, _args.hat, _args.imgt, _args.out,
                                                   __oneF=_args.oneF, __twoF=_args.twoF, __threeF=_args.threeF, __fourF=_args.fourF,
                                                   __Ggroup=_args.Ggroup, __Pgroup=_args.Pgroup,
                                                   __f_NoCaption=_args.NoCaption, __leave_NotFound=_args.leave_NotFound)

                _args.chped = myNomenCleaner.getResult()


                if _args.chped == -1:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Failed to process HPED('{}') by NomenCleaner.".format(_args.hped))
                    sys.exit()

            else:

                print(std_ERROR_MAIN_PROCESS_NAME + "No HLA type information was given.\n"
                                                    "Please check '--chped', '--hped', or '--rhped' arguments again.")
                sys.exit()

            # Priority : *.chped > *.hped > *.rhped




            ########## < [2] bMarkerGenerator > ##########

            mybMarkers = HATK_bMarkerGenertor(_args.chped, _args.out, _args.hg, _args.dict_AA, _args.dict_SNPS,
                                              _variants=_args.variants, __save_intermediates=_args.save_intermediates,
                                              _p_src="bMarkerGenerator/src")

            _args.variants = mybMarkers.getReuslt()

            if _args.variants == -1:
                print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate Binary Markers for '{}'".format(_args.chped))
                sys.exit()
            else:
                print(std_MAIN_PROCESS_NAME + "bMarkerGenerator result(Prefix) : \n{}".format(mybMarkers.getReuslt()))




            ########## < [3] Association Test (Logistic Regression) > ##########

            # HARD_CODING = 'tests/_0_WholeImplementation/20190510_58C_NBS_test/58C_NBS_test'
            # _args.variants = HARD_CODING

            myLogistcRegression = HATK_LogisticRegression(_args.variants, _args.out,
                                                          _phe=_args.pheno, _phe_name=_args.pheno_name,
                                                          _covar=_args.covar, _covar_name=_args.covar_name,
                                                          _condition=_args.condition,
                                                          _condition_list=_args.condition_list,
                                                          _ref_allele=_args.reference_allele)

            _args.assoc_result = [myLogistcRegression.getResult()] # `_args.assoc_result` is supposed to be a list.

            if _args.assoc_result[0] == -1:
                print(std_ERROR_MAIN_PROCESS_NAME + "Failed to perform logistic regression test on the Binary Markers('{}').".format(_args.variants))
                sys.exit()
            else:
                print(std_MAIN_PROCESS_NAME + "Logistic Regression result : \n{}".format(myLogistcRegression.getResult()))




            ########## < [4] Manhattan Plot > ##########

            myManhattan = HATK_Manhattan(_args.assoc_result, _args.out+".manhattan", _args.hg,
                                         _point_col=_args.point_color, _top_color=_args.top_color,
                                         _point_size=_args.point_size, _yaxis_unit=_args.yaxis_unit,
                                         _p_src="HLA_Manhattan/src", _p_data="HLA_Manhattan/data", _HLA=None)

            if myManhattan.getResult() == -1:
                print(std_WARNING_MAIN_PROCESS_NAME + "Failed to plot HLA Manhattan plot.")
            else:
                print(std_MAIN_PROCESS_NAME + "Manhattan result : \n{}".format(myManhattan.getResult()))




            ########## < [5] Heatmap Plot > ##########

            # Plotting Heatmaps over 8 HLA genes

            Heatmap_status = {_hla : -1 for _hla in HLA_names}

            for i in range(0, len(HLA_names)):

                t_HLA = HLA_names[i]
                t_maptable = _args.maptable[t_HLA]

                myHeatmap = HATK_Heatmap(t_HLA, _args.out + ".HLA_{}.heatmap".format(t_HLA), t_maptable,
                                         _args.assoc_result[0], __save_intermediates=_args.save_intermediates,
                                         _p_src="HLA_Heatmap/src", _p_data="HLA_Heatmap/data")

                if myHeatmap.getResult() == -1:
                    print(std_WARNING_MAIN_PROCESS_NAME + "Failed to plot Heatmap of HLA {} gene.".format(t_HLA))
                else:
                    Heatmap_status[t_HLA] = myHeatmap.getResult()

            print(std_MAIN_PROCESS_NAME + "Heatmap result : ")
            for i in range(0, len(HLA_names)):
                print(" {} : {}".format(HLA_names[i], Heatmap_status[HLA_names[i]]))





        elif f_sum == 1:

            ##### Partial Implementation #####

            if _args.imgt2seq:

                ### IMGT2Seq
                myIMGT2Seq = HATK_IMGT2Seq(_args.imgt, _args.hg, _args.out,
                                           _args.oneF, _args.twoF, _args.threeF, _args.fourF, _args.Ggroup, _args.Pgroup,
                                           _no_indel=_args.no_indel, _multiprocess=_args.multiprocess,
                                           _save_intermediates=_args.save_intermediates,
                                           _imgt_dir=_args.imgt_dir)

                print(std_MAIN_PROCESS_NAME + "IMGT2Seq result : \n{}".format(myIMGT2Seq))

            elif _args.bmarkergenerator:

                ### bMarkerGenerator
                mybMarkers = HATK_bMarkerGenertor(_args.chped, _args.out, _args.hg, _args.dict_AA, _args.dict_SNPS,
                                                  _variants=_args.variants, __save_intermediates=_args.save_intermediates,
                                                  _p_src="bMarkerGenerator/src")

                print(std_MAIN_PROCESS_NAME + "bMarkerGenerator result(Prefix) : \n{}".format(mybMarkers.getReuslt()))

            elif _args.hla2hped:

                ### HLA2HPED
                myHLA2HPED = HATK_HLA2HPED(_args.rhped, _args.out, _args.platform)

                print(std_MAIN_PROCESS_NAME + "HLA2HPED result : \n{}".format(myHLA2HPED.getResult()))

            elif _args.nomencleaner:

                ### NomenCleaner
                myNomenCleaner = HATK_NomenCleaner(_args.hped, _args.hat, _args.imgt, _args.out,
                                                   __oneF=_args.oneF, __twoF=_args.twoF, __threeF=_args.threeF, __fourF=_args.fourF,
                                                   __Ggroup=_args.Ggroup, __Pgroup=_args.Pgroup,
                                                   __f_NoCaption=_args.NoCaption, __leave_NotFound=_args.leave_NotFound)

                print(std_MAIN_PROCESS_NAME + "NomenCleaner result : \n{}".format(myNomenCleaner.getResult()))

            elif _args.logistic:

                ### Logistic Regression
                myLogistcRegression = HATK_LogisticRegression(_args.variants, _args.out,
                                                              _phe=_args.pheno, _phe_name=_args.pheno_name,
                                                              _covar=_args.covar, _covar_name=_args.covar_name,
                                                              _condition=_args.condition, _condition_list=_args.condition_list,
                                                              _ref_allele=_args.reference_allele)

                print(std_MAIN_PROCESS_NAME + "Logistic Regression result : \n{}".format(myLogistcRegression.getResult()))

            elif _args.omnibus:

                ### Omnibus Test
                myOminubus = HATK_OmibusTest(_args.out, _args.fam, _phe=_args.pheno, _phe_name=_args.pheno_name,
                                             _covar=_args.covar, _covar_name=_args.covar_name,
                                             _condition=_args.condition, _condition_list=_args.condition_list,
                                             _bgl_phased=_args.phased, _aa=_args.aa)

            elif _args.metaanalysis:

                ### Meta Analysis
                myMeta = HATK_MetaAnalysis(_args.s1_logistic_result, _args.s2_logistic_result, _args.out,
                                           _study1_m=_args.s1_bim, _study2_m=_args.s2_bim)
                print(std_MAIN_PROCESS_NAME + "MetaAnalysis result : \n{}".format(myMeta))

            elif _args.heatmap:

                ### HLA Heatmap
                myHeatmap = HATK_Heatmap(_args.HLA, _args.out, _args.maptable, _args.ar,
                                         __save_intermediates=_args.save_intermediates, _p_src="HLA_Heatmap/src",
                                         _p_data="HLA_Heatmap/data")

                print(std_MAIN_PROCESS_NAME + "Heatmap result : \n{}".format(myHeatmap.getResult()))

            elif _args.manhattan:

                ### HLA Manhattan
                myManhattan = HATK_Manhattan(_args.ar, _args.out, _args.hg,
                                             _point_col=_args.point_color, _top_color=_args.top_color,
                                             _point_size=_args.point_size, _yaxis_unit=_args.yaxis_unit,
                                             _p_src="HLA_Manhattan/src", _p_data="HLA_Manhattan/data", _HLA=_args.HLA)

                if myManhattan.getResult() == -1:
                    print(std_WARNING_MAIN_PROCESS_NAME + "Failed to plot HLA Manhattan plot.")
                else:
                    print(std_MAIN_PROCESS_NAME + "Manhattan result : \n{}".format(myManhattan.getResult()))


        elif f_sum > 1:

            # When more than 1 module is asked to be implemented.
            # Currently, it will be classified to wrong implementation.
            print(std_ERROR_MAIN_PROCESS_NAME + "Each module must be implemented separately.")
            sys.exit()
        else:
            pass
