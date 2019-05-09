#-*- coding: utf-8 -*-

import os, sys, re


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
        Preprocessing arguments
        (1) _args.input
        (2) _args.out
        (3) _args.imgt2sequence, ..., etc. (Module flags)
        """

        # (1) _args.input

        if _args.input:

            """
            * --input := the common prefix to point out next files
                (1) *.hped or *.chped
                (2) normal_SNPs
                (3) phenotype file(*.phe)
                (4) covariate file(*.covar)
                (5) condition (*.condvars)
                (6) reference allele (*.refallele)
            """
            if not _args.hped:
                _args.hped = _args.input + ".hped" if os.path.exists(_args.input + ".hped") else None
            if not _args.chped:
                _args.chped = _args.input + ".chped" if os.path.exists(_args.input + ".chped") else None

            # `normal_SNPs` => just as prefix.
            if not _args.variants:
                _args.variants = _args.input if os.path.exists(_args.input + ".bed") and os.path.exists(_args.input + ".bim") and os.path.exists(_args.input + ".fam") else None
            if not _args.pheno:
                _args.pheno = _args.input + ".phe" if os.path.exists(_args.input + ".phe") else _args.input + ".pheno" if os.path.exists(_args.input + ".pheno") else None
            if not _args.covar:
                _args.covar = _args.input + ".covar" if os.path.exists(_args.input + ".covar") else _args.input + ".cov" if os.path.exists(_args.input + ".cov") else None

            if not _args.condition_list:
                _args.condition_list = _args.input + ".condvars" if os.path.isfile(_args.input + ".condvars") else None
            if not _args.reference_allele:
                _args.reference_allele = _args.input + ".refallele" if os.path.isfile(_args.input + ".refallele") else _args.input + ".ref" if os.path.exists(_args.input + ".ref") else None



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
                        _args.logistic, _args.omnibus, _args.meta_analysis,
                        _args.heatmap,
                        _args.manhattan]

        f_sum = sum(flag_MODULEs)
        # print("f_sum : {}".format(f_sum))


        if f_sum == 0:

            ##### Whole Implementation #####
            pass

        elif f_sum == 1:

            ##### Partial Implementation #####

            if _args.imgt2seq:

                ### IMGT2Seq
                from IMGT2Seq.IMGT2Seq import HATK_IMGT2Seq

                myIMGT2Seq = HATK_IMGT2Seq(_args.imgt, _args.hg, _args.out,
                                           _no_indel=_args.no_indel, _multiprocess=_args.multiprocess,
                                           _save_intermediates=_args.save_intermediates,
                                           _imgt_dir=_args.imgt_dir)

                print(std_MAIN_PROCESS_NAME + "IMGT2Seq results : \n{}".format(myIMGT2Seq))

            elif _args.bmarkergenerator:

                ### bMarkerGenerator
                from bMarkerGenerator.bMarkerGenerator import HATK_bMarkerGenertor

                mybMarkers = HATK_bMarkerGenertor(_args.chped, _args.out, _args.hg, _args.dict_AA, _args.dict_SNPS,
                                                  _variants=_args.variants, __save_intermediates=_args.save_intermediates,
                                                  _p_src="bMarkerGenerator/src")

                print(std_MAIN_PROCESS_NAME + "bMarkerGenerator result(Prefix) : \n{}".format(mybMarkers.getReuslts()))

            elif _args.hla2hped:

                ### HLA2HPED
                from HLA2HPED.HLA2HPED import HATK_HLA2HPED

                myHLA2HPED = HATK_HLA2HPED(_args.rhped, _args.out, _args.platform)

                print(std_MAIN_PROCESS_NAME + "HLA2HPED result : \n{}".format(myHLA2HPED.getResults()))

            elif _args.nomencleaner:

                ### NomenCleaner
                from NomenCleaner.NomenCleaner import HATK_NomenCleaner

                myNomenCleaner = HATK_NomenCleaner(_args.iat, _args.imgt, _args.out,
                                                   _args.hped, _args.hped_G, _args.hped_P,
                                                   _args.oneF, _args.twoF, _args.threeF, _args.fourF,
                                                   _args.G_group, _args.P_group, _args.old_format,
                                                   __f_NoCaption=_args.NoCaption, __leave_NotFound=_args.leave_NotFound)

                print(std_MAIN_PROCESS_NAME + "NomenCleaner result : \n{}".format(myNomenCleaner.getResults()))

            elif _args.logistic:

                ### Logistic Regression
                pass
            elif _args.omnibus:

                ### Omnibus Test
                pass
            elif _args.meta_analysis:

                ### Meta Analysis
                pass
            elif _args.heatmap:

                ### HLA Heatmap
                from HLA_Heatmap.heatmap import HATK_Heatmap

                myHeatmap = HATK_Heatmap(_args.HLA, _args.out, _args.maptable, _args.assoc_result,
                                         __as4field=_args.as4field, __save_intermediates=_args.save_intermediates,
                                         _p_src="HLA_Heatmap/src", _p_data="HLA_Heatmap/data")

                print(std_MAIN_PROCESS_NAME + "Heatmap results : \n{}".format(myHeatmap.getResults()))

            elif _args.manhattan:

                ### HLA Manhattan
                from HLA_Manhattan.manhattan import HATK_Manhattan

                myManhattan = HATK_Manhattan(_args.assoc_result, _args.out, _args.hg,
                                             _point_col=_args.point_color, _top_color=_args.top_color,
                                             _point_size=_args.point_size, _yaxis_unit=_args.yaxis_unit,
                                             _p_src="HLA_Manhattan/src", _p_data="HLA_Manhattan/data")

                print(std_MAIN_PROCESS_NAME + "Manhattan results : \n{}".format(myManhattan.getResults()))


        elif f_sum > 1:

            # When more than 1 module is asked to be implemented.
            # Currently, it will be classified to wrong implementation.
            print(std_ERROR_MAIN_PROCESS_NAME + "Each module must be implemented separately.")
            sys.exit()
        else:
            pass