# MetaAnalysis

## (1) Introduction


## (2) Usage examples

```
$ python HATK.py \
    --metaanalysis \
    --s1-logistic-result example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.assoc.logistic \
    --s2-logistic-result example/MetaAnalysis/RESULT_EXAMPLE_wtccc_filtered_NBS_CD.hatk.300+300.chr6.hg18.assoc.logistic \
    --out MyMetaAnalysis/RESULT_EXAMPLE_wtccc_filtered.RAvsCD

    
```

*.bim files of each Logistic regression result can be used to apply flipped relationship of each shared markers.

```
$ python HATK.py \
    --metaanalysis \
    --s1-logistic-result example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.assoc.logistic \
    --s1-bim example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.bim \
    --s2-logistic-result example/MetaAnalysis/RESULT_EXAMPLE_wtccc_filtered_NBS_CD.hatk.300+300.chr6.hg18.assoc.logistic \
    --s2-bim example/MetaAnalysis/RESULT_EXAMPLE_wtccc_filtered_NBS_CD.hatk.300+300.chr6.hg18.bim \
    --out MyMetaAnalysis/RESULT_EXAMPLE_wtccc_filtered.RAvsCD2

    
```
