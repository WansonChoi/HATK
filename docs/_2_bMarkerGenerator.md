# bMarkerGenerator

## Usage example

(1) 
```
python3 HATK.py \
    --bmarkergenerator \
    --variants example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300 \
    --chped example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.imgt3320.4field.chped \
    --out tests/_2_b_MarkerGenerator/20190829_RA_300_300/20190702_wtccc_filtered_58C_RA.hatk.58C_RA.300+300 \
    -hg 19 \
    --dict-AA tests/_1_IMGT2Sequence/20190829_hg19_imgt3320/HLA_DICTIONARY_AA.hg19.imgt3320 \
    --dict-SNPS tests/_1_IMGT2Sequence/20190829_hg19_imgt3320/HLA_DICTIONARY_SNPS.hg19.imgt3320
```

(2) Just HLA markers. (No SNP variant markers)

Don't include '--variant' argument in above useage example.

```
python3 HATK.py \
    --bmarkergenerator \
    --chped example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.imgt3320.4field.chped \
    --out tests/_2_b_MarkerGenerator/20190829_RA_300_300/20190702_wtccc_filtered_58C_RA.hatk.58C_RA.300+300 \
    -hg 19 \
    --dict-AA tests/_1_IMGT2Sequence/20190829_hg19_imgt3320/HLA_DICTIONARY_AA.hg19.imgt3320 \
    --dict-SNPS tests/_1_IMGT2Sequence/20190829_hg19_imgt3320/HLA_DICTIONARY_SNPS.hg19.imgt3320
```
