# Logistic Regression

## (1) Introduction

Logistic regression is one of the most widely used association test methods for the `Genome-Wide Association Test(GWAS)`. HATK provides the module to implement a logistic regression that is optimized to perform `HLA fine-mapping`.

<br>

## (2) Usage Example

```
$ python HATK.py \
    --logistic \
    --input example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --pheno-name RA \
    --pheno example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.phe \
    -o MyLogisticReg/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18
```