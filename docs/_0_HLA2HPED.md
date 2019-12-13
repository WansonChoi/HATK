# HLA2HPED

## (1) What is HLA2HPED

`HLA2HPED` is a module to merge and convert the output from other HLA imputation or inference software to `HLA PED(HPED)` file so that researchers can use `HATK` directly using the output of other HLA software. Currently, this module supports **xHLA**, **HIBAG**, and **Axiom** software and it is planned to make HLA2HPED to cover more HLA imputation/inference software output.

`HLA PED(HPED)` file format is similar to the PLINK ped file but consists of 22(6 + 8*2) columns. Left 6 columns are exactly same as PLINK ped file ('Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype'). Other 16 columns are Individual's HLA type information (2 HLA alleles for each HLA gene) of **8 HLA genes(A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1; in this order)**. If some allele values are **not accessible(NA)**, then fill them with just '0'.


![README-HLA2HPED-HPED](./img/README_0_HLA2HPED-HPED.png)
(an example of HPED file.)

<br>

## (2) Usage Examples.

1. xHLA

The output of xHLA is assumed to be multiple JSON files each of which represents a single sample's HLA type information. The next photo is an example of the JSON output file of xHLA. Multiple JSON files like this will be given to HLA2HPED with the argument '--rhped'.

![README-HLA2HPED-xHLA_example](./img/README_0_HLA2HPED-xHLA_example.png)

```
$ python HATK.py \
    --hla2hped \
    --platform xHLA \
    --out MyHLA2HPED_xHLA/test \
    --rhped example/HLA2HPED/xHLA/test{,2,3}.json
```

The 1st column, Family_ID, will be set by the sample id, the 2nd column value by default.

> The usage like the last line '--rhped example/HLA2HPED/xHLA/test{,2,3}.json' is called **brace expansion** in Unix/Linux. That line is exactly same as '--rhped example/HLA2HPED/xHLA/test.json example/HLA2HPED/xHLA/test2.json example/HLA2HPED/xHLA/test3.json'.

<br>

2. HIBAG

HIBAG output is usually generated as files split by HLA. The next photo is an example of one output file which contains HLA-A imputation.

Because HATK focuses on the **8 HLA genes(A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1)**, HLA2HPED will work by taking 8 output files of each HLA gene. **Please pass the files of the 8 HLA genes with that order**.

![README-HLA2HPED-HIBAG_example](./img/README_0_HLA2HPED-HIBAG_example.png)

```
$ python HATK.py \
    --hla2hped \
    --platform HIBAG \
    --out MyHLA2HPED_HIBAG/MyHIBAG \
    --rhped example/HLA2HPED/HIBAG/HIBAG_TestResult.{A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1}.out

```


If some files of HLA genes are not accessible, then allocate 'NA' to that position of the '--rhped'.

```
$ python HATK.py \
    --hla2hped \
    --platform HIBAG \
    --out MyHLA2HPED_HIBAG/MyHIBAG \
    --rhped NA example/HLA2HPED/HIBAG/HIBAG_TestResult.B.out NA NA NA NA NA NA

# Using only HLA-B gene information
```

<br>

3. Axiom

The output of Axiom is also generated as split files of respective HLA. Likewise HIBAG, **Please pass the files of the 8 HLA genes with the above order**.

```
$ python HATK.py \
    --hla2hped \
    --platform AXIOM \
    --out MyHLA2HPED_Axiom/MyAxiom \
    --rhped example/HLA2HPED/Axiom/AxiomHLA_4dig_{A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1}_Results.txt
```
