# Manhattan Plot

## (1) Introduction
The Manhattan plot is an indispensable method to visualize the result of the association test. HATK provides the module to plot this.

<br>

## (2) Usage Example

```
$ python HATK.py \
    --manhattan \
    -ar example/RESULT_EXAMPLE/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.assoc.logistic \
    -imgt 3320 \
    -hg 18 \
    -o MyManhattan/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18
```


## (3) Result Example
![Manhattan_example](img/README_5-1_Manhattan_example.png)

- Yellow: HLA marker (ex. HLA_A_0101)
- Red: Amino acid marker (ex. AA_A_9_30018537_F, AA_A_-15_30018338)
- Dark grey: Intragenic DNA base pair marker (ex. SNPS_A_30018457, SNPS_A_30018461_A)
- Ligtht grey: Intergenic dbSNP marker (ex. rs41557221)
