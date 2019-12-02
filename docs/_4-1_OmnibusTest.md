# Omnibus Test

## (1) What is Ominibus Test?

Under construction.

## (2) Setting condition and Phasing

Under construction.

## (3) How to perform phasing.

Under construction.

We strongly recommend using Beagle(v3.0.1) software to phase target data.

> Phasing output generated from BEAGLE with the version 4.x.x and 5.x.x can be used in this Omnibus Test. However, BEAGLE with those version doesn't allow markers to have same base position while HATK generates markers which share same base position. In other words, **It is recommended to just use BEAGLE with the version 3.x.x to generate the phased beagle file of the HATK output**.

## (4) Usage Examples

When no condition(i.e. when interested in only one amino acid locus) is given, phasing isn't needed. However

```
$ python3 HATK.py \
    --omnibus \
    --fam example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.fam \
    --phased example/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.bgl.phased \
    -o MyOmnibusTest/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --pheno example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.phe \
    --pheno-name RA

```

The input '\*.bgl.phased' file is first processed to '\*.aa' file, which is genuine input of Omnibus Test. Once the '\*.aa' file of the input '\*.bgl.phased' has been generated, then you don't need to repeat this process. The user can skip this with '--aa' argument.

```
$ python3 HATK.py \
    --omnibus \
    --fam example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.fam \
    --aa example/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.aa \
    -o MyOmnibusTest/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --pheno example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.phe \
    --pheno-name RA

```
This 2nd example will generate the same output as that of the 1st example.



If you want to set another marker(possibly another amino acid position) as covariate in the omnibus test, you can set this with '--condition' argument. If you want to pass more than 1 condition, then pass them as 'comma-separated (single) string'.

```
$ python3 HATK.py \
    --omnibus \
    --fam example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.fam \
    --aa example/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.aa \
    -o MyOmnibusTest/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --pheno example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.phe \
    --pheno-name RA \
    --condition AA_A_-15
    # --condition AA_A_-15,AA_DRB1_166 (when passing multiple coditions.)

```

(Tip) To sort the output of the omnibus test on P-value, Use the below bash command.
```
$ sort -gk 5 MyOmnibusTest/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18.RA.NA.omnibus
```