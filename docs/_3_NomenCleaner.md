# NomenCleaner

## Usage example.

NomenCleaner works based on one output file from IMGT2Seq, 'Integrated Allele Table(*.iat)' file.


```
python3 HATK.py \
    --nomencleaner \
    -iat /Users/wansun/Git_Projects/HATK/tests/_1_IMGT2Sequence/20190303_hg19_imgt3320/HLA_INTEGRATED_ALLELE_TABLE.hg19.imgt3320.iat \
    -hped /Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/wtccc_filtered_58C_NBS_RA_T1D.hped \
    -o /Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/wtccc_filtered_58C_NBS_RA_T1D \
    --4field \
    -imgt 3320

```