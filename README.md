# HATK


## (1) What is HATK?

HATK is ... .


## (2) Installation

We strongly recommend below two ways to install HATK.

1. Use 'Anaconda'(https://www.anaconda.com/). (or 'Miniconda'(https://docs.conda.io/en/latest/miniconda.html))

2. Use Docker(https://www.docker.com/).


```
docker pull wschoibhlab/hatk
```


## (3) Usage

2. run this image (in the project directory.)

```
docker run -v `pwd`:/HATK wschoibhlab/hatk \
    python3 HATK.py \
    --variants example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300 \
    --hped example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.hped \
    --pheno example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.phe \
    --pheno-name RA \
    --out MyHLAStudy/MyHLAStudy_wtccc_filtered_58C_RA.hatk.58C_RA.300+300 \
    -imgt 3320 \
    -hg 19 \
    --imgt-dir example/IMGTHLA3320 \
    --multiprocess
```
## (4) Citation

## (5) Lincense






<!-- comment 
## \<History\>

2nd Repository for HATK project.

(2018. 8. 2.)
Remote repository has been moved from Bitbucket to Github.


(2018. 12. 19.)
The branch 'b_20181219' has been created to
	(1) introduce logging system,
	(2) optimize and enhance the general performance,

and etc. 


(2019. 1. 10.)
The core engine modules("HLAtoSequences.py", "encodeVariants.py", "encodeHLA.py") are reworked urgently to solve the memroy usage problem(It was found to use maximum 64G RAM apporximately maybe due to Pandas).

The rework was primarily done in the work with Yang Luo in the repository of "MakeReference_v2" and the finalized rework outputs are applied to this project.
-->