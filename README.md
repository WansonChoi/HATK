# HATK_2nd

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


[Using Dockers]

1. pull hatk docker image.

docker pull docker pull wschoibhlab/hatk


2. run the image in docker container. (in the project directory)

docker run -v `pwd`:`pwd` -w `pwd` wschoibhlab/hatk \
    python3 HATK.py \
    --input example/wtccc_filtered_58C_NBS_RA_T1D \
    --hped example/wtccc_filtered_58C_NBS_RA_T1D.hped \
    --out tests/_0_wholeProcess/20190625/20190626_example \
    --pheno-name RA \
    -imgt 3320 \
    -hg 19 \
    --imgt-dir IMGT2Seq/IMGTHLA3320 \
    --multiprocess


