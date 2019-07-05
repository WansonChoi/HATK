# HATK


## (1) What is HATK?

Fine mapping of human leukocyte antigen (HLA) genes driving diseases through their alleles or amino acid residues has been challenging. Using HLA information from HLA typing, imputation, or inference, HATK expands HLA alleles to amino acid sequences using the most recent IMGT/HLA database and prepares a dataset suitable for fine-mapping analysis. Our software also provides useful functionalities such as various association tests, visualization tools, and nomenclature conversion.

<br>
<br>


## (2) Installation

We strongly recommend using either 'Anaconda(or Miniconda)' or 'Docker' to set up HATK.

**(2-1) Docker**

1. Install Docker(https://docs.docker.com/)

    Note that different steps are followed to install Docker depending on your Operating System and sudo privilege is required.

2. Pull the pre-built HATK docker image file from Docker Hub. (https://hub.docker.com/r/wschoibhlab/hatk)
    ```
    docker pull wschoibhlab/hatk
    ```
<br>

**(2-2) Anaconda(or Miniconda)**


1. install Anaconda or Miniconda
    - Anaconda : (https://www.anaconda.com/)
    - Miniconda : (https://docs.conda.io/en/latest/miniconda.html)

2. Create a new independent work environment.

    By using 'HATK.yml' file in the project folder, Create a new work environment. 

    ```
    conda env create -f HATK.yml
    ```

    For more detailed explanation about Anaconda's managing environment, Please check this reference(https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually)

3. Download each dependent software in 'dependency/' folder.

    In the project folder, make a folder named 'dependency'. (with right exactly this name)
    ```
    mkdir dependency
    ```

    And download each dependent software to this directory.

    - beagle.jar
    - beagle4.jar
    - beagle2vcf.jar
    - beagle2linkage.jar
    - linkage2beagle.jar
    - vcf2beagle.jar
    - plink(v1.9)

    The copyright of 1 ~ 6 belongs to B. Browning (https://faculty.washington.edu/browning/beagle/b4_1.html) and 7 belongs to Purcell's laboratory (https://www.cog-genomics.org/plink/1.9/general_usage#cite)

<br>
<br>


## (3) Usage

**(3-1) Docker** 

run that pulled image (in the project directory.)


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

Note that, in case of using docker, output directory must be created inside of the project folder.

<br>

**(3-2) Anaconda**

run below command

```
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