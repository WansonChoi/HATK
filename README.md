# HATK


## (1) What is HATK?

HATK is ...(under construction)

1. What is HATK
2. What output does it generates?
3. What input is required to generate the output.
<br>
<br>


## (2) Installation

We strongly recommend using 'Anaconda(or Miniconda)' to set up HATK. HATK supports OS X and Linux environment(ex. Ubuntu) and currently dosen't support Windows.



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

    And download below software to this directory.

    - beagle.jar(http://faculty.washington.edu/browning/beagle/b3.html - 'Old version'; **Choose the version "3.0.4"**)
    <!-- - beagle4.jar (https://faculty.washington.edu/browning/beagle/b4_1.html#download) -->
    - beagle2vcf.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - beagle2linkage.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - linkage2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - vcf2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - plink(v1.9) - (https://www.cog-genomics.org/plink2)

    The copyright of 1 ~ 5 belongs to B. Browning (https://faculty.washington.edu/browning/beagle/b4_1.html) and 6 belongs to Purcell's laboratory (https://www.cog-genomics.org/plink/1.9/general_usage#cite).

4. Install R language(https://www.r-project.org/) and below R packages.

    - gplots
    - RColorBrewer
    - shape


<br>
<br>


## (3) Usage example

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