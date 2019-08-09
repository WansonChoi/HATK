# HATK (HLA Analysis Toolkit)


## (1) Introduction


Human leukocyte antigen (HLA) genes encode the major histocompatibility complex (MHC) protein, which controls human immune responses and affects the susceptibility to various diseases. Identifying which allele or amino acid position of the HLA gene is driving the disease is called HLA fine-mapping, which is an indispensable analysis in studies of autoimmune diseases. However, for researchers who want to conduct HLA fine-mapping, it can be a burden because there are multiple technical problems that need to be solved such as acquiring HLA sequence information from the IPD-IMGT/HLA official database, fitting HLA allele name in standard nomeclature system, preparing panel for association test, huge amount of text-preprocessing and etc. HATK provides a collection of tools that not only can solve those technical problems but also can help researchers to analyze the fine-mapping result.

* HATK란? : HLA region에서 association test를 수행할 수 있게 해주고, 이를 통해 얻은 association signal을 fine-mapping할 수 있게 도와주는 tool.
* HLA Finemapping을 하는데 있어 맞닥드리게 될 문제점들.
* 그리고 HATK가 제공하는 이에 대한 Solution들.

## (1.5) Technical chanllenges
#### [Preparing HLA information]
HATK works primarily based on 'Plink' file format. Here, we defined one more file format, HPED(HLA PED) file format, which is similar to Plink ped file but consists of 6 + 8 columns. Left 6 columns are exactly same as Plink ped file ('Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype') and other 8 columns are Individual's HLA diploid genotype(2 HLA alleles for each individual) over 8 HLA genes(A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1).

#### [Too many HLA alleles]



#### [Continuously updating Database]


#### [Generating panel for Association Test]


<!-- HATK is ...(under construction)

1. What is HATK
2. What output does it generates?
3. What input is required to generate the output. -->
<br>
<br>


## (2) Installation

We strongly recommend using 'Anaconda(or Miniconda)' to set up HATK. HATK supports OS X and Linux environment(ex. Ubuntu) and currently dosen't support Windows.



1. install Anaconda or Miniconda
    - Anaconda : (https://www.anaconda.com/)
    - Miniconda : (https://docs.conda.io/en/latest/miniconda.html)
<br>
<br>
2. Create a new independent work environment.

    By using 'HATK.yml' file in the project folder, Create a new virtual environment. 

    ```
    conda env create -f HATK.yml
    ```

    This command will generate a new work environment named 'HATK', which handles dependent Python packages or libraries of HATK, independent to your system. For more detailed explanation about Anaconda's managing environment, Please check this reference(https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually)

    After a new environment succuessfully installed, activte it.
    ```
    conda activate HATK
    ```

    HATK can be implemented in this environment. When you want to go back to your system setting, then
    ```
    conda activate base
    ```
    or
    ```
    conda deactivate HATK
    ```

    If you want to remove this virtual environment forever in your system, then

    ```
    conda remove -n HATK
    ```
    <br>

3. Download each dependent software in 'dependency/' folder.

    Though Anaconda can facilitate installing necessary Python packages or libraries with less endeavour, there are some other software which you must install manually by yourself.

    In the project folder, make a folder named 'dependency'. (**with right exactly this name**)
    ```
    mkdir dependency
    ```

    And download below software to this directory.

    - beagle.jar (http://faculty.washington.edu/browning/beagle/b3.html - 'Old version'; **Choose the version "3.0.4"**)
    <!-- - beagle4.jar (https://faculty.washington.edu/browning/beagle/b4_1.html#download) -->
    - beagle2vcf.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - beagle2linkage.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - linkage2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - vcf2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - plink(**v1.9**) - (https://www.cog-genomics.org/plink2)

    The copyright of 1 ~ 5 belongs to B. Browning (https://faculty.washington.edu/browning/beagle/b4_1.html) and that of 6 belongs to Purcell's laboratory (https://www.cog-genomics.org/plink/1.9/general_usage#cite).
<br>
<br>

4. Install R language(https://www.r-project.org/) and below R packages.

    - gplots (https://cran.r-project.org/web/packages/gplots/index.html)
    - RColorBrewer (https://cran.r-project.org/web/packages/RColorBrewer/index.html)
    - shape (https://cran.r-project.org/web/packages/shape/index.html)


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