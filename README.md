# HLA Analysis Toolkit (HATK)

## (1) Introduction

`HATK`(**HLA Analysis Tool-Kit**) is a collection of tools and modules to perform `HLA fine-mapping` analysis, which is to identify which HLA allele or amino acid position of the HLA gene is driving the disease. HLA fine-mapping analysis is an indispensable analysis in studies of autoimmune diseases.

In `GWAS`(Genome-wide Association Test) and its fine-mapping analysis, researchers can obtain candidate causal variants of the target disease. However, the association test performed on the variants in the HLA(Human Leukocyte Antigen) region, chromosome 6p21, usually shows unreliable results because this region bears an outlandish polymorphism. For example, it is common for the variants in this region to have tri-allelic variants or even more while typical variants in most of the genomic regions are usually bi-allelic. Also, it has been found that there are at least 7,000 alleles for only one HLA-B gene. Consequently, Performing conventional association test based on SNP array panel usually generate inaccurate association test result in the HLA region.

On the other hand, the `IPD-IMGT/HLA`, which is a specialist database, provides the official and most detailed information of the HLA region. Being updated roughly 4 times a year, They keep and manage whole HLA allele information. They name those many alleles based on the nomenclature defined by the '`WHO Nomenclature Committee For Factors of the HLA System`’. Also, they provide each HLA allele's (1) amino acid and (2) DNA sequence information. To use these useful data, Exact HLA allele information of patients is required and researchers may have to employ expensive HLA typing technologies. However, thanks to the recent development of HLA imputation and inference technologies, researchers now can obtain hundreds to thousands of patients’ HLA allele information and detour the cost issue of using HLA typing service.

Ultimately, HATK aims to perform an association test targeted to the HLA region. Based on patients’ HLA type information and its corresponding Amino acid and DNA sequence information distributed by the IMGT-HLA database, HATK builds a marker panel including not only the typical intergenic genomic variants(i.e. SNPs) markers but also variants of HLA region. Also, HATK provides the additional association test method so that researchers can analyze the signals arising in the amino acid sequence position.


<br>
<br>


## (2) Backgrounds

## (3) Installation

First of all, Download this project somewhere directory of your Linux(or OS_X) system.

```
$ git clone git@github.com:WansonChoi/HATK.git
$ cd HATK
```

We strongly recommend using 'Anaconda(or Miniconda)' to set up HATK. HATK supports OS X and Linux environments(ex. Ubuntu) and currently dosen't support Windows.



1. install Anaconda or Miniconda

    - Anaconda : (https://www.anaconda.com/)
    - Miniconda : (https://docs.conda.io/en/latest/miniconda.html)
<br>

2. Create a new independent Python virtual environment.

	By using 'HATK_LINUX.yml'('HATK_OSX.yml' in case of using OS X) file in the project folder, Create a new Python virtual environment.
    
	```
	conda env create -f HATK_LINUX.yml
	```
	
	The above command will generate a new Python virtual environment named 'HATK', which contains dependent Python packages, Java and R libraries, independent to your original Python system. For more detailed explanation about Anaconda's managing Python virtual environment, Please check this reference(https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually).

	If the new virtual environment is succuessfully installed, then activate it.

	```
	conda activate HATK
	```

	HATK will be implemented in this virtual environment. When you want to go back to your original Python system setting, then

	```
	conda activate base
	```
	or    
	```
	conda deactivate HATK
	```

	If you want to remove this virtual environment for HATK forever in your Anaconda, then

	```
	conda env remove -n HATK
	```
<br>
<br>

3. Download each dependent software in 'dependency/' folder.

    Though Anaconda handles installing necessary Python packages and R libraries with less endeavour, there are some other software which you must prepare manually by yourself (primarily due to copyright issue).

    In the project folder, make a folder named 'dependency'. (**with right exactly this name**)
    ```
    mkdir dependency
    ```

    And download below software to this directory and change their name as below.

    - beagle.jar (http://faculty.washington.edu/browning/beagle/b3.html - 'Old version'; **Choose the version "3.0.4"**)
    <!-- - beagle4.jar (https://faculty.washington.edu/browning/beagle/b4_1.html#download) -->
    - beagle2vcf.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - beagle2linkage.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - linkage2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - vcf2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html  - 'File conversion utilities')
    - plink - (https://www.cog-genomics.org/plink2 - v1.9beta)

    The copyright of 1 ~ 5 belongs to B. Browning (https://faculty.washington.edu/browning/beagle/b4_1.html) and that of 6 belongs to Purcell's laboratory (https://www.cog-genomics.org/plink/1.9/general_usage#cite).
<br>
<br>

<!-- 4. Install R language(https://www.r-project.org/) and below R packages.

    You **MUST** check wether these R packages are successfully installed or not. If not, Generating heatmap plot will fail.

    - gplots (https://cran.r-project.org/web/packages/gplots/index.html)
    - RColorBrewer (https://cran.r-project.org/web/packages/RColorBrewer/index.html)
    - shape (https://cran.r-project.org/web/packages/shape/index.html) -->


<br>
<br>


## (4) Usage example

```
python3 HATK.py \
    --variants example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --hped example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.hped \
    --2field \
    --pheno example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.phe \
    --pheno-name RA \
    --out MyHLAStudy/MyHLAStudy_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    -imgt 3320 \
    -hg 18 \
    --imgt-dir example/IMGTHLA3320 \
    --multiprocess 2
```

This command will implement (1) IMGT2Seq, (2) NomenCleaner, (3) bMarkerGenerator, (4) HLA_Analyzer(Association Test - logistic regression), (5) Manhattan Plot and (6) Heatmap Plot, which are the minimal components for HLA fine-mapping to be done.

On the other hand, as mentioned above, each module of HATK can be implemented repectively. **The README files of each of those modules are prepared in 'docs/' folder.** Those files include more detailed explanation and respective usage examples.

## (5) Citation

## (6) Lincense






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
