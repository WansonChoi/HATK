# HATK (HLA Analysis Toolkit)


## (1) Introduction


Human leukocyte antigen (HLA) genes encode the major histocompatibility complex (MHC) protein, which controls human immune responses and affects the susceptibility to various diseases. Identifying which allele or amino acid position of the HLA gene is driving the disease is called HLA fine-mapping, which is an indispensable analysis in studies of autoimmune diseases. However, for researchers who want to conduct HLA fine-mapping, it can be a burden because there are multiple technical problems that need to be solved such as acquiring HLA sequence information from IPD-IMGT/HLA(https://www.ebi.ac.uk/ipd/imgt/hla/) database, fitting HLA allele name in standard nomeclature system, preparing panel for association test, huge amount of text-preprocessing and etc. HATK provides a collection of tools that not only can solve those technical problems but also can help researchers to analyze the fine-mapping result.


## (1.5) Why do we need HATK? (Current Technical challenges)
#### [No standard file format for HLA information]

- HLA information을 담을 수 있는 파일 포맷이 standard가 없어서 되게 번거로움
 
There is no officially consented file format for patients' HLA allele information. Indeed, Various HLA related software use their own file format. Here, we introduce a new file format, HPED(HLA PED), which is similar to Plink ped file but consists of 22(6 + 8*2) columns. Left 6 columns are exactly the same as Plink ped file ('Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype') and other 16 columns are Individual's HLA diploid (unphased) genotypes(2 HLA alleles for each HLA gene) of 8 HLA genes(A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1; in this order). Those files in respective format can be converted to HPED format so that HLA allele information can be dealt with and stored in more uniform and systematic way.


#### [Too many HLA alleles and Evolving Nomenclature for them]

- HLA region은 특이해서 HLA gene에 들어갈 수 있는 allele들의 못해도 100개가 넘음.(인용. )
- 이렇게 많고 다양한 HLA allele들을 nominate할 수 있는 nomenclature system이 1개 이상있음. 
- 또 Major한 change가 있었음.(ex. 2-field(4-5 digits) without colons => Maximum 4-field with colons)

- 1년동안에도 두어번씩 업데이트되면서 계속 새로운 HLA allele들이 추가되기도 하고 기존의 allele들이 사라지기도 함.

Genes in HLA region are well known for their extreme polymorphism. For example, More than 6 thousands of alleles have been found in HLA-B gene so far based on the statistics of the IPD-IMGT/HLA(https://www.ebi.ac.uk/ipd/imgt/hla/stats.html). To handle those many alleles, nomenclature of 4-digit form without separator (ex. A*3301) was introduced for the first time in 1987. However, as more and more new HLA alleles were found, that nomenclature became not enough to embrace all alleles, which is called 'rollover' problem. So, around 2010, 'The WHO Nomenclature Committee for Factors of the HLA System' met and updated the old nomenclature(http://hla.alleles.org/nomenclature/naming_2010.html). The practical matter for researcher is that both old and updated nomeclature are used in HLA research area. Researches may confront a situation where HLA allele information is given in either the old or updated nomenclature. Also, they might need to convert a given HLA allele name within the old nomenclature to the one in the updated nomenclature, or vice versa.



#### [Continuously updating Database]

- HLA allele nomenclature과 더불어 Amino acid/DNA sequence들도 계속 조금씩 변함.
- 필요한 버전에 맞게 가져다 쓸 수 있게 하고 싶었음. (ex. 최신버전 반영)

As mentioned above, new HLA alleles are being found continuously. Some HLA alleles are excluded in the database because some alleles are found to be the same one. Some alleles are assigned extra fields. To include all those changes, the database is continously updated. Furthermore, the update happens in a short period relative to that of human genome build versions(ex. hg18, hg19, hg38). Researchers may want to choose the specific or latest version of IPD-IMGT/HLA database.


#### [Laborious work in data preparation for Association Test]

- Association Test를 하기 위해서는 Marker panel이 필요함. 필연적으로 redundant한 text processing이 수반되는데 이걸 HATK가 해결해줌. 이렇게 준비된 Marker panel로 가장 많이 활용되는 plink로 logistic regression으로 association test result를 구해줌. plink logistic regression이 아닌 statistical method를 활용하고 싶다면 해당 marker panel로 원하는데로 활용하면 됨.

Though researchers acquired HLA information such as sequences or alleles, they must prepare marker panel using those information so as to perform HLA fine-mapping, which is obviously laborious and redundant work. Furthermore, they need proper statistical methods for HLA fine-mapping and effective visualization tools to analyze the result of the test.

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

2. Create a new independent virtual Python environment.

	By using 'HATK_LINUX.yml'('HATK_OSX.yml' in case of using OS X) file in the project folder, Create a new virtual environment. 
    
	```
	conda env create -f HATK_LINUX.yml
	```
	
	The above command will generate a new virtual environment named 'HATK', which contains dependent Python packages, Java and R libraries, independent to your original Python system. For more detailed explanation about Anaconda's managing virtual environment, Please check this reference(https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually).

	If the new environment is succuessfully installed, then activate it.

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

	If you want to remove this virtual environment forever in your Anaconda, then

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
    - plink(**v1.9**) - (https://www.cog-genomics.org/plink2)

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


## (3) Usage example

```
python3 HATK.py \
    --variants example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --hped example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.hped \
    --2field
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
