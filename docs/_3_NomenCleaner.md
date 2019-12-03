# NomenCleaner


- 1987년에 최초로 HLA allele들의 Nomenclature가 정의된 이후로 계속 새로운 HLA allele들이 발견되면서 기존의 Nomenclature도 계속 변화를 겪음.
- 또한, 특별한 목적을 가지고 도입된 부가적인 Nomenclature들도 있음.(P-group, G-group).
- 따라서, 연구자들은 HLA allele정보들이 어떤 Nomenclature에 담겨져 받게 될지 모름.
- 또, 특정 Nomenclature에 담긴 HLA allele name들을 다른 Nomenclature 규칙을 따르는 이름으로 바꿔야할 필요성도 맞닥드리게 될 수 있음.
- Nomencleaner는 이에 대한 solution임.


## (1) WHO HLA Nomenclature

### (1-1) Old vs. Updated Nomenclature

Since the 'WHO Nomenclature Committee for Factors of the HLA System' defiend a nomenclature for HLA alleles for the first time in 1987, there have been several changes in the nomenclature. Especially, the committee did major update to the nomeclature in April 2010. (http://hla.alleles.org/nomenclature/naming_2010.html) 


In the old nomenclature, HLA allele name is supposed to have 

1. Only 2 digits in each field.
2. No field separator.

On the other hand, in the updated one, HLA allele name has

1. 2 to 3 digits in each field
2. Field separators
3. Maximum 4 fields.


(사진)
(http://hla.alleles.org/nomenclature/naming.html)


The IPD-IMGT/HLA database currently uses the updated one and researchers should use it. However, there are some practical reasons that hinder its usage.

- It can be burdensome for preceding researchers to move from the old to updated nomenclature.
- HLA typing technologies, which take a critical role in HLA research, still face a challenge in unambiguously determining full resolution of a single HLA allele, i.e. determinig all fields exactly. Consequently, though HLA allele names are based on the updated nomenclature, it is quite common to skip 3rd and 4th fields in order to aggregate multiple possible answers.

### (1-2) G and P codes
To deal with the above 2nd reason, the committee introduced G and P codes in the updated nomenclature. 

* G-group : HLA alleles that share identical nucleotide sequences for the exons encoding the peptide binding domains.
* P-group : HLA alleles that encode for identical peptide binding domains.

Because most HLA typing technologies focus on resolving alleles that encode differences within the peptide binding domains, those classification criteria can effectively represent an ambigous result of typing strategies.

(http://hla.alleles.org/nomenclature/naming_2010.html - '3. Reporting of ambiguous HLA allele typing')


In summary of '(1) WHO nomenclature', researchers should be able to recognize not only the old and updated nomenclatures but also additional G and P codes.

1. Old nomenclature (Without separator; defined in 1987 for the first time.)
2. Updated nomenclature (Maximum 4 fields with seperator; currently used by the IPD-IMGT/HLA)
3. G-group
4. P-group



## (2) What NomenCleaner does

NomenCleaner basically takes HLA alleles in either the updated, G-group or P-group nomenclature in HPED file('*.hped').

Then NomenCleaner converts their names to be in the new nomenclature which user requests (ex. 4-field, 2-field, G-group or P-group). In NomenCleaner, '4-field' allude to the updated HLA nomenclature which the IPD-IMGT/HLA currently uses.


=============(2019. 09. 17.)


## (3) Inevitable approximated mapping.

- 많은 사람들이 updated to old는 그렇다치지만 old to updated는 납득이 안갈거임. 이 부분에 관한 내용.
- 한 가지 명확한 사실은, 2-field(4-5 digits) 정보만 주어졌을 때, 3번째, 4번째 field에 대한 정보를 알아낼 수는 없음.
- 그러나, IPD-IMGT/HLA 데이터베이스 정보는 무조건 updated standard nomenclature로만 정보가 공개되기 때문에, sequence정보를 활용하기 위해서는 어떤식으로든 4-field형태의 HLA allele name이 준비되어야 함.
- 이 부분에 대해 NomenCleaner는 다음과 같은 approximated mapping을 수행함.

- 설명...
- 이렇게 approximated mapping된 4-field HLA allele name들을 기반으로 P-group과 G-group을 찾음( 혹은, 뒤에 suffix field들을 제거해서 다시 separator 구분이 되는 old 를 만듬.)
- 결론은, old to updated로 가는 NomenCleaner는 approximated된 결과물이란걸 꼭 염두해 놓고 후속 연구를 해라.





## (4) Usage Example.

NomenCleaner requires 'Integrated Allele Table(*.iat)' file from IMGT2Seq.


```
python3 HATK.py \
    --nomencleaner \
    -hat example/RESULT_EXAMPLE/HLA_ALLELE_TABLE.imgt3320.hat \
    --hped example/wtccc_filtered_58C_RA.hatk.58C_RA.300+300.hped \
    -o MyNomenCleaner/RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.58C_RA.300+300.chr6.hg18 \
    --2field \
    -imgt 3320

```