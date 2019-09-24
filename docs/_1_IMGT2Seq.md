# IMGT2Seq

## (1) Introduction

The IPD-IMGT/HLA database contains major HLA information such as all HLA allele names, DNA or Amino acid sequences. IMGT2Seq preprocesses those information into new utilizable files(dictionary or table form) so that subsequenct steps of HATK can work on them. Consequently, IMGT2Seq generates 5 files. 

For example, 

1. HLA_DICTIONARY_AA.hg19.imgt3320.txt
2. HLA_DICTIONARY_AA.hg19.imgt3320.map
3. HLA_DICTIONARY_SNPS.hg19.imgt3320.txt
4. HLA_DICTIONARY_SNPS.hg19.imgt3320.map
5. HLA_INTEGRATED_ALLELE_TABLE.hg19.imgt3320.iat

1 is a dictionary of amino acid residue sequences ranging from the 1st exon to the last exon. 2 is genomic position information of those amino acid sequences. 3, 4 files are same as 1, 2 files but as to the DNA base pair sequences not only including exons but also introns between them. bMarkerGenerator will use those 1 - 4 files in preparing marker panel for associastion test.

5 is for NomenCleaner. All HLA allele names in the database are processed into 5 file. NomenCleaner works by transforming a given HLA allele name to the name within the updated nomenclature by searching this 5 file. 

When it comes to the genomic positions of 2 and 4 files, **Note that they are not perfectly compatible with those of Human genome build version(ex. hg18, hg19, hg38).** It's mainly because the IPD-IMGT/HLA database is continuously updated at least 1 or 2 times a year while Human Genome build is updated in a longer period than that. So, we made IMGT2Seq share only the start position of the 1st exon of each HLA gene based on IGV(https://software.broadinstitute.org/software/igv/). In other words, Consecutively incremented(or decremented) genomic positions are assigned to the rest of the positions except the start position. As a result, **Please be aware that it's highly undesirable to use the result of HATK in a research where genomic position information has to be strictly and thoroughly considered.**


## (2) Special characters

Usesrs might get confused due to some characters contained in the result of IMGT2Seq such as '*(asterisk)', '.(dot)', 'x' that look bizarre at first. To understand these characters, users should get the concept of 'Official Reference Sequence' and 'Virtual Sequence' defined by the IPD-IMGT/HLA. Official and detailed explanation about this can be found in the next link(https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html).

Briefly, there are the official reference seqeucnes for all HLA genes, e.g. A\*01:01:01:01 for HLA-A, B\*07:02:01:01 for HLA-B, which are used as a reference sequence. Meanwhile, for each HLA allele, all the individual sequence entries are submitted to the IPD-IMGT/HLA database and a virtual sequence is created based on them. This virtual sequence is then aligned to the official reference sequence. Insertions and deletions(Indels) coming out during this alignment are marked as '.(dot)'. 

Researchers would better distinguish between 'indels(dots) in the official reference sequence' and 'indels(dots) in the virtual sequence' because the former indels are not assigned numbering(See 'Numbering of the Sequence Alignment' section - https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html). In HATK, the spots of the former indels are processed to 'Z' or 'z', where 'Z' means insertions and 'z' does deletions. These spots are processed into markers in the form of 'INDEL_~' later in bMarkerGenerator. On the other hand, dots(indels) in the virtual sequence mean "Deletion in the virtual sequence".

'*(asterisk)' represents 'unknown at any point in the alignment'. In IMGT2Seq, this characters are processed to 'x' by us.




## (3) Usage example

```
python3 HATK.py \
    --imgt2seq \
    -hg 19 \
    -imgt 3320 \
    -o tests/_1_IMGT2Sequence/20190829_hg19_imgt3320/20190829.hg19.imgt3320 \
    --imgt-dir example/IMGTHLA3320 \
    --multiprocess 
```

A path to the folder of downloaded IPD-IMGT/HLA database must be given to '--imgt-dir' argument. Be aware which version of IMGT database('-imgt') and Human Genome version('-hg') you are going to use.