# te-detector
Detect transposable elements in long DNA reads

#### Workflow
##### Step1: Install LAST
http://last.cbrc.jp/


##### Step2: Align reads to reference with repeat-masking(to reduce running time) using LAST

Check [here](https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md) for a detailed recipie

`windowmasker -mk_counts -in reference.fa > reference.wmstat`

`windowmasker -ustat reference.wmstat -outfmt fasta -in reference.fa > reference-wm.fa`

`lastdb -P8 -uNEAR -R11 -c mydb-wm reference-wm.fa`

`last-train -P8 -Q0 mydb-wm reads.fa > reads-wm.par`

`lastal -P8 -p reads-wm.par mydb-wm reads.fa | ../bin/last-split -m1 > reads-wm.maf`
 
 Here we use '-m1' option to keep all alignments: not omitting any alignment whose probability of having the wrong genomic locus is higher than one threshold.
 
 ##### Step3: Detect insertions and extract relevant reads
 
 `python tedet_insertion.py reads-wm.maf readList`
 
 readList is a text file with names of reads with insertions.
 
 `sh extract.sh readList reads.fa > reads_redo.fa`
 
 ##### Step4: Re-align reads_redo.fa carefully without repeat-masking to reference
 
`lastdb -P8 -uNEAR -R11 -c mydb reference.fa`

`last-train -P8 -Q0 mydb reads_redo.fa > reads_redo.par`

`lastal -P8 -p reads_redo.par mydb reads_redo.fa | last-split -m1 > reads_redo.maf`

##### Step5: Re-detect insertions and check whether they can be aligned to annotated TEs

`python tedet_anno.py reads_redo.maf rmsk.txt out`

Here you can get three output file with prefix 'out' (modifiy as you want)

out_insertion:

Log all detected insertions in such format:

targetChr | targetSite | readname | insertionStart | insertionEnd | tsd_length | tsd_seq | leftFlanking | rightFlanking

out_positive_full

out_positive_partial:

