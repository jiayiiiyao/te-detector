# te-detector
Detect transposable elements in long DNA reads.

Charaterize TE structures such as target-site duplication and 3'transduction.


#### Workflow
##### Step1: Install LAST
http://last.cbrc.jp/


##### Step2: Align reads to reference with repeat-masking(to reduce running time) using LAST

Check [here](https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md) for a detailed recipie

```
windowmasker -mk_counts -in reference.fa > reference.wmstat

windowmasker -ustat reference.wmstat -outfmt fasta -in reference.fa > reference-wm.fa

lastdb -P8 -uNEAR -R11 -c mydb-wm reference-wm.fa

last-train -P8 -Q0 mydb-wm reads.fa > reads-wm.par

lastal -P8 -p reads-wm.par mydb-wm reads.fa | ../bin/last-split -m1 > reads-wm.maf
```
 
 Here we use '-m1' option to keep all alignments: not omitting any alignment whose probability of having the wrong genomic locus is higher than some threshold
 
 
 ##### Step3: Collect insertion signals 
 
 ```
 python tedet_insertion.py [reads-wm.maf] [insertions]
 ```
 
 >Output format:
 
 >target-chr | target-site | readname | start-in-read | insertion-length | overlap/gap
 
 *overlap/gap: >0 gap; <0 overlap
 
 
 ##### Step4: Cluster insertion signals by coordinate and length
 
```
python cluster-insertion.py [insertions] [cluster-size] [out]
```
Log only clusters with size >= cluster-size (supporting reads)

>Output format:
 
>target-chr | target-site-cluster | average-length | target-site | readname | start-in-read | insertion-length | overlap/gap
 

#### Option: If you want to compare two datasets, and want to extract those specific insertions only in one dataset, you can use:

```
python diff-cluster.py [cluster1] [cluster2] [cluster1-cluster2]
```

The output includes insertion clusters in cluster1 but not in cluster2.



##### Step5: Identify transposable element insertions

```
python tedet_anno.py [insertion-cluster] [rmsk.txt] [out]
```


- **out_positive_full: (Full TE insertions are logged)**

> Row 1:

```
targetChr | targetSite | donorChr | donorStart | donorEnd | readname | readStart | readEnd | insertionStart | insertionEnd
```


> Row 2:

```
RepeatName | RepeatClass | RepeatFamily | RepeatStart | RepeatEnd
```


> Row 3:

```
 tsdLength | tsdSeq | leftFlankingSeq | rightFlankingSeq | 3'transductionLength | transductionSeq
```



- **out_positive_partial: (Partial TE insertions due to 5' truncation are logged)**

Same format as out_positive_full

