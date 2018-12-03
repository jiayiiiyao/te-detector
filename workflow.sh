#Step1: Align with repeat-masking(reference vs. full_reads)
##Prepare a reference genome

../soft/ncbi-blast-2.7.1+/bin/windowmasker -mk_counts -in hg38.fa > hg38.wmstat
../soft/ncbi-blast-2.7.1+/bin/windowmasker -ustat hg38.wmstat -outfmt fasta -in hg38.fa > hg38-wm.fa

../bin/lastdb -P8 -uNEAR -R11 -c mydb_wm hg38-wm.fa

##Substitution and gap rates

../bin/last-train -P8 -Q0 mydb_wm FAB45271.fa > FAB45271_wm.par

##Aligning DNA sequences

../bin/lastal -P8 -p FAB45271_wm.par mydb_wm FAB45271.fa | ../bin/last-split -m1 > FAB45271_wm.maf

#Step2: Find insertions and extract relevant reads
##Find insertions and log readnames
python tedet_insertion.py FAB45271_wm.maf > insertions.txt

## Extract reads
sh extract insertions.txt FAB45271.fa > FAB45271_redo.fa

#Step3: Realign carefully without repeat-masking (reference vs. extracted_reads)
../bin/lastdb -P8 -uNEAR -R11 -c mydb hg38.fa
../bin/last-train -P8 -Q0 mydb FAB45271_redo.fa > FAB45271_redo.par
../bin/lastal -P8 -p FAB45271_redo.par mydb FAB45271_redo.fa | ../bin/last-split -m1 > FAB45271_redo.maf

#Step4: Check candidate insertions
python tedet_anno.py FAB45271_redo.maf rmsk.txt > FAB45271_TE
