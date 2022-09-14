#!/bin/bash
\grep 'hg38::chr' hg38_fair+new_CAGE_peaks_phase1and2.bed | awk -F '\t' 'BEGIN {OFS = FS} {gsub("chr","",$1)} {print $1, $2+1, $3, $4, $5, $6, $8}' > ONLY_hg38_fair+new_CAGE_peaks_phase1and2.bed
bedtools intersect -a Homo_sapiens.GRCh38.95.utr_table.dedup.bed -b ONLY_hg38_fair+new_CAGE_peaks_phase1and2.bed -wo -s | sort -k1,1 -k2,2n -V > Homo_sapiens.GRCh38.95.utr_table.cages.tmp2
grep '5_utr' Homo_sapiens.GRCh38.95.utr_table.cages.tmp2 > Homo_sapiens.GRCh38.95.utr_table.utr5.cages.tmp2
./mysql_cages_table_loader.py Homo_sapiens.GRCh38.95.utr_table.cages.tmp2 
mv b.gz > Homo_sapiens.GRCh38.95.cages.gz
