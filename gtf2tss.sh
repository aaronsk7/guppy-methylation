#!/bin/bash

GTF=Poecilia_reticulata.Guppy_female_1.0_MT.103.gtf
G=Poecilia_reticulata.Guppy_female_1.0_MT.dna_sm.toplevel.fa.g
OUT=Poecilia_reticulata.Guppy_female_1.0_MT.103_tss.gtf

grep -w 'exon_number "1"' $GTF \
| awk '$7=="+"' \
| cut -f1,4,5,9 \
| cut -d '"' -f-2,6,12 \
| sed 's/ensembl/uncharacterised/' \
| awk '!arr[$2]++' \
| sed 's/gene_id "//' \
| sed 's/"/\t/g' \
| awk '{OFS="\t"} {print $1,$2,$2+1,"+",$4,$5,$6} ' > tmp_fwd.bed

grep -w 'exon_number "1"' $GTF \
| awk '$7=="-"' \
| cut -f1,4,5,9 \
| cut -d '"' -f-2,6,12 \
| sed 's/ensembl/uncharacterised/' \
| awk '!arr[$3]++' \
| sed 's/gene_id "//' \
| sed 's/"/\t/g' \
| awk '{OFS="\t"} {print $1,$3-1,$3,"-",$4,$5,$6} ' > tmp_rev.bed

cat tmp_fwd.bed tmp_rev.bed > tss.bed

bedtools slop -i tss.bed -g $G -b 1500 > tss_ext.bed

mv tss_ext.bed $OUT

rm tmp_fwd.bed tmp_rev.bed
