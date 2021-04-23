#!/bin/bash

# calculate the read length of fastq files after read trimming with skewer

for FQ in *fastq ; do

  AVLEN=$(sed -n '2~4p' $FQ | awk '{print length($1)}' \
  | numaverage)

  echo $FQ $AVLEN

done | tee readlen.txt
