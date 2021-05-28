#!/bin/bash

insertlen(){
  BAM=$1
  ( samtools sort -n $BAM \
  | bedtools bamtobed -bedpe -i  \
  | awk '{print $6-$2}' > $BAM.insertlen.txt )  2>/dev/null
}
export -f insertlen

parallel -j 16 insertlen ::: *bam
