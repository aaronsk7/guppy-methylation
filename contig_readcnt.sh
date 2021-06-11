#!/bin/bash

# this script count number of uniquely mapped reads on each main chromosome
# and spike in contig like pUC19 and lambda

for BAM in *bam ; do
  for CONTIG in $(samtools idxstats $BAM | cut -f1 | tac | sed 1d | grep -v "KK" ) ; do
    CNT=$(samtools view -c -q40 $BAM $CONTIG )
    echo $BAM $CONTIG $CNT
  done
done > contig_readcnt.txt
