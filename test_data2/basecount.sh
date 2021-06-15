#!/bin/bash

# The idea of thisscript is to calculate the number of bases removed at each QC step

# Starting bases

for FQ in *.fastq.gz  ; do
  CNT=$(zcat $FQ | sed -n '2~4p' | awk '{print length($1)}' | numsum )
  echo STARTING $FQ $CNT
done


# Bases kept after skewer

for FQ in *fastq ; do
  CNT=$(sed -n '2~4p' $FQ | awk '{print length($1)}' | numsum )
  echo AFTERSKEWER $FQ $CNT
done


# bases kept after removing unmapped and duplicate reads

for BAM in T*bam ; do
  CNT=$(samtools view -f 1 -F 1284 -q 40 $BAM \
  | cut -f10 | awk '{print length($1)}' | numsum )
  echo AFTERMAPPED $BAM $CNT
done

# Bases removed due to overlapping

for BAM in T*bam ; do
  CNT=$(samtools view -h -f 0x2 -F 1284 -q 40 T09-F-N.bam | samtools sort -n \
  | bamToBed -bedpe 2>/dev/null | awk 'NF==10' \
  | awk '$3>$5 {print $3 - $5}' | numsum )
  echo OVERLAPPING $BAM $CNT
done

