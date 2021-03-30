#!/bin/bash

REF=../ref/Poecilia_reticulata.Guppy_female_1.0_MT.dna_sm.toplevel_withLambdaPUC19.fa

for FQZ1 in $(ls *_R1_*fastq.gz) ; do

  FQZ2=$(echo $FQZ1 | sed 's/_R1_/_R2_/' )
  FQ1=$(echo $FQZ1 | sed 's/fastq.gz/fastq-trimmed-pair1.fastq/')
  FQ2=$(echo $FQZ1 | sed 's/fastq.gz/fastq-trimmed-pair2.fastq/')
  BASE=$(echo $FQZ1 | cut -d '_' -f1)
  BAM=$BASE.bam
  VCF=$BAM.vcf
  VCFZ=$VCF.gz
  CGBED=$VCFZ.cg.bed
  CXBED=$VCFZ.cx.bed

  skewer -t 16 -q 20 $FQZ1 $FQZ2

  biscuit  align -@30 $REF $FQZ1 $FQZ2 \
  | samblaster \
  | samtools sort -o $BAM -O BAM -

  samtools index $BAM

  samtools flagstat $BAM > $BAM.stats.txt

  biscuit pileup -o $VCF $REF $BAM

  bgzip -f $VCF

  tabix -p vcf $VCFZ

  biscuit vcf2bed -t cg $VCFZ > $CGBED

  biscuit vcf2bed -t c $VCFZ > $CXBED

done
