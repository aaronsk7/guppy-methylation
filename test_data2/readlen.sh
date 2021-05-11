#!/bin/bash

for FQ in *fastq ;
  do LEN=$(sed -n '2~4p' $FQ | awk '{print length($1)}' | numaverage)
  echo $FQ $LEN
done > readlen.txt

