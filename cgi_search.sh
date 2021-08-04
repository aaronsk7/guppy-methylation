#!/bin/bash

SEQ=$1
OUT=$(echo $SEQ | sed 's/.fa$/_out.txt/')
GFF=$(echo $SEQ | sed 's/.fa$/_out.gff/')

cpgplot -sequence $SEQ -window 100 -minlen 200 -minoe 0.6 -minpc 50 -outfile $OUT -outfeat $GFF -noplot
