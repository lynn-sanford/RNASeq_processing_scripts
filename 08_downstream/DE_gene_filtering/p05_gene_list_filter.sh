#!/bin/bash

# Takes each DESeq2 result text file, puts the column headers in a new file, then filters all padj < 0.05 DE genes

for TXT in /Users/lysa8537/2018_RNASeq/accumulated_data/07_diff_ex/map2_glutrepfiltered/DESeq*.txt; do
head -n 1 $TXT > p05_gene_list_${TXT:89:9}.txt
tail -n +2 $TXT | sort -k7 | awk -v x=0.05 '$7<=x {print $0}' >> p05_gene_list_${TXT:89:9}.txt
done
