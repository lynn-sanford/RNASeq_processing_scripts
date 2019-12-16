These scripts are generally designed to work together

Run the scripts from the command line with the syntax included in the script headings

p05_gene_list_filter.sh comes first, running on DESeq2 result txt files to filter out all DE genes below a threshold padj
list_compare_foldchange.py comes next, if you want to compare two sets of significant genes to see which ones overlap
TPM_gene_list_extract.py can then be used on the comparison files to filter the TPM file to only the genes you care about
