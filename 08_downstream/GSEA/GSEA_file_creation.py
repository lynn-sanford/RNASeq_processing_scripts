# Generally you shouldn't use filtered data for GSEA, so the gene_list file should be all the genes in the GTF file
# Make this file from the TPM file with the commands:
#  tail -n +2 all_samples_TPM.txt > all_samples_TPM_no_header.txt
#  awk '{print $1}' all_samples_TPM_no_header.txt > complete_gene_list.txt

# Syntax: python GSEA_file_creation.py gene_list TPM_file outfile

import sys

# Check if files are there

if len(sys.argv) < 3:
  print("Not enough files, rerun program")
  sys.exit()

# Take in filenames from command line

file1 = sys.argv[1]
file2 = sys.argv[2]
outfile = sys.argv[3]

data1 = []
data2 = []

# Read in gene_list 

with open(file1) as f:
  for line in f:
    split_cols = line.split('\n')
    data1.append(split_cols[0])

# Process TPM file

i = 0

with open(file2) as f:
  for line in f:
    if i > 0:
      split_cols = line.split("\n")
      split_cols = split_cols[0].split("\t")
      data2.append(split_cols)
    else:
      header = line.split('\n')
      header = header[0].split("\t")
    i = i + 1

# Find matches and record indeces of matches

TPM_gene_list = list(zip(*data2))[0]
gene_row = []

for i in range(len(data1)):
  gene_row.append(TPM_gene_list.index(data1[i]))

# Make GCT file

with open(outfile, 'wt') as f:
  f.write('#1.2\n')
  f.write(str(len(data1)))
  f.write('\t')
  f.write(str(len(header)))
  f.write('\n')
  f.write('Gene\tDescription\t')
  f.write('\t'.join(header[:]))
  f.write('\n')
  for i in range(len(gene_row)):
    f.write(TPM_gene_list[gene_row[i]].upper())
    f.write('\t')
    f.write(TPM_gene_list[gene_row[i]])
    f.write('\t')
    f.write('\t'.join(data2[gene_row[i]][1:]))
    f.write('\n')

