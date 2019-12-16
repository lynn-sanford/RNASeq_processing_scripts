# This script is designed to work on a union gene list created by the list_compare_foldchange.py script to filter TPM values
# Syntax: python TPM_gene_list_extract.py gene_list TPM_file outfile

import sys
import statistics
from statistics import mean

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

# Read in gene_list data

i = 0

with open(file1) as f:
  for line in f:
    if i > 0:
      split_cols = line.split('\n')
      split_cols = split_cols[0].split("\t")
      data1.append(split_cols)
    else:
      header_list = line.split('\n')
      header_list = header_list[0].split("\t")
    i = i + 1

# Process TPM file

i = 0

with open(file2) as f:
  for line in f:
    if i > 0:
      split_cols = line.split("\n")
      split_cols = split_cols[0].split("\t")
      data2.append(split_cols)
    else:
      header_tpm = line.split('\n')
      header_tpm = header_tpm[0].split("\t")
    i = i + 1

# Find matches and record indeces of matches

TPM_gene_list = list(zip(*data2))[0]
gene_row = []

for i in range(len(data1)):
  gene_row.append(TPM_gene_list.index(data1[i][0]))


# Calculate average TPMs and up/down reg and store in new structure

avgs = []

for i in range(len(gene_row)):
  cont_avg = mean(list(map(float,data2[gene_row[i]][1:6])))
  glut_avg = mean(list(map(float,data2[gene_row[i]][6:11])))
  glzn_avg = mean(list(map(float,data2[gene_row[i]][11:16])))
  gtpa_avg = mean(list(map(float,data2[gene_row[i]][16:21])))
  if glzn_avg > gtpa_avg:
    updown = 'up'
  else:
    updown = 'down'
  avg_obj = list(map(str,[cont_avg,glut_avg,glzn_avg,gtpa_avg,updown]))
  avgs.append(avg_obj)

# Make header line

header = ['Gene', 'Cont_avg', 'Glut_avg', 'GlZn_avg', 'GTPA_avg','Up/Down','']
for i in range(len(header_list)-1):
  header.append(header_list[i+1])

header.append('')
for i in header_tpm:
  header.append(i)

# Make readable list of genes, TPMs, average TPMs, up or down reg, p-values for comparisons

with open(outfile, 'wt') as f:
  f.write('\t'.join(header))
  f.write('\n')
  for i in range(len(gene_row)):
    f.write(TPM_gene_list[gene_row[i]].upper())
    f.write('\t')
    f.write('\t'.join(avgs[i]))
    f.write('\t\t')
    f.write('\t'.join(data1[i][1:]))
    f.write('\t\t')
    f.write('\t'.join(data2[gene_row[i]][1:]))
    f.write('\n')

