# This script is easiest to run while in the directory with the files
# This runs on the p05 filtered DESeq2 results generated by the p05_gene_list_filter.sh script
# Load python/3.6.3 module before running

# Syntax: python list_compare_foldchange.py file1 file2 outfile

import sys

# Check if files are there

if len(sys.argv) < 3:
  print("Not enough files, rerun program")
  sys.exit()

# Take in filenames from command line

file1 = sys.argv[1]
file2 = sys.argv[2]
outfile = sys.argv[3]

# Initialize variables

data1 = []
data2 = []

#Import files and parse gene lists

i = 0

with open(file1) as f:
  for line in f:
    if i > 0:
      split_cols = line.split("\t")
      if len(split_cols) == 8:
        data1.append([split_cols[0], split_cols[2]])
        header1 = file1[14:23]
      else:
        split_cols = line.split("\n")
        comp_split = split_cols[0].split("\t")
        data1.append(comp_split)
    else:
      header1 = line.split('\n')
      header1 = header1[0].split("\t")
      if len(header1) > 1:
        header1 = header1[:]
        header1 = '\t'.join(header1)
        split_cols = line.split("\t")
        if len(split_cols) < 8:
          header1 = header1.split("ene\t")
          if len(header1) > 1:
            header1 = header1[1]
      else:
        header1 = []
    i = i + 1

i = 0

with open(file2) as f:
  for line in f:
    if i > 0:
      split_cols = line.split("\t")
      if len(split_cols) == 8:
        data2.append([split_cols[0], split_cols[2]])
        header2 = file2[14:23]
      else:
        split_cols = line.split("\n")
        comp_split = split_cols[0].split("\t")
        data2.append(comp_split)
    else:
      header2 = line.split('\n')
      header2 = header2[0].split("\t")
      if len(header2) > 1:
        header2 = header2[1:]
        header2 = '\t'.join(header2)
      else:
        header2 = []
    i = i + 1

# Determine iteration variables

index_iter = len(data1)
index_search = len(data2)
ref_list = list(zip(*data1))[0]
search_list = list(zip(*data2))[0]

# Find matches and record indeces of matches

matching = []

for i in range(index_iter):
  for j in range(index_search):
    if (search_list[j] == ref_list[i]):
      matching.append([i,j])

# Generate new list from matching

if len(header1) > 0:
  if len(header2) > 0:
    matching_list = ['\t'.join(['Gene', header1, header2])]
  else:
    matching_list = ['\t'.join(['Gene', header1])]
else:
  matching_list = ['Gene']

matching_list[0] = (matching_list[0] + '\n')

for i in range(len(matching)):
  append_line = []
  for j in range(len(data1[0])):
    append_line.append(data1[matching[i][0]][j])
  for k in range((len(data2[0]) - 1)):
    append_line.append(data2[matching[i][1]][k + 1])
  append_line = '\t'.join(append_line)
  append_line = (append_line + '\n')
  matching_list.append(append_line)

# Write out list

with open(outfile, 'wt') as f:
  for i in range(len(matching_list)):
    f.write(matching_list[i])


