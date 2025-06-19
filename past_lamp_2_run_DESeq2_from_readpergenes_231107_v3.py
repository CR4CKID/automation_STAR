#!/usr/bin/env python3
#=============================================================================#
# Author	: DaeHee Kim
# Date		: 2023-11-07
# Usage		: python3 run_DESeq2_from_readpergenes_231107.py -t [treat] -c [control] -o [ouput_prefix]
# Example	: 
# Description	: 
#=============================================================================#

import argparse
import os
import glob
from collections import defaultdict
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store', help='treat')
parser.add_argument('-c', action='store', help='control')
parser.add_argument('-o', action='store', help='output_prefix')

args = parser.parse_args()
group_a = args.t
group_b = args.c
output_prefix = args.o

if group_a[-1] != "/":
    group_a += "/"
if group_b[-1] != "/":
    group_b += "/"

readpergenes_a = sorted(glob.glob(f"{group_a}*"))
readpergenes_b = sorted(glob.glob(f"{group_b}*"))

genes_merged = defaultdict(list)
meta = [["id", "group"]]
sample_name = []

for i, reads in enumerate(readpergenes_a):
    name = os.path.basename(reads).split("_STAR_ReadsPerGene")[0]
    sample_name.append(name)
    genes_merged["id"].append(f"A{i+1}")
    meta.append([f"A{i+1}","A"])
    with open(reads, "r") as f:
        lines = f.readlines()

    for line in lines[4:]:
        line = line.strip().split("\t")
        genes_merged[line[0]].append(line[1])
        
    print(f"{i+1} / {len(readpergenes_a)} {group_a[:-1]} merged",end="\r")
print()

for i, reads in enumerate(readpergenes_b):
    name = os.path.basename(reads).split("_STAR_ReadsPerGene")[0]
    sample_name.append(name)
    genes_merged["id"].append(f"B{i+1}")
    meta.append([f"B{i+1}","B"])
    with open(reads, "r") as f:
        lines = f.readlines()
        
    for line in lines[4:]:
        line = line.strip().split("\t")
        genes_merged[line[0]].append(line[1])
        
    print(f"{i+1} / {len(readpergenes_b)} {group_b[:-1]} merged",end="\r")
print()

a_name = group_a.split("/")[-2]
b_name = group_b.split("/")[-2]

expression_file_name = f"{a_name}_{b_name}_merged.txt"
expression_meta_name = f"{a_name}_{b_name}_merged.meta"

with open(expression_file_name, "w") as f:
    f.write("\t".join(["id"]+genes_merged["id"])+"\n")
    del(genes_merged["id"])
    
    for gene_name in sorted(genes_merged.keys()):
        f.write("\t".join([gene_name]+genes_merged[gene_name])+"\n")
        
with open(expression_meta_name, "w") as f:
    f.write("\t".join(meta[0])+"\n")
    
    for m in meta[1:]:
        f.write("\t".join(m)+"\n")
        
with open("run_DESeq2.R", "w") as f:
    f.write('library("DESeq2")\n')
    f.write(f'x <- read.delim("{expression_file_name}",header=TRUE,row.names="id",sep="\t")\n')
    f.write('roundx <- round(x)\n')
    f.write(f'coldata <-read.delim("{expression_meta_name}",header=TRUE,row.names="id",sep="\t")\n')
    f.write('dds <- DESeqDataSetFromMatrix(countData =roundx, colData = coldata, design =~group)\n')
    f.write('dds <- estimateSizeFactors(dds)\n')
    f.write('normalized_counts <- counts(dds, normalized=TRUE)\n')
    f.write(f'write.table(data.frame("id"=rownames(normalized_counts),normalized_counts), file="{output_prefix}_expression.deseq2_norm.txt", sep="\t", quote=F, col.names=TRUE, row.names=FALSE)\n')
    f.write('des <- DESeq(dds)\n')
    f.write('res <- results(des, contrast=c("group","A","B"))\n')
    f.write(f'write.table(data.frame("id"=rownames(res),res), file="{output_prefix}_expression.deg.txt", sep="\t", quote=F, col.names=TRUE, row.names=FALSE)\n')
    
os.system(f"Rscript run_DESeq2.R")

sample_name_line = "\t".join(sample_name)

# with open(expression_file_name, "r") as f:
#     lines = f.readlines()
# with open(f"named_{expression_file_name}", "w") as f:
#     f.write(f"id\t{sample_name_line}\n")
#     for line in lines[1:]:
#         f.write(line)
# with open(f"{output_prefix}_expression.deseq2_norm.txt", "r") as f:
#     lines = f.readlines()
# with open(f"named_{output_prefix}_expression.deseq2_norm.txt", "w") as f:
#     f.write(f"id\t{sample_name_line}\n")
#     for line in lines[1:]:
#         f.write(line)
data_1 = pd.read_csv(expression_file_name, sep="\t", names= sample_name, skiprows=1)
data_2 = pd.read_csv(f"{output_prefix}_expression.deseq2_norm.txt", sep="\t", names= sample_name, skiprows=1)

with pd.ExcelWriter(f'{output_prefix}_readcount_data.xlsx') as writer:
    data_1.to_excel(writer, sheet_name="raw")
    data_2.to_excel(writer, sheet_name="normalized")
    
    worksheet = writer.sheets['raw']
    # 열 너비 조정
    for column_cells in worksheet.columns:
        # 각 열의 최대 텍스트 길이 계산
        max_length = 0
        column = column_cells[0].column_letter  # 열 이름(A, B, C, ...)
        for cell in column_cells:
            try:
                # 셀의 값의 길이 계산
                if cell.value:
                    max_length = max(max_length, len(str(cell.value)))
            except:
                pass
        # 약간의 여유 공간을 주기 위해 2를 더함
        adjusted_width = max_length + 2
        worksheet.column_dimensions[column].width = adjusted_width
    
    worksheet = writer.sheets['normalized']
    for column_cells in worksheet.columns:
        # 각 열의 최대 텍스트 길이 계산
        max_length = 0
        column = column_cells[0].column_letter  # 열 이름(A, B, C, ...)
        for cell in column_cells:
            try:
                # 셀의 값의 길이 계산
                if cell.value:
                    max_length = max(max_length, len(str(cell.value)))
            except:
                pass
        # 약간의 여유 공간을 주기 위해 2를 더함
        adjusted_width = max_length + 2
        worksheet.column_dimensions[column].width = adjusted_width