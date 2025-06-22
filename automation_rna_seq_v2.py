#!/usr/bin/env python3
#=============================================================================#
# Author    : DaeHee Kim
# Date      : 2025-05-30
# Usage     : RNA-seq 자동화 개선버전
# Example   : 
# Description   : 
#=============================================================================#

"""

1. GSM을 SRR로 변환
2. SRR을 prefetch
3. SRA 파일 fasterq-dump
4. STAR로 align, fastq 파일 삭제, bam 파일 삭제
5. ReadsPerGene 파일 각 GSM 별로 이동
6. Deseq2 분석

신경 써야 할 점
- SRR이 없는 경우는 어떻게 할 것인가?
- - SRR이 없는 경우는 제외하고 진행

- SRR이 여러개인 경우는 어떻게 할 것인가?
- - SRR이 여러개인 경우는 합쳐서 진행

- Fastq 파일 single-end와 paired-end 구분

"""



import glob, random, os, time, argparse, requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store', help='GSM included GEO file')
parser.add_argument('-c', action='store', help='core to use', default= 32)
parser.add_argument('-g', action='store', help='genome (hs, mm, sc, dm, mg, xl)', default= "mm")
args = parser.parse_args()
file = args.f
core = int(args.c)
genome = args.g

if genome == "mm":
    star_index = "/program/STAR_index/mm10"
elif genome == "hs":
    star_index = "/program/STAR_index/hg38"
elif genome == "sc":
    star_index = "/program/STAR_index/r64"
elif genome == "dm":
    star_index = "/program/STAR_index/BDGP6"
elif genome == "mg":
    star_index = "/program/STAR_index/mg10"
elif genome == "xl":
    star_index = "/program/STAR_index/xl10"
    
if not os.path.isdir("1_prefetch"):
    os.mkdir("1_prefetch")
if not os.path.isdir("2_fastq"):
    os.mkdir("2_fastq")
if not os.path.isdir("3_fastq_combined"):
    os.mkdir("3_fastq_combined")
if not os.path.isdir("4_aligned"):
    os.mkdir("4_aligned")
if not os.path.isdir("5_ReadsPerGenes_GSM"):
    os.mkdir("5_ReadsPerGenes_GSM")
if not os.path.isdir("6_deseq2_results"):
    os.mkdir("6_deseq2_results")
if not os.path.isdir(f"./6_deseq2_results/merged"):
    os.mkdir(f"./6_deseq2_results/merged")
if not os.path.isdir(f"./6_deseq2_results/norm"):
    os.mkdir(f"./6_deseq2_results/norm")
if not os.path.isdir(f"./6_deseq2_results/deseq2"):
    os.mkdir(f"./6_deseq2_results/deseq2")

with open("log", "w") as f:
    pass 
    
#=============================================================================#

# 출력 및 로그 작성용 함수
def log_writer(text, line=False):
    with open("log", "a") as f:
        if line:
            f.write("#"*80)
            f.write("\n")
        else:
            f.write("\n")
        #print(text)
        f.write(text)
        f.write("\n")
        
#=============================================================================#

# 시작 파일은 GSE \t GSM,GSM \t GSM,GSM 으로 구성
start_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(start_time)

df = pd.read_excel(file,header=0, names = ['GSE', 'Treat','Control'])

grouped_data = defaultdict(list)
for name, group in df.groupby('GSE'):
    grouped_data[name] = group.to_dict(orient='records')

#=============================================================================#

# GSM을 SRR로 변환하는 함수
def gsm_to_srr(gsm):
    srr_arr = []
    treat_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?sp=runinfo&acc=" + gsm
    try:
        response = requests.get(treat_url)
    except:
        log_writer(f"Gathering {gsm} info failed", True)

    if response.status_code == 200:
        lines = response.text.strip().split("\n")
    else:
        print("Failed to retrieve data: Status code", response.status_code)
        
    for line in lines[1:]:
        line = line.strip().split(",")
        srr_num = line[0]
        srr_arr.append(srr_num)
    return srr_arr


#=============================================================================#

# 겹치는 GSE 구분
gse_counts = df['GSE'].value_counts()
gse_index = df.groupby('GSE').cumcount() + 1

def make_gse_unique(row):
    if gse_counts[row['GSE']] == 1:
        return row['GSE']
    else:
        return f"{row['GSE']}_{row['gse_index']}"

df['gse_index'] = gse_index
df['GSE_unique'] = df.apply(make_gse_unique, axis=1)
df = df.drop(columns=['gse_index'])

#=============================================================================#

# GSM to SRR 변환
gse_gsm_dict = df.groupby('GSE_unique').agg({
    'Treat': lambda x: x.iloc[0].split(','),
    'Control': lambda x: x.iloc[0].split(',')
    
}).to_dict(orient='index')

gsm_to_srr_i = 0
gsm_to_srr_dict = defaultdict(list)
for gse, cond_dict in gse_gsm_dict.items():
    gsms = cond_dict['Treat'] + cond_dict['Control']
    
    for gsm in gsms:
        gsm_to_srr_i += 1
        print(f"Converting {gsm_to_srr_i} GSM to SRR                   ", end="\r")
        srr = gsm_to_srr(gsm)
        if srr:
            gsm_to_srr_dict[gsm] += (srr)
        else:
            log_writer(f"Failed to convert {gsm} to SRR", True)


#=============================================================================#

def prefetch_srr(srr):
    os.system(f"prefetch {srr} -O ./1_prefetch --max-size u >> log 2>> log")
    time.sleep(1)
    
def fasterq_dump_srr(srr):
    os.system(f"fasterq-dump ./1_prefetch/{srr} -e 8 -O ./2_fastq >> log 2>> log")
    time.sleep(1)    

total_gse_count = len(gse_gsm_dict)
gse_i = 0
for gse, gsm_dict in gse_gsm_dict.items():
    treat_gsm = gsm_dict['Treat']
    control_gsm = gsm_dict['Control']
    gse_i += 1
    
    # SRR prefetch
    for gsm in treat_gsm + control_gsm:
        srr_list = gsm_to_srr_dict[gsm]
        for srr in srr_list:
            print(f"Prefetching {srr} for {gsm} ({gse} {gse_i} / {total_gse_count})                  ", end="\r")
            prefetch_srr(srr)
    
    # 현재 prefetch된 SRA 파일과 fastq 파일 확인
    current_prefecthed_arr = glob.glob(f"./1_prefetch/*/*.sra")
    current_prefecthed_arr = [os.path.basename(path)[:-4] for path in current_prefecthed_arr]
    current_fastq_arr = glob.glob(f"./2_fastq/*.fastq")
    current_fastq_arr = [os.path.basename(path).replace('.', '_').split("_")[0] for path in current_fastq_arr]
    
    # SRA fasterq-dump
    for gsm in treat_gsm + control_gsm:
        srr_list = gsm_to_srr_dict[gsm]
        for srr in srr_list:
            if srr in current_prefecthed_arr:
                if not srr in current_fastq_arr:
                    print(f"Fasterq-dump {srr} for {gsm} ({gse} {gse_i} / {total_gse_count})                  ", end="\r")
                    fasterq_dump_srr(srr)
                    
    # SRR 합치기
    current_combined_fastq_arr = glob.glob(f"./3_fastq_combined/*.fastq")
    current_combined_fastq_arr = [os.path.basename(path).replace('.', '_').split("_")[0] for path in current_combined_fastq_arr]

    
    for gsm in treat_gsm + control_gsm:
        if gsm in current_combined_fastq_arr:
            log_writer(f"{gsm} already combined fastq file exists", True)
            continue
        cat_command_1 = "cat"
        cat_command_2 = "cat"
        cat_command_single = "cat"
        srr_list = []
        for srrs in gsm_to_srr_dict[gsm]:
            
            srr_files = sorted(glob.glob(f"2_fastq/{srrs}*"))
            if len(srr_files) == 2:
                is_it_paired = True
                cat_command_1 += f" {srr_files[0]}"
                cat_command_2 += f" {srr_files[1]}"
                srr_list.append(srr_files[0])
                srr_list.append(srr_files[1])
            elif len(srr_files) == 1:
                is_it_paired = False
                cat_command_single += f" {srr_files[0]}"
                srr_list.append(srr_files[0])
            elif len(srr_files) == 3:
                is_it_paired = True
                cat_command_1 += f" {srr_files[1]}"
                cat_command_2 += f" {srr_files[2]}"
                srr_list.append(srr_files[1])
                srr_list.append(srr_files[2])
                log_writer(f"{gsm} : {srr} has paired and sigle" )
            else:
                log_writer(f"{gsm} : {srr} not found" )
            #print(cat_command_1, cat_command_2, cat_command_single)
            
        
        if is_it_paired:
            cat_command_1 += f" > ./3_fastq_combined/{gsm}_1.fastq"
            cat_command_2 += f" > ./3_fastq_combined/{gsm}_2.fastq"
            
            try:
                print(f"Combining fastq {gsm} ({gse} {gse_i} / {total_gse_count})                               ", end="\r")
                os.system(cat_command_1)
                os.system(cat_command_2)
                combining_success = True
            except:
                log_writer(f"{gsm} fastq combining failed")

        else:
            cat_command_single += f" > ./3_fastq_combined/{gsm}.fastq"
            
            try:
                print(f"Combining fastq {gsm} ({gse} {gse_i} / {total_gse_count})                               ", end="\r")
                os.system(cat_command_single)
                log_writer(f"{gsm} fastq combining end")
                combining_success = True
            except:
                log_writer(f"{gsm} fastq combining failed")
                
        if combining_success:
            for srr in srr_list:
                os.system(f"rm {srr}")

    # STAR로 align
    current_aligned_fastq_arr = glob.glob(f"./4_aligned/*ReadsPerGene.out.tab")
    current_aligned_fastq_arr = [os.path.basename(path).replace('.', '_').split("_")[0] for path in current_aligned_fastq_arr]
    
    for gsm in treat_gsm + control_gsm:
        if gsm in current_aligned_fastq_arr:
            log_writer(f"{gsm} already aligned", True)
            continue
        fastq_files = sorted(glob.glob(f"3_fastq_combined/{gsm}*"))
        if len(fastq_files) == 2:
            read_file_in = f"{fastq_files[0]} {fastq_files[1]}"
        else:
            read_file_in = f"{fastq_files[0]}"

        STAR_com = f"STAR --runThreadN {core} --quantMode GeneCounts --outSAMstrandField intronMotif --genomeDir {star_index} --readFilesIn {read_file_in} --outFileNamePrefix ./4_aligned/{gsm}_STAR_"
        
        try:
            print(f"STAR aligning {gsm} ({gse} {gse_i} / {total_gse_count})                               ", end="\r")
            os.system(STAR_com)
            align_success = True
        except:
            log_writer(f"{gsm} aligning failed")
            
        if align_success:
            for fastq in fastq_files:
                os.system(f"rm {fastq}")

    # 그룹별로 GSM 분리
    if not os.path.isdir(f"./5_ReadsPerGenes_GSM/{gse}"):
        os.mkdir(f"./5_ReadsPerGenes_GSM/{gse}")
    if not os.path.isdir(f"./5_ReadsPerGenes_GSM/{gse}/treat"):
        os.mkdir(f"./5_ReadsPerGenes_GSM/{gse}/treat")
    if not os.path.isdir(f"./5_ReadsPerGenes_GSM/{gse}/control"):
        os.mkdir(f"./5_ReadsPerGenes_GSM/{gse}/control")
        
    for gsm in treat_gsm:
        current_reads_per_gene = glob.glob(f"./4_aligned/{gsm}_STAR_ReadsPerGene.out.tab")
        if len(current_reads_per_gene) == 1:
            os.system(f"cp {current_reads_per_gene[0]} ./5_ReadsPerGenes_GSM/{gse}/treat/{gsm}_ReadsPerGene.out.tab")
        else:
            log_writer(f"{gsm} ReadsPerGene file not found", True)
    for gsm in control_gsm:
        current_reads_per_gene = glob.glob(f"./4_aligned/{gsm}_STAR_ReadsPerGene.out.tab")
        if len(current_reads_per_gene) == 1:
            os.system(f"cp {current_reads_per_gene[0]} ./5_ReadsPerGenes_GSM/{gse}/control/{gsm}_ReadsPerGene.out.tab")
        else:
            log_writer(f"{gsm} ReadsPerGene file not found", True)
            
    readpergenes_a = sorted(glob.glob(f"5_ReadsPerGenes_GSM/{gse}/treat/*_ReadsPerGene.out.tab"))
    readpergenes_b = sorted(glob.glob(f"5_ReadsPerGenes_GSM/{gse}/control/*_ReadsPerGene.out.tab"))

    genes_merged = defaultdict(list)
    meta = [["id", "group"]]
    sample_name = []

    for i, reads in enumerate(readpergenes_a):
        name = os.path.basename(reads).split("_ReadsPerGene")[0]
        sample_name.append(name)
        genes_merged["id"].append(f"A{i+1}")
        meta.append([f"A{i+1}","A"])
        with open(reads, "r") as f:
            lines = f.readlines()

        for line in lines[4:]:
            line = line.strip().split("\t")
            genes_merged[line[0]].append(line[1])
            
        print(f"{i+1} / {len(readpergenes_a)} {gse} treat merged",end="\r")
    print()

    for i, reads in enumerate(readpergenes_b):
        name = os.path.basename(reads).split("_ReadsPerGene")[0]
        sample_name.append(name)
        genes_merged["id"].append(f"B{i+1}")
        meta.append([f"B{i+1}","B"])
        with open(reads, "r") as f:
            lines = f.readlines()
            
        for line in lines[4:]:
            line = line.strip().split("\t")
            genes_merged[line[0]].append(line[1])
            
        print(f"{i+1} / {len(readpergenes_b)} {gse} control merged",end="\r")
    print()


    expression_file_name = f"6_deseq2_results/merged/{gse}_merged.txt"
    expression_meta_name = f"6_deseq2_results/merged/{gse}_merged.meta"

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
        f.write(f'write.table(data.frame("id"=rownames(normalized_counts),normalized_counts), file="6_deseq2_results/norm/{gse}_expression.deseq2_norm.txt", sep="\t", quote=F, col.names=TRUE, row.names=FALSE)\n')
        f.write('des <- DESeq(dds)\n')
        f.write('res <- results(des, contrast=c("group","A","B"))\n')
        f.write(f'write.table(data.frame("id"=rownames(res),res), file="6_deseq2_results/deseq2/{gse}_expression.deg.txt", sep="\t", quote=F, col.names=TRUE, row.names=FALSE)\n')
        
    os.system(f"Rscript run_DESeq2.R")
    
#=============================================================================#