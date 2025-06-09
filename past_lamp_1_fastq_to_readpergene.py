#!/usr/bin/env python3
#=============================================================================#
# Author    : DaeHee Kim
# Date      : 2024-08-01
# Usage     : 램프 성의 표시용 summary
# Example   : lamp_auto_DEG_2.py -d 5_DEG/ -g da_th.gmt -s 2 -l 0
# Description   : 메타분석후 파일을 바로 사용
#=============================================================================#

import glob, random, os
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import math

import argparse
from collections import defaultdict
import os
import glob
import time

#=============================================================================#

parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store', help='SRR included GEO file')
parser.add_argument('-c', action='store', help='core to use')
parser.add_argument('-g', action='store', help='genome (hs, mm, sc, dm, mg, xl)')
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

if not os.path.isdir("1_gz"):
    os.mkdir("1_gz")
if not os.path.isdir("2_fastq"):
    os.mkdir("2_fastq")
if not os.path.isdir("3_aligned"):
    os.mkdir("3_aligned")
if not os.path.isdir("4_ReadsPerGenes"):
    os.mkdir("4_ReadsPerGenes")
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

start_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(start_time)

df = pd.read_excel(file,header=None, names = ['name', 'group', 'fastq_1','fastq_2'])

grouped_data = defaultdict(list)
for name, group in df.groupby('name'):
    grouped_data[name] = group.to_dict(orient='records')
    
#=============================================================================#

# 다운로드 및 압축 해제
already_existed = set(f.split(".")[0] for f in sorted(os.listdir(f"2_fastq")))
already_stared = set(f.split("_STAR_")[0] for f in sorted(os.listdir(f"3_aligned")))

total_key_len = len(grouped_data.keys())
i = 1
for key, value in grouped_data.items():
    print(f"\n\n** {i} / {total_key_len} processing\n\n")
    data = value[0] # name, group, fastq_1, fastq_2
    fastq_1_file_name = data['fastq_1'].split("/")[-1].split(".")[0]
    fastq_2_file_name = data['fastq_2'].split("/")[-1].split(".")[0]

    if data['fastq_1'] != "X":
        os.system(f"wget -nc -P 1_gz {data['fastq_1']}")
        os.system(f"wget -nc -P 1_gz {data['fastq_2']}")

    if f"{key}_1" not in already_existed:
        os.system(f"cp 1_gz/{key}_1* 2_fastq")
        os.system("gzip -f -d 2_fastq/*gz")
    if f"{key}_2" not in already_existed:
        os.system(f"cp 1_gz/{key}_2* 2_fastq")
        os.system("gzip -f -d 2_fastq/*gz") 

    if key in already_stared:
        print(f"{key} already aligned")
    else:
        STAR_com = f"STAR --runThreadN {core} --quantMode GeneCounts --outSAMstrandField intronMotif --genomeDir {star_index} --readFilesIn 2_fastq/{key}_1.fastq 2_fastq/{key}_2.fastq --outFileNamePrefix ./3_aligned/{key}_STAR_ >> log 2>> log"
        os.system(STAR_com)
        
    group_name = value[0]['group']
    if not os.path.isdir(f"4_ReadsPerGenes/{group_name}"):
        os.mkdir(f"4_ReadsPerGenes/{group_name}")
    os.system(f"cp 3_aligned/{key}_STAR_ReadsPerGene.out.tab 4_ReadsPerGenes/{group_name}")
    
    i += 1
    
finish_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(f"Start Time : {start_time}")
print(f"End Time : {finish_time}")