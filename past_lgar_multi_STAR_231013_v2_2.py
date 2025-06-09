#!/usr/bin/env python3
#=============================================================================#
# Author	: DaeHee Kim
# Date		: 2023-10-13
# Usage		: python3 multi_STAR_231013.py -f [SRR 포함 GSE 파일] -c [코어 수] -g [생물]
# Example	:  python3 multi_STAR_231013.py -f young_SRR.txt -c 50 -g mm
# Description	: 
#=============================================================================#

import argparse
from collections import defaultdict
import os
import glob
import time

parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store', help='SRR included GEO file')
parser.add_argument('-c', action='store', help='core to use')
parser.add_argument('-g', action='store', help='genome (hs, mm, sc)')
args = parser.parse_args()
file = args.f
core = int(args.c)
genome = args.g

if genome == "mm":
    star_index = "/program/STAR_index/mm10"
elif genome == "hs":
    star_index = "/program/STAR_index/hg38"
elif genome == "sc":
    star_index = "/program/STAR_index/Saccharomyces_cerevisiae" 
elif genome == "dm":
    star_index = "/program/STAR_index/dm6"

    
GSEs = defaultdict(list)
with open(file, 'r', errors='replace') as f:
    lines = f.readlines()

for line in lines:
    line = line.strip().split("\t")
    GSEs[line[2]].append(line[-6])

srr_list = []
for srrs in GSEs.values():
    srr_list += srrs

#=============================================================================#

start_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(start_time)

srr_nums = len(srr_list)
gsm_nums = len(GSEs.keys())

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

for cur_gsm_index, gsm in enumerate(GSEs.keys()):
    
    finished_files = sorted(glob.glob(f"4_aligned/{gsm}*"))
    if len(finished_files) > 4:
        print(f"{gsm} already existed")
        continue
    
    
# prefetch
    if not os.path.isdir("1_prefetch"):
        os.mkdir("1_prefetch")

    for i, srr in enumerate(GSEs[gsm]):
        log_writer(f"{srr} prefetch start", True)
        print(f"{i+1} / {len(GSEs[gsm])}\tCurrent prefetch GSM : {gsm}\t{cur_gsm_index+1} / {gsm_nums}", end="\r")
        os.system(f"prefetch {srr} -O ./1_prefetch --max-size u >> log 2>> log")
        log_writer(f"{srr} prefetch end\n")
    print(f"{gsm} Prefetch Done                                          ")
        
    #=============================================================================#

    # fasterq-dump 실행
    if not os.path.isdir("2_fastq"):
        os.mkdir("2_fastq")

    for i, srr in enumerate(GSEs[gsm]):
        fqd_path = f"./1_prefetch/{srr}"
        log_writer(f"{srr} fasterq-dump start", True)
        print(f"{i+1} / {len(GSEs[gsm])}\tCurrent fasterq GSM : {gsm}\t{cur_gsm_index+1} / {gsm_nums}", end="\r")
        os.system(f"fasterq-dump {fqd_path} -e 16 -O ./2_fastq >> log 2>> log")
        log_writer(f"{srr} fasterq-dump end\n")
    print(f"{gsm} Fasterq-dump Done                                          ")

    #=============================================================================#

    # SRR 합치기

    if not os.path.isdir("3_fastq_combined"):
        os.mkdir("3_fastq_combined")

    cat_command_1 = "cat"
    cat_command_2 = "cat"
    cat_command_single = "cat"

    srr_list = []
    
    for srrs in GSEs[gsm]:
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
            
    if is_it_paired:
        cat_command_1 += f" > ./3_fastq_combined/{gsm}_1.fastq"
        cat_command_2 += f" > ./3_fastq_combined/{gsm}_2.fastq"
        
        log_writer(f"{gsm} fastq combining start", True)
        try:
            print(f"Current fastq combining : {gsm}\t{cur_gsm_index+1} / {gsm_nums}", end="\r")
            os.system(cat_command_1)
            os.system(cat_command_2)
            log_writer(f"{gsm} fastq combining end")
            combining_success = True
        except:
            log_writer(f"{gsm} fastq combining failed")

    else:
        cat_command_single += f" > ./3_fastq_combined/{gsm}.fastq"
        
        log_writer(f"{gsm} fastq combining start", True)
        try:
            print(f"Current fastq combining : {gsm}\t{cur_gsm_index+1} / {gsm_nums}", end="\r")
            os.system(cat_command_single)
            log_writer(f"{gsm} fastq combining end")
            combining_success = True
        except:
            log_writer(f"{gsm} fastq combining failed")
            
    if combining_success:
        for srr in srr_list:
            os.system(f"rm {srr}")
                
    print(f"{gsm} Fastq Combining Done                                          ")

    #=============================================================================#

    # Align
    if not os.path.isdir("4_aligned"):
        os.mkdir("4_aligned")

    fastq_files = sorted(glob.glob(f"3_fastq_combined/{gsm}*"))
    if len(fastq_files) == 2:
        read_file_in = f"{fastq_files[0]} {fastq_files[1]}"
    else:
        read_file_in = f"{fastq_files[0]}"

    STAR_com = f"STAR --runThreadN {core} --quantMode GeneCounts --outSAMstrandField intronMotif --genomeDir {star_index} --readFilesIn {read_file_in} --outFileNamePrefix ./4_aligned/{gsm}_STAR_"
    
    log_writer(f"{gsm} aligning start", True)
    try:
        print(f"Current fastq combining : {gsm}\t{cur_gsm_index+1} / {gsm_nums}")
        os.system(STAR_com)
        log_writer(f"{gsm} aligning end")
        align_success = True
    except:
        log_writer(f"{gsm} aligning failed")
        
    if align_success:
        for fastq in fastq_files:
            os.system(f"rm {fastq}")

    print(f"{gsm} Align Done")

#=============================================================================#

if not os.path.isdir("5_ReadsPerGenes"):
    os.mkdir("5_ReadsPerGenes")

os.system("cp 4_aligned/*ReadsPerGene.out.tab 5_ReadsPerGenes/")

#=============================================================================#

finish_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(f"Start Time : {start_time}")
print(f"End Time : {finish_time}")

log_writer(f"Start Time : {start_time}", True)
log_writer(f"End Time : {finish_time}")

#=============================================================================#