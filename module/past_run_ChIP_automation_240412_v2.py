#!/usr/bin/env python3
#=============================================================================#
# Author	: DaeHee Kim
# Date		: 2024-04-12
# Usage		: run_ChIP_automation_240412_v2.py -f 0_filtering/test -c 30 -g mm
# Example	: 
# Description	: 
#=============================================================================#

from collections import defaultdict
import argparse, requests, os, glob, time, psutil

parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store', help='SRR included GEO file')
parser.add_argument('-c', action='store', help='core to use')
parser.add_argument('-g', action='store', help='genome (hs, mm, sc)')
args = parser.parse_args()
file = args.f
core = int(args.c)
genome = args.g
#=============================================================================#

# genome 정보에 맞는 인덱스 파일 위치 할당
if genome == "sc":
    bt_index = "/home/rlarl0240/Index/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Bowtie2Index/genome"
    genome = 12157105
elif genome == "hs":
    bt_index = "/home/rlarl0240/Index/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
elif genome == "mm":
    bt_index = "/home/rlarl0240/Index/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
elif genome == "dm":
    bt_index = "/home/rlarl0240/Index/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
    
#=============================================================================#

start_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(start_time)

if not os.path.isdir("1_prefetch"):
    os.mkdir("1_prefetch")
if not os.path.isdir("2_fastq"):
    os.mkdir("2_fastq")
if not os.path.isdir("3_trimmed"):
    os.mkdir("3_trimmed")
if not os.path.isdir("4_sam"):
    os.mkdir("4_sam")
if not os.path.isdir("5_bam"):
    os.mkdir("5_bam")
if not os.path.isdir("6_macs3"):
    os.mkdir("6_macs3")
    
pair_dict = dict()
processed_gsms = dict()

#=============================================================================#

# 출력 및 로그 작성용 함수
def log_writer(text, line=False):
    with open("log", "a") as f:
        if line:
            f.write("#"*80)
            f.write("\n")
        else:
            f.write("\n")
        f.write(text)
        f.write("\n")
        
#=============================================================================#
already_existed = set(f.split("_")[0] for f in sorted(os.listdir(f"6_macs3")))
bam_existed = set(f.split(".")[0] for f in sorted(os.listdir(f"5_bam")))
print(already_existed)
print(bam_existed)

treat_input_dict = dict()
with open(file, "r", errors='replace') as f:
    lines = f.readlines()
for line in lines:
    line = line.strip().split("\t")
    treat_gsm = line[2]
    input_gsm = line[27]
    pair = line[15]
    
    if treat_gsm in already_existed:
        continue
    elif treat_gsm in bam_existed:
        processed_gsms[treat_gsm] = input_gsm
        if pair == "PAIRED":
            pair_dict[treat_gsm] = True
        else:
            pair_dict[treat_gsm] = False
        continue

    treat_input_dict[treat_gsm] = input_gsm

#=============================================================================#

for treat_gsm, input_gsm in treat_input_dict.items():
    print(treat_gsm, input_gsm)
    treat_srrs, input_srrs = [], []
    
    treat_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?sp=runinfo&acc=" + treat_gsm
    try:
        response = requests.get(treat_url)
    except:
        log_writer(f"Gathering {treat_gsm} info failed", True)
        continue
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
    else:
        print("Failed to retrieve data: Status code", response.status_code)
        
    for line in lines[1:]:
        line = line.strip().split(",")
        srr_num = line[0]
        treat_srrs.append(srr_num)
        
    input_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?sp=runinfo&acc=" + input_gsm
    try:
        response = requests.get(input_url)
    except:
        log_writer(f"Gathering {input_gsm} info failed", True)
        continue
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        
    for line in lines[1:]:
        line = line.strip().split(",")
        srr_num = line[0]
        input_srrs.append(srr_num)
    
    #=============================================================================#
    
    for srr in treat_srrs+input_srrs:
        # prefetch
        log_writer(f"{srr} prefetch start", True)
        os.system(f"prefetch {srr} -O ./1_prefetch --max-size u >> log 2>> log")
        log_writer(f"{srr} prefetch end\n")

    #=============================================================================#
    
        #fasterq-dump
        fqd_path = f"./1_prefetch/{srr}"
        log_writer(f"{srr} fasterq-dump start", True)
        os.system(f"fasterq-dump {fqd_path} -O ./2_fastq >> log 2>> log")
        log_writer(f"{srr} fasterq-dump end\n")
    
    #=============================================================================#

        #trimming
        srr_files = sorted(glob.glob(f"2_fastq/{srr}*"))
        srr_list = []
        if len(srr_files) == 2:
            is_it_paired = "--paired "
            srr_list.append(srr_files[0])
            srr_list.append(srr_files[1])
        elif len(srr_files) == 1:
            is_it_paired = ""
            srr_list.append(srr_files[0])
        elif len(srr_files) == 3:
            is_it_paired = "--paired "
            srr_list.append(srr_files[1])
            srr_list.append(srr_files[2])
        
        trim_cmd = f"trim_galore {is_it_paired}--basename {srr} {' '.join(srr_list)} -o 3_trimmed -j 4 >> log 2>> log"
        
        log_writer(f"{srr} trimming start", True)
        os.system(trim_cmd)
        log_writer(f"{srr} trimming end\n")
        
#=============================================================================#



    # treat bowtie2
    print(treat_srrs)
    srr_1_list = []
    srr_2_list = []
    
    if len(treat_srrs) == 0:
        log_writer(f"{treat_gsm} {input_gsm} not processed", True)
        pair_dict[treat_gsm] = "skip"
        continue
    
    for srr in treat_srrs:
        srr_files = sorted(glob.glob(f"3_trimmed/{srr}*fq"))
        if len(srr_files) == 2:
            is_it_paired = True
            srr_1_list.append(srr_files[0])
            srr_2_list.append(srr_files[1])
        elif len(srr_files) == 1:
            is_it_paired = False
            srr_1_list.append(srr_files[0])
        
    if is_it_paired:
        bw_cmd = f"bowtie2 -q -p {core} -x {bt_index} -1 {' '.join(srr_1_list)} -2 {' '.join(srr_2_list)} -S ./4_sam/{treat_gsm}.sam >> log 2>> log"
        pair_dict[treat_gsm] = True
    else:
        bw_cmd = f"bowtie2 -q -p {core} -x {bt_index} -U {' '.join(srr_1_list)} -S ./4_sam/{treat_gsm}.sam >> log 2>> log"
        pair_dict[treat_gsm] = False
    
    log_writer(f"{treat_gsm} align start", True)
    log_writer(bw_cmd)
    os.system(bw_cmd)
    log_writer(f"{treat_gsm} align end\n")

    # input bowtie2
    print(input_srrs)
    srr_1_list = []
    srr_2_list = []
    for srr in input_srrs:
        srr_files = sorted(glob.glob(f"3_trimmed/{srr}*fq"))
        if len(srr_files) == 2:
            is_it_paired = True
            srr_1_list.append(srr_files[0])
            srr_2_list.append(srr_files[1])
        elif len(srr_files) == 1:
            is_it_paired = False
            srr_1_list.append(srr_files[0])
        
    if is_it_paired:
        bw_cmd = f"bowtie2 -q -p {core} -x {bt_index} -1 {' '.join(srr_1_list)} -2 {' '.join(srr_2_list)} -S ./4_sam/{input_gsm}.sam >> log 2>> log"
    else:
        bw_cmd = f"bowtie2 -q -p {core} -x {bt_index} -U {' '.join(srr_1_list)} -S ./4_sam/{input_gsm}.sam >> log 2>> log"
    
    log_writer(f"{input_gsm} align start", True)
    log_writer(bw_cmd)
    os.system(bw_cmd)
    log_writer(f"{input_gsm} align end\n")

#=============================================================================#

    os.system("rm 2_fastq/*")
    os.system("rm 3_trimmed/*")
    
#=============================================================================#

# convert sam to bam
sam_files = sorted(glob.glob(f"4_sam/*sam"))

while sam_files:
    processes = psutil.process_iter()
    samtools_running_count = 0
    for process in processes:
        try:
            name = process.name()
            user = process.username()
        except:
            pass
        else:
            if name == 'samtools' and user == "rlarl0240":
                samtools_running_count += 1
                
    #for i in range(int(core) - samtools_running_count):
    for i in range(8 - samtools_running_count):
        # 8개가 최대일듯
        if len(sam_files) == 0:
            break
        
        sam = sam_files.pop()
        sam = os.path.basename(sam).split(".")[0]
        print(sam)
        

        try:
            os.system(f"samtools view -bS 4_sam/{sam}.sam > 5_bam/{sam}.bam &")
        except:
            log_writer(f"{sam} converting fail")
            continue
        
        log_writer(f"{sam} converting start")
        

    print(f"Current converting sams : {samtools_running_count}")
    time.sleep(5)
    
while True:
    processes = psutil.process_iter()
    for process in processes:
        try:
            name = process.name()
            user = process.username()
        except:
            pass
        else:
            if name == 'samtools' and user == "rlarl0240":
                print("Waiting for converting")
                time.sleep(5)
                break
    else:
        print("Converting End!")
        break
    
#=============================================================================#

os.system("rm 4_sam/*")

#=============================================================================#

# macs3
IP_list = list(treat_input_dict.keys())

while IP_list:
    processes = psutil.process_iter()
    macs_running_count = 0
    for process in processes:
        try:
            name = process.name()
            user = process.username()
        except:
            pass
        else:
            if name == 'macs3' and user == "rlarl0240":
                macs_running_count += 1
                
    for i in range(int(core) - macs_running_count):
        
        if len(IP_list) == 0:
            break
        
        treat_gsm = IP_list.pop()
        input_gsm = treat_input_dict[treat_gsm]
        
        treat_gsm_file = f"5_bam/{treat_gsm}.bam"
        input_gsm_file = f"5_bam/{input_gsm}.bam"

        try:
            if pair_dict[treat_gsm] == True:
                bam_format = "BAMPE"
            elif pair_dict[treat_gsm] == 'skip':
                continue
            else:
                bam_format = "BAM"
        except:
            continue

        #macs3_cmd = f"macs3 callpeak -f {bam_format} -t {treat_gsm_file} -c {input_gsm_file} --broad -g {genome} --outdir 6_macs3 -n {treat_gsm}_{input_gsm} --broad-cutoff 0.1 -B --nomodel --extsize 147 >> log 2>> log &"
        macs3_cmd = f"macs3 callpeak -f {bam_format} -t {treat_gsm_file} -c {input_gsm_file} --broad -g {genome} --outdir 6_macs3 -n {treat_gsm}_{input_gsm} --broad-cutoff 0.1 -B >> log 2>> log &"

        try:
            os.system(macs3_cmd)
        except:
            log_writer(f"{treat_gsm}.{input_gsm} macs3 fail")
            continue
        
        log_writer(f"{treat_gsm}.{input_gsm} macs3 start")
        

    print(f"Current running macs : {macs_running_count}")
    time.sleep(5)
 
while True:
    processes = psutil.process_iter()
    for process in processes:
        try:
            name = process.name()
            user = process.username()
        except:
            pass
        else:
            if name == 'macs3' and user == "rlarl0240":
                print("Waiting for macs3")
                time.sleep(5)
                break
    else:
        print("Macs3 End!")
        break
    
#=============================================================================#
    
    # os.system("rm 2_fastq/*")
    # os.system("rm 3_trimmed/*")
    # os.system("rm 4_sam/*")
    # os.system("rm 5_bam/*")
    
#=============================================================================#

finish_time = time.strftime('%Y-%m-%d %H:%M:%S')
print(f"Start Time : {start_time}")
print(f"End Time : {finish_time}")

log_writer(f"Start Time : {start_time}", True)
log_writer(f"End Time : {finish_time}")

#=============================================================================#