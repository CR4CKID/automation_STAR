# automation_STAR : the automated RNA-seq pipeline for identify differentially expressed genes (DEGs)

<br>

## Prerequisites

To run this project, make sure you have the following Python or R packages installed (or newer):

| Package  | Version (>=) | Description                                                                 |
|----------|--------------|-----------------------------------------------------------------------------|
| `numpy`  | >=2.0.0       | Core library for numerical computing and array operations                  |
| `pandas` | >=2.2.2       | Powerful data structures for data analysis and manipulation                 |
| `prefetch` | >=3.1.1       |	NCBI SRA Toolkit utility for downloading sequencing data from the Sequence Read Archive (SRA) |
| `fasterq-dump`   | >=3.1.1      |SRA Toolkit utility for converting SRA files to fastq format |
| `STAR`  | >=2.7.10a      | Ultrafast universal RNA-seq aligner for mapping sequencing reads to a reference genome |
| `DESeq2`  | >=1.34.0      | R package for differential gene expression analysis based on the negative binomial model |


## Usage
| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `-f`     | ✔        | –       | Directory containing DESeq2 result files |
| `-g`     | ✔        | `mm` | genome (hs, mm, sc, dm, mg, xl)) |
| `-c`     | ✘        | `32`    | Number of CPU cores |

<br>
*Datasets must be prepared in the format of DESeq2 output results <br><br>

Simple example :
<pre lang="markdown"> hStouffer.py -d dataset_directory -o meta_dataset_fat </pre>

More detail example :
<pre lang="markdown"> hStouffer.py -d dataset_directory -o meta_dataset_max -r 10000 -c 64 -l max </pre>

<br>


## Contact
Daehee Kim : rlarl0240@knu.ac.kr <br>
Jun-yeong Lee : junyeong@knu.ac.kr <br>
Laboratory of Genome Architecture and Regulation, School of Life Sciences, College of Natural Sciences, Kyungpook National University
