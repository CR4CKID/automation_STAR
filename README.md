# automation_STAR : the automated RNA-seq pipeline for identify differentially expressed genes (DEGs)

<br>

## Prerequisites

To run this project, make sure you have the following Python or R packages installed (or newer):

| Package  | Version (>=) | Description                                                                 |
|----------|--------------|-----------------------------------------------------------------------------|
| `numpy`  | >=2.0.0       | Core library for numerical computing and array operations                  |
| `pandas` | >=2.2.2       | Powerful data structures for data analysis and manipulation                 |
| `prefetch` | >=3.1.1       |             |
| `fasterq-dump`   | >=3.1.1      |                            |
| `STAR`  | >=2.7.10a      |    |

<br>
Prerequisities can be simply install using requirements.txt

<pre lang="markdown"> pip install -r requirements.txt </pre>
<br>

## Usage
| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `-d`     | ✔        | –       | Directory containing DESeq2 result files |
| `-p`     | ✘        | `False` | p-value threshold (e.g., 3 = 1e-3) |
| `-o`     | ✔        | –       | Output file name |
| `-m`     | ✘        | `False` | Group size |
| `-r`     | ✘        | `1000`  | Number of repetitions |
| `-n`     | ✘        | `False` | Apply scaling (`True` / `False`) |
| `-t`     | ✘        | `"stouffer"` | Method for p-value combination:<br>`fisher`, `pearson`, `tippett`, `stouffer`, `mudholkar_george`, `median`, `percentile_70` |
| `-c`     | ✘        | `32`    | Number of CPU cores |
| `-w`     | ✘        | `True`  | Apply p-value capping (`True` / `False`) |
| `-l`     | ✘        | `False` | Max or min cutoff threshold |
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
