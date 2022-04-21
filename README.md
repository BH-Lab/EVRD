# Iterative scrutiny pipeline for EVRD

The EVRD (Eukaryotic Viral Reference Database) is a collection of all viruses,
viral-related,and viral-like sequences from eukaryotes. The main method is to 
"purify" the GenBank and UniProt through host and vector sequence filtering 
and virus database cross validation. The fragments of host and vector sequences 
incorrectly annotated as viruses and the virus sequences incorrectly annotated 
across families were removed. This directory mainly includes downloading SRA 
data (sra.sh), establishing original virus database index (nT_db.smk, aa_db.smk)
and virus database processing process (snakfile-aa.smk, snakfile-nt.smk). The status of EVRD paper is under review.
(Chen J, Yan X, Sun Y, Ren Z, Yan G, Wang G, Liu Y, Zhao Z, Liu Y, Tu C, He B. 2022. De-heterogeneity of the eukaryotic viral 
reference database (EVRD) improves the accuracy and efficiency of viromic analysis. [bioRxiv link](https://doi.org/10.1101/2022.03.03.482774) )


## Users' Guide
This contains instructions to run the EVRD purification pipeline. 
This script consists of 4 parts in total                       
- Part I.   Host genomes iterative scrutiny                             
- Part II.  Vector sequence iterative scrutiny                         
- Part III. Annotation cross scrutiny                                 
- Part IV.  Viral metagenome cross scrutiny                             
Finally,you need to review the final documents produced by each part


## Installation
We recommend to use conda to build environment and run our script under it. 

1.Build the environment:

```bash
conda create -n EVRD_env -y
conda activate EVRD_env
conda install -c bioconda -c conda-forge blast diamond seqkit megahit snakemake sratools -y 
```
2.Download the code from our xxx
```bash
git clone https://github.com/BH-Lab/EVRD
```


## Usage
```bash

usage: snakemake -s snakefile-nt.smk  <options> [<args>]

Available options are:
    --printshellcmds, -p  Print out the shell commands that will be executed.
    --cores [N], --jobs [N], -j [N]
                        Use at most N cores in parallel (default: 1). If N is
                        omitted, the limit is set to the number of available
                        cores.
    -h, --help            show this help message and exit

```

There are 4 snakefiles used in our project, including:

1. snakefile-aa.smk & snakefile-nt.smk: The workflow of the purification process is described in snakefile-aa.smk and snakefile-nt.smk.
2. NT_db.smk & AA_db.smk: Nucleotide/Amino acid sequences was indexed by Blast/Diamond.
3. while.py: Merge the results of two files with the same ID into one file



If users want to test the code, the config files need to be modified:

- InputConfig.txt 
    - diamond_db: path to Amino acid virus database index in your workstation/PC/cluster
    - diamond_nr: path to nr index in your workstation/PC/cluster
    - diamond_NVPC: path to NVPC index in your workstation/PC/cluster
    - Blast_db: path to Nucleotide virus database index in your workstation/PC/cluster
    - Blast_nt: path to nt index in your workstation/PC/cluster

- ToolsConfig.txt
    - blast: path to blast
    - diamond: path to diamond in your workstation/PC/cluster
    - seqkit: path to seqkit in your workstation/PC/cluster


## Algorithm overview
<img src="D:\桌面\论文图表\Aditional File 2_00.png"/> 

**Overview of bioinformatic pipeline (left column) that is composed of** 
**heterogeneity scrutiny (upper), and finalization and assessment (lower), the right column indicates** 
**the sequence number in databases corresponding to each treatment. The genomic ID of hosts used** 
**in Host genomes scrutiny is shown in the bracket next to the species.**

##Questions/Comments
Any questions or comments should be addressed to [Dr. He](heb-001001@163.com) or Github issue.



