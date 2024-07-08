# LorBin

## Introduction

LorBin is a deep learning-based binner suitable for contigs assembled from long reads.

## System Requirements
LorBin is developed in Python 3.10 with modules and external tools.  
### Hardware requirements
LorBin requires only a standard computer with enough RAM to support the in-memory operations.  
### OS requirements
LorBin is supported and tested in Linux systems.
## Installation
### Install LorBin via source code
You can install LorBin from the source code. After installing Anaconda (or miniconda), first, obtain LorBin:  
git clone https://github.com/  
Then create an environment to run LorBin.  
```
cd path_to_LorBin
conda env create -f lorbin_env.yaml
conda activate comebin_env
```
## Preprocessing
The preprocessing steps aim to obtain a single FASTA file and generate bam files as input to our program.
### Concatenate the input contigs to a single FASTA file
If you have several FASTA files and want to bin them at once, you can use concat.py to concatenate the input contigs to a single FASTA file.
```
python concat.py -o path_to_all_contigs -fa path_to_FASTA_file1 path_to_FASTA_file2 ...
```
Your contig headers must be unique. Then, concat.py will use '-' to connect the sample index and contigs header.
### Generate bam files
```
minimap2 -a path_to_all_contigs path_to_FASTQ_file1 | samtools view -h -b -S | samtools view -b -F 4 | samtools sort -@ 20 > FILE1.mapped.sorted.bam
minimap2 -a path_to_all_contigs path_to_FASTQ_file1 | samtools view -h -b -S | samtools view -b -F 4 | samtools sort -@ 20 > FILE2.mapped.sorted.bam
...
```

## How to run LorBin
### Run LorBin via source code
#### Activate the environment
```
conda activate lorbin_env
cd path_to_LorBin
```
#### Binning
You can use subcommand 'bin' to bin the contigs.
```angular2html
python lorbin.py bin -o outputdir -fa path_to_all_contigs -b FILE1.mapped.sorted.bam FILE1.mapped.sorted.bam --multi
```
```angular2html
usage: lorbin.py bin [-h] -o OUTPUT -fa FASTA [--bin_length BIN_LENGTH] -b BAM [BAM ...] [--num_process NUM_PROCESS] [--evaluation EVALUATION] [-a AKEEP] [--multi]

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory (will be created if non-existent)
  -fa FASTA, --fasta FASTA
                        Path to the input fasta file.
  --bin_length BIN_LENGTH
                        minimum bin size in bps (Default: 80000)
  -b BAM [BAM ...], --bam BAM [BAM ...]
                        Path to the input BAM(.bam) file.
  --num_process NUM_PROCESS
                        Number of threads used (default: 10)
  --evaluation EVALUATION
                        Evaluation model used(no_markers, markers110, markers35, default: nomarkers
  -a AKEEP, --akeep AKEEP
                        The cut-off parameters of re-clustering decision model(0~1, default:0.6)
  --multi               Set True if uses multi samples
```
### Only generate data
If you only need the kmer and abundance data, you can use subcommand 'generate_data'.
```angular2html
python lorbin.py generate_data -o outputdir -fa path_to_all_contis -b FILE1.mapped.sorted.bam FILE1.mapped.sorted.bam
```
```angular2html
usage: lorbin.py generate_data [-h] -o OUTPUT -fa FASTA [--bin_length BIN_LENGTH] -b BAM [BAM ...] [--num_process NUM_PROCESS]

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory (will be created if non-existent)
  -fa FASTA, --fasta FASTA
                        Path to the input fasta file.
  --bin_length BIN_LENGTH
                        minimum bin size in bps (Default: 80000)
  -b BAM [BAM ...], --bam BAM [BAM ...]
                        Path to the input BAM(.bam) file.
  --num_process NUM_PROCESS
                        Number of threads used (default: 10)
```
### Only cluster used embedding.csv
If you have the embedded features of contigs and only want to use the two-stage clustering to bin, you can use subcommand 'cluster' to get the binning result. Note that the embedded features should be stored in .csv files and the index should be the contigs header.
```angular2html
python lorbin.py cluster -o outputdir -fa path_to_all_contigs --embedding path/embedding.csv
```
```angular2html
usage: lorbin.py cluster [-h] -o OUTPUT -fa FASTA [--bin_length BIN_LENGTH] [--evaluation EVALUATION] [-a AKEEP] [--multi] --embeddingdir EMBEDDINGDIR

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory (will be created if non-existent)
  -fa FASTA, --fasta FASTA
                        Path to the input fasta file.
  --bin_length BIN_LENGTH
                        minimum bin size in bps (Default: 80000)
  --evaluation EVALUATION
                        Evaluation model used(no_markers, markers110, markers35, default: nomarkers
  -a AKEEP, --akeep AKEEP
                        The cut-off parameters of re-clustering decision model(0~1, default:0.6)
  --multi               Cluster uses more samples
  --embeddingdir EMBEDDINGDIR, -e EMBEDDINGDIR
                        The path of embedding csv file used in clustering
```