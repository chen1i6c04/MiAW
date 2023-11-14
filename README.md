# MiAW
MiSeq Analysis Workflow of CDC central lab

## Installation
```
conda env create -f requirement.yaml
```
## Usage

### Basic
```
python miaw.py -1 <R1.fastq> -2 <R2.fastq> -o <output_dir>
```

### Classify with Kraken2
```
python miaw.py -1 <R1.fastq> -2 <R2.fastq> -o <output_dir> --kraken2_db <database>
```
* Can find pre-build database in [here](https://benlangmead.github.io/aws-indexes/k2)

### Check Genome Completeness
```
python miaw.py -1 <R1.fastq> -2 <R2.fastq> -o <output_dir> --busco_db <database>
```