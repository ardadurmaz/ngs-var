# ngs-var

A comprehensive Next-Generation Sequencing (NGS) pipeline for DNA and RNA analysis, integrating trimming, alignment, variant calling, and post-processing with a customizable workflow and toolset.

## Installation 

    git clone https://github.com/ardadurmaz/ngs-var.git

## Setup

1. [Download](#components) the appropriate components for each process.
2. Specify paths in the [configuration file](ngs.config).
3. Enter sample details in the [run file](sample_run.tsv).

## Components

| **Category**  | **SNVIndel** | **RNASeq** |
|:-------------:|:------------:|:----------:|
| **Trimmers**  | [Fastp](https://github.com/OpenGene/fastp)<br>[Skewer](https://github.com/relipmoc/skewer)<br>[Cutadapt](https://github.com/marcelm/cutadapt)<br>[AfterQC](https://github.com/OpenGene/AfterQC) | [Fastp](https://github.com/OpenGene/fastp)<br>[Skewer](https://github.com/relipmoc/skewer)<br>[Cutadapt](https://github.com/marcelm/cutadapt)<br>[AfterQC](https://github.com/OpenGene/AfterQC) |
| **Aligners**  | [BWA](https://github.com/lh3/bwa) | [STAR](https://github.com/alexdobin/STAR) |
| **Callers**   | [GATK](https://github.com/broadinstitute/gatk)<br>[Strelka](https://github.com/Illumina/strelka)<br>[Manta](https://github.com/Illumina/manta) | [Salmon](https://github.com/COMBINE-lab/salmon) |
| **Utilities** | [Samtools](https://github.com/samtools/samtools)<br>[Sambamba](https://github.com/biod/sambamba)<br>[CNVkit](https://github.com/etal/cnvkit)<br>[Rscript](https://cran.r-project.org/bin/windows/base/) | - |

## Usage

### Arguments

- `--dir`: Name of working directory
- `--config`: Configuration file
- `--threads`: Number of threads to use
- `--title`: Title of the run
- `--log`: Name of log file
- `--verbose`: Verbosity
- `--dry`: Dry-Run
- `--clear`: Clear workspace
- `--upload`: Upload intermediate results to S3 storage (aws must be set for non-root user)
- `--exome`: Whole Exome Data
- `--bqsr`: Do BQSR
- `--trimmer`: Name of the trimmer (`Cutadapt`, `Fastp`, `Skewer`, `AfterQC`)
- `--aligner`: Name of the aligner (`BWA`, `STAR`)
- `--tool`: Name of the caller (`HaplotypeCaller`, `Strelka2`, `MuTect2`)
- `--workflow`: Type of analysis to run (`GermlineSNVIndel`, `SomaticSNVIndel`, `RNASeq`)
- `--in_file`: Targets file containing sample ids and associated fastq files (required)
- `--bed`: Bed file for regions
- `--cbed`: Bgzip compressed and tabix indexed bed file for regions

### Example
```bash
python3 -u ngs_main.py \
  --in_file /path/to/sample_run.tsv \
  --dir /path/to/ngs_WD \
  --config /path/to/ngs.config \
  --threads 16 \
  --title "Sample_Run" \
  --log /path/to/log.txt \
  --verbose \
  --dry \
  --clear \
  --upload \
  --exome \
  --bqsr \
  --trimmer fastp \
  --aligner BWA \
  --tool HaplotypeCaller \
  --workflow GermlineSNVIndel \
  --bed /path/to/hg38.bed \
  --cbed /path/to/capture.bed
```
