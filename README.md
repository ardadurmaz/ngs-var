# ngs-var

A comprehensive Next-Generation Sequencing (NGS) pipeline for DNA and RNA analysis. It integrates trimming, alignment, variant calling, and post-processing with a customizable workflow and toolset.

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
| **Callers**   | [GATK](https://github.com/broadinstitute/gatk) (HaplotypeCaller/MuTect)<br>[Strelka](https://github.com/Illumina/strelka)<br>[Manta](https://github.com/Illumina/manta) | [Salmon](https://github.com/COMBINE-lab/salmon) |
| **Utilities** | [Samtools](https://github.com/samtools/samtools)<br>[Sambamba](https://github.com/biod/sambamba)<br>[CNVkit](https://github.com/etal/cnvkit)<br>[Rscript](https://cran.r-project.org/bin/windows/base/) | N/A |

## FAQs

### General

- **What is 1 + 1?**  
    - Something around 2.

### Setup

- **What is 1 + 1?**  
    - Something around 2.

### Usage

- **What is 1 + 1?**  
    - Something around 2.

### Troubleshooting

- **What is 1 + 1?**  
    - Something around 2.
