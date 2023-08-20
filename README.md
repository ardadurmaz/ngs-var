# ngs-var
A comprehensive Next-Generation Sequencing (NGS) pipeline for DNA and RNA analysis, integrating trimming, alignment, variant calling, and post-processing with a customizable workflow and toolset.

# Installation 
 ```sh
 git clone https://github.com/ardadurmaz/ngs-var.git
```
# Setup
- Download the appropriate <a href="#Components">components</a> for each process
- Specify paths in the <a href="ngs.config">configuration file</a>
- Enter sample details in the <a href="sample_run.tsv"> run file</a>

# Components
| Category   | SNVIndel                                                                                              | RNASeq                                                        |
|------------|-------------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| **Trimmers**  | [Fastp](https://github.com/OpenGene/fastp), [Skewer](https://github.com/relipmoc/skewer), [Cutadapt](https://github.com/marcelm/cutadapt), [AfterQC](https://github.com/OpenGene/AfterQC)                               | - |
| **Aligners**  | [BWA](https://github.com/lh3/bwa)                                                                                           | [STAR](https://github.com/alexdobin/STAR)                                                  |
| **Callers**   | [GATK](https://github.com/broadinstitute/gatk) (HaplotypeCaller/MuTect), [Strelka](https://github.com/Illumina/strelka), [Manta](https://github.com/Illumina/manta)                              | [Salmon](https://github.com/COMBINE-lab/salmon)                                                |
| **Utilities** | [Samtools](https://github.com/samtools/samtools), [Sambamba](https://github.com/biod/sambamba), [CNVkit](https://github.com/etal/cnvkit), [Rscript](https://cran.r-project.org/bin/windows/base/)                                 | N/A                                                           |


# FAQs
## General
### What is 1 + 1?
Something around 2.
## Setup
### What is 1 + 1?
Something around 2.
## Usage
### What is 1 + 1?
Something around 2.
## Troubleshooting
### What is 1 + 1?
Something around 2.
