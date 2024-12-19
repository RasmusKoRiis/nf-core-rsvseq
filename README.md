# sarsseq :high_brightness:

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

This pipeline processes FASTQ files from Nanopore sequencing of RSV A and B, generating consensus sequences and analyzing them for mutations and sequencing statistics The main steps include:

- Alignment and consensus sequencing with IRMA.
- Consensus sequence analysis with Nextclade.
- Mutation calling
- Generation of a comprehensive report in CSV format.
- Output of sequences in a multiple FASTA file.

The pipeline consist of four different worflows listed bellow:

1) RSV A and B FASTQ analysis (human)
  Alignment of FASTQ and mutation analysis 
2) RSV A and B FASTA analysis (human-fasta) (under development)
   Mutation analysis 

## Compatibility

- **Operating System**: Linux
- **Dependencies**: Docker and Nextflow

## Usage

### Sample Sheet Preparation

Prepare a sample sheet (CSV or TSV*) in the `assets` folder with the following format:
* TSV file is not not compulsory

```
PCR-PlatePosition,SequenceID,Barcode,KonsCt
A1*,sampleID,barcodeID,ct-value*
```
*not compulsory

Each row lists a sample to be analyzed. Samples not listed in the sheet will be excluded from the analysis.

### Directory Structure

#### For FASTQ-analysis
Ensure your directory structure is as follows:

```
./
  |-data
         |-barcode3
               |-XXXX_pass_barcode03_XXXX.fastq.gz
               |-YYYY_pass_barcode03_YYYY.fastq.gz
  |-nf-core-rsv
               |-assets
                     |-samplesheet.csv
                     |-samplesheet.tsv
               |-...
```

### Running the Pipeline

Navigate to the `nf-core-rsv` folder and execute the following command with default parameters:

#### SARS-CoV-2 FASTQ analysis

```bash
nextflow run main.nf -profile docker --runid runid_name   --input samplesheet.csv --outdir ../outdir_name
```

### Important Parameters

- `--input` (default: `assets/samplesheet.csv`): Path to the samplesheet.
- `--samplesDir` (default: `../data`): Directory containing the FASTQ files in the structure given above.

All parameters are detailed in the `nextflow.config` file.

## Pipeline Output

The output includes:

- Consensus sequences.
- Mutation calls.
- Sequencing statistics (coverage, quality parameters).
- A report in CSV format.
- A multiple FASTA file of sequences that passed quality filters.


## Credits

rsvseq was originally written by Rasmus Kopperud Riis.

