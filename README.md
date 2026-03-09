# rsvseq :high_brightness:

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

```bash
nextflow run main.nf -profile docker --runid runid_name   --input samplesheet.csv --outdir ../outdir_name
```

### Important Parameters

- `--input` (default: `assets/samplesheet.csv`): Path to the samplesheet.
- `--samplesDir` (default: `../data`): Directory containing the FASTQ files in the structure given above.
- `--primerdir` (default: `assets/primer/`): Legacy primer layout used when no primer scheme is set.
- `--primer_schemes_dir` (default: `assets/primer_schemes/`): Root directory for versioned primer schemes.
- `--primer_scheme` (default: `null`): Version folder name under `assets/primer_schemes/RSVA|RSVB/` (for example `V1`, `V2`).

All parameters are detailed in the `nextflow.config` file.

### Primer Scheme Selection

If `--primer_scheme` is not set, the pipeline uses:

- `assets/primer/a/reference.fasta`
- `assets/primer/a/primer.fasta`
- `assets/primer/b/reference.fasta`
- `assets/primer/b/primer.fasta`

If `--primer_scheme` is set, the pipeline uses:

- `assets/primer_schemes/RSVA/<SCHEME>/...`
- `assets/primer_schemes/RSVB/<SCHEME>/...`

Supported file names inside each scheme folder are:

- `reference.fasta` or `RSVA.reference.fasta` / `RSVB.reference.fasta`
- `primer.fasta` or `RSVA.primer.fasta` / `RSVB.primer.fasta`
- If FASTA is missing, the pipeline will generate it automatically from `primer.bed` / `RSVA.primer.bed` / `RSVB.primer.bed` plus the reference FASTA.

Example:

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --outdir ../outdir_name \
  --primer_scheme V1
```

## Pipeline Output

The output includes:

- Consensus sequences.
- Mutation calls.
- Sequencing statistics (coverage, quality parameters).
- A report in CSV format.
- A multiple FASTA file of sequences that passed quality filters.


## Credits

rsvseq was originally written by Rasmus Kopperud Riis.
