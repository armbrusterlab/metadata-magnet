# Metadata Magnet pipeline
<img width="1645" height="337" alt="Metadata Magnet long" src="https://github.com/user-attachments/assets/39370886-ab25-4f8e-a956-8a542dbd2484" />

## Table of Contents
- [Introduction](#introduction)
- [Pipeline tutorial](#pipeline-tutorial)
- [Pipeline overview](#pipeline-overview)
- [GUI screenshots](#gui-screenshots)
- [Acknowledgements](#acknowledgements)

## Introduction
This pipeline facilitates analysis of the ways in which bacterial protein alleles differ across environmental sources.

As of the time of writing, the [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) tool does not include metadata for bacterial proteins. Additionally, [NCBI calculated ortholog datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/) are currently available only for vertebrates and insects. However, by using this pipeline, you can find a set of orthologs for your bacterial protein of interest and collect metadata for these orthologs.

A variety of utilities are provided: you may filter your ortholog dataset in various ways, automatically categorize isolation source metadata into groups such as "host" or "natural environment", generate figures to visualize differences between categories, and run statistical tests on allele properties.

## Pipeline tutorial
To run the pipeline, please refer to the tutorial linked [here](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/pipeline_tutorial.md).

## Pipeline overview
Stages of the pipeline:
1. Obtain a BLAST database of bacterial protein sequences. You may download a prebuilt database such as nr, or use make_blastdb_from_genomes.sh to assemble one from a genome database downloaded using download_databases.sh or directly using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download).
2. In Nextflow: BLAST for potential orthologs of your bacterial protein sequence of interest, filter the dataset, retrieve the corresponding metadata, and optionally align the sequences. The main outputs of this stage are a **FASTA** and tab-separated **metadata table** for the ortholog dataset.
3. Optionally, further filter the dataset by filtering the metadata. This step may be necessary if you have specific filtering requirements not accounted for by the pipeline. We provide some helper functions in filter_metadata.py, but usually, filtering can be done in the command line.
4. In R Shiny GUI: Plug in the **FASTA** and **metadata table** to produce figures and run statistical tests.
<img width="750" height="983" alt="image" src="https://github.com/user-attachments/assets/e0f6993d-b71e-4b17-8bbb-d53fcc2c9ecb" />

## GUI screenshots
Data input tab:
<img width="2439" height="1351" alt="image" src="https://github.com/user-attachments/assets/f048d26e-cbe4-4e0c-9d0f-10ac6803e79c" />

Figure generation tab:
<img width="3000" height="1747" alt="image" src="https://github.com/user-attachments/assets/6f426a2b-4eb9-4d33-8355-b1c697bb6df1" />

Statistical test tab:
<img width="3000" height="1747" alt="image" src="https://github.com/user-attachments/assets/f0fc154c-7af6-4d6c-8be1-ebb796b3c058" />

## Acknowledgements
* Advisor: Dr. Catherine Armbruster
* Thesis committee:
   * Dr. Irene Kaplow
   * Dr. Phillip Compeau
* Assistance with bioinformatic analysis: Dr. Arkadiy Garber
* Open-source projects used for this pipeline:
   * [bit](https://github.com/AstrobioMike/bit)
   * [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
   * [pynteny](https://github.com/Robaina/Pynteny)
   * [seqtk](https://github.com/lh3/seqtk)
