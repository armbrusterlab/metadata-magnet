# Metadata Magnet pipeline
<img width="1645" height="337" alt="Metadata Magnet long" src="https://github.com/user-attachments/assets/39370886-ab25-4f8e-a956-8a542dbd2484" />

## Table of Contents
- [Introduction](#introduction)
- [Pipeline tutorial](#pipeline-tutorial)
- [Pipeline overview](#pipeline-overview)
- [GUI screenshots](#gui-screenshots)
- [Acknowledgements](#acknowledgements)

## Introduction
This pipeline facilitates comparative analysis of bacterial proteins across environmental sources. Using this pipeline, you can retrieve a set of homologs for your bacterial protein of interest and collect the corresponding metadata. Isolation source metadata is automatically categorized according to an ontology, greatly reducing the time required to curate data and to bin it for statistical analysis.

A variety of utilities are provided for pipeline outputs: you may filter your homolog dataset in various ways, automatically categorize isolation source metadata into groups such as "clinical" or "natural environment", generate figures to visualize differences between categories, and run statistical tests on allele properties.

## Pipeline tutorial
To run the pipeline, please refer to the tutorial linked [here](https://github.com/armbrusterlab/metadata-magnet/blob/main/pipeline_tutorial.md).  
Although this pipeline is designed to obtain homologs of your target bacterial protein sequence, there is a [workaround](https://github.com/armbrusterlab/metadata-magnet/blob/main/run_with_orthologs.md) if you are specifically interested in orthologs of your protein. This document doubles as a guide for running a few of the pipeline scripts as stand-alones.  

## Pipeline overview
Stages of the pipeline:
1. Obtain a BLAST database of bacterial protein sequences. You may download a prebuilt database such as nr, or use make_blastdb_from_genomes.sh to assemble one from a genome database downloaded using download_databases.sh or directly using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download).
2. In Nextflow: BLAST for potential homologs of your bacterial protein sequence of interest, filter the dataset, retrieve the corresponding metadata, and optionally align the sequences. The main outputs of this stage are a **FASTA** and tab-separated **metadata table** for the homolog dataset.
3. Optionally, further filter the dataset by filtering the metadata. This step may be necessary if you have specific filtering requirements not accounted for by the pipeline. We provide some helper functions in filter_metadata.py, but usually, filtering can be done in the command line.
4. In R Shiny GUI: Plug in the **FASTA** and **metadata table** to produce figures and run statistical tests.
<img width="2205" height="3461" alt="image" src="https://github.com/user-attachments/assets/f8a89180-6ac3-4f16-92ce-4dbb92baffa7" />


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
