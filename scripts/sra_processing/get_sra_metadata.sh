#!/bin/bash
# Intended to be run by Nextflow, so files and dirs are created in the current dir
# $1: NCBI query

echo "Fetching metadata from NCBI esearch..."
esearch -db sra -query "$1" | efetch -format runinfo | cut -d ',' -f 1,7,8 | awk 'NR==1{print;next} {print| "sort -k1,1n -t\",\""}' > metadata_esearch.csv

echo "Fetching metadata from pysradb..."
pysradb metadata $(cat "metadata_esearch.csv" | cut -f 1 -d "," | tail -n +2) --detailed > metadata_pysradb.tsv

echo "Fetching taxonomy analysis metadata..."
mkdir -p "taxonomy_analysis/"
for run in $(tail -n +2 "metadata_pysradb.tsv" | cut -f 1 | less); do
    echo "Processing $run"
    wget -O "taxonomy_analysis/${run}.json" "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy?&acc=${run}&cluster_name=public"
    sleep 0.75
done