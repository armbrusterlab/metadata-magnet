#!/bin/bash

# informal script for the sole purpose of filtering my wspF_wspR synteny searches to emulate wspE_wspF_wspR searches because Pynteny wasn't cooperating on that

summary=$1 # synteny summary file
sign=$2 # "pos" or "neg"

outname="${summary%.*}_wspE_filtered.tsv" # just save it to the same dir as the input file...
echo "genome_id	contig	organism	isolation_source	titles	locus_tag	protein_id	sequence	seq_tech	assembly_method	genome_coverage	isolation_site" > $outname

log=$(dirname $summary)/wspE_filtering.log
echo "match_count	genome_id	contig	wspF_locus_tag	wspF_locus_num	wspE_locus_num	strand" > $log

# I think the hmmer outputs are identical regardless of sign, but just to be sure
#synteny_outputs="/home/kcw2/metadata-magnet/nextflow/work/7d/cff95a6665b91271e3a92093fc08dc" # work dir from last successful Nextflow run
# wait, we have to search in a directory for which we DID search for wspE, otherwise HMMer wouldn't have looked for wspE at all
synteny_outputs="/home/kcw2/data/blast_outputs/wspF/synteny_search_outputs_updated" # the issue is that I've changed the set of genomes searched since then
if [[ $sign = "pos" ]]; then
  datadir="${synteny_outputs}/posForward/genomes"
else
  datadir="${synteny_outputs}/negBackward/genomes"
fi

echo "Searching ${datadir} for wspE hits"

IFS=$'\n' # rows in synteny summary file are newline separated
for line in $(cat "$summary" | tail -n +2); do 
  # get info from synteny summary file
  genome_id=$(echo $line | cut -f 1)
  contig=$(echo $line | cut -f 2)
  locus_tag=$(echo $line | cut -f 6)

  
  # get the locus number for wspE, and adjust according to $sign
  wspF_locus=$(cat "${datadir}/${genome_id}/hmmer_outputs/hmmer_output_PIRSF000876.txt" | grep "$locus_tag" | head -n 1 | cut -d " " -f 1 | cut -d "_" -f 4)
  #echo "genome $genome_id locus tag $locus_tag locus number $wspF_locus"
  # this is the wspF HMM
  # I assume there's only one line per locus tag, but just to be sure, I'm using head.
  
  # adjust the locus number according to $sign
  if [[ $sign = "pos" ]]; then
    wspE_locus=$((wspF_locus - 1)) # for posForward, subtract 1
  else
    wspE_locus=$((wspF_locus + 1)) # for negBackward, add 1
  fi
  
  # go into the wspE HMM file and look for a match with the appropriate locus and sign (strand)
  wspE_file="${datadir}/${genome_id}/hmmer_outputs/hmmer_output_PTHR43395.txt"
  match_count=$(cat "$wspE_file" | cut -f 1 -d " " |  awk -F'_' -v n="$wspE_locus" -v c="$contig" -v s="$sign" '$3 == c && $4 == n && $7 == s' | wc -l)
  if [[ $match_count > 0 ]]; then # write the synteny summary line to the output file. Since we filter by locus, contig, and sign, there should be exactly 1 hit, if any.
    echo "$match_count	$genome_id	$contig	$locus_tag	$wspF_locus	$wspE_locus	$sign" >> $log
    echo $line >> $outname
  fi
done

echo "Finished writing to $outname"