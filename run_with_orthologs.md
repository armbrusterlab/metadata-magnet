In a [normal run](https://github.com/armbrusterlab/metadata-magnet/blob/main/pipeline_tutorial.md#example-run-setup-nextflow-r-shiny-gui) of this pipeline, you would obtain homologs of your target protein sequence via BLAST search. 
However, if you'd like to restrict your search to [_orthologs_](https://bio.libretexts.org/Bookshelves/Microbiology/Microbiology_(Boundless)/07%3A_Microbial_Genetics/7.13%3A_Bioinformatics/7.13C%3A_Homologs_Orthologs_and_Paralogs) of your target gene, we provide an alternative entry point to the pipeline. 
We haven't yet updated the Nextflow pipeline to be compatible with this alternative entry point, so for now you will need to run the steps manually.

Table of contents:
- [Obtain FASTA from OrthoDB](#obtain-fasta-from-orthodb)
- [Convert FASTA to expected .blast format](#convert-fasta-to-expected-blast-format)
- [Filtering](#filtering-optional)
- [Metadata retrieval](#retrieve-metadata)
- [Metadata processing](#process-metadata)
- [Convert back to FASTA](#convert-metadata-file-to-fasta)
- [Align FASTA](#alignment-optional)
- [Generate figures and statistical tests](#plug-into-gui)

# Retrieve orthologs
## Obtain FASTA from OrthoDB
First, search the name of your target protein (case-sensitive) on the [OrthoDB website](https://www.orthodb.org/?).
<img width="2496" height="1415" alt="image" src="https://github.com/user-attachments/assets/f5c46f2a-8812-4b23-9124-1c014865061b" />  

Expand the group you are interested in, and click "View protein FASTA". Right-click to Save As, and upload the file to the server on which you have cloned this pipeline. 
<img width="2471" height="479" alt="image" src="https://github.com/user-attachments/assets/4ade2341-617e-4aaf-809e-818cb29a716e" />

## Convert FASTA to expected .blast format
This step produces a file with the same output format as [run_local_blast.sh](https://github.com/armbrusterlab/metadata-magnet/blob/main/scripts/blast_processing/run_local_blast.sh). Run the following Python script in the bash command line. The first argument is the name of the FASTA downloaded from OrthoDB, and the second is the name of the output BLAST file.
```bash
conda activate metadata-magnet-env # set up from YAML in metadata-magnet/nextflow/envs

cd metadata-magnet/nextflow/example_data

my_blast="mucA_orthoDB_converted.blast"
python metadata-magnet/scripts/setup/orthoDB_fasta_to_blast.py "mucA_orthoDB.fasta" "$my_blast"
```
This step yields a .blast file which is compatible with scripts from this pipeline.

# Filtering (optional)
Most filtering scripts will not work, as they rely on the evalue provided by the BLAST search, which is not available if the data was obtained thorugh OrthoDB instead. However, you may still perform synteny-aware filtering if desired.

## Synteny search
Please refer to documentation on [synteny search inputs](https://github.com/armbrusterlab/metadata-magnet/blob/main/synteny_search_inputs.md), which includes an explanation of the intersection between synteny search hits and BLAST hits. It is recommended to run this in a tmux window, as it may take several hours depending on input dataset size.
```bash
conda activate pynteny-env # set up from YAML in metadata-magnet/nextflow/envs

# perform synteny search
bash metadata-magnet/scripts/synteny_wrapper.sh -g ${genomeDBsynteny} -i ${syntenyInput} -L ${hmmsList} -d ${hmmsDir} -m ${hmmsMetadata}

# intersect with .blast file made from OrthoDB fasta
bash "metadata-magnet/scripts/synteny_search/intersect_blast_and_synteny_pairwiseBlastApproach.sh" \
            -b ${my_blast} \
            -s ${syntenyInput} \
            -m \
            -p ${intersectPident} \
            -q ${intersectQcovs} \
            -c "locus"
```

# Retrieve metadata
You may skip this step if you performed synteny search, as the synteny summary file doubles as a metadata file.
Local genome database search is not possible because there are no genome accessions in the OrthoDB, so metadata will need to be retrieved via NCBI esearch.  
```bash
conda activate metadata-magnet-env # set up from YAML in metadata-magnet/nextflow/envs
python metadata-magnet/scripts/blast_processing/blast2gen.py "${my_blast}" "${my_blast%.*}_metadata.blast" "myemail@mail.com"
```

# Process metadata
This step annotates the metadata file with isolation source categories and subcategories. In this example, the -o (output filename) option is omitted, so the input file is overwritten with a version that features the new columns.
```bash
conda activate metadata-magnet-env # set up from YAML in metadata-magnet/nextflow/envs
categoryFile="metadata-magnet/nextflow/data/category_keywords.txt"
subcategoryFile="metadata-magnet/nextflow/data/subcategory_keywords.txt"

bash metadata-magnet/scripts/metadata_processing_wrapper.sh -f "${my_blast%.*}_metadata.blast" -c ${categoryFile} -s ${subcategoryFile}
```

# Convert metadata file to FASTA
Although the input file was a FASTA, its headers are not compatible with the GUI we have provided to produce figures and statistical tests. The easiest way to obtain a compatible FASTA is to convert the metadata file back to FASTA.  
You may skip this step if you performed synteny intersection, as that script produces a FASTA.
```bash
conda activate metadata-magnet-env # set up from YAML in metadata-magnet/nextflow/envs
bash metadata-magnet/scripts/blast_processing/convert_blast_to_fasta.sh -b "${my_blast}" -r -g -c "locus"
```
The output file has the same name as the input file, but instead has a .fasta extension.

# Alignment (optional)
You may use any alignment program, but for this example I will use MUSCLE/super5.
```bash
conda activate metadata-magnet-env # set up from YAML in metadata-magnet/nextflow/envs
my_fasta="mucA_orthoDB_converted.fasta"

seqCount=\$(grep '^>' ${my_fasta} | wc -l) # check for lines starting with >
echo "There are \$seqCount sequences in the input FASTA."

if [ "\$seqCount" -le 1 ]; then # can't align fasta with 0 or 1 sequences
    echo "Fewer than 2 sequences; not running MUSCLE."
    cp ${my_fasta} "${my_fasta%.*}_aligned.fasta"
elif [ "\$seqCount" -le 1000 ]; then # <= 1000 seqs: use full MUSCLE algorithm, since the size is reasonably small
    muscle -align ${input_file} -output "${my_fasta%.*}_aligned.fasta"
else # use super5 for larger datasets
    muscle -super5 ${input_file} -output "${my_fasta%.*}_aligned.fasta"
fi
```

# Plug into GUI
This step remains exactly as described [here](https://github.com/armbrusterlab/metadata-magnet/blob/main/pipeline_tutorial.md#analyze-nextflow-outputs-in-the-gui). The input files are the new FASTA (which may or may not be aligned) and the metadata file.
