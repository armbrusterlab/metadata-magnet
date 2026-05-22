from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from pathlib import Path

def extract_gene(gb_file, outname, overwrite = True, asm = None, target = None, target_type = None, buffer_upstream=0, buffer_downstream=0, translate=False):
    '''
    gb_file = path to GenBank.
    To outname, write a FASTA either containing just the gene, or if either target or target_type are None, all genes in FASTA.
    If overwrite, clears the contents of outname; otherwise continue to build upon FASTA located at outname.

    asm: assembly accession (i.e. name of dir to look in when searching genome_db for GBFF to retrieve metadata from later); a failsafe for inference from gb_file path

    target = query name; target_type indicates where in the record to search for the target
    (e.g. target = "mucA", target_type = "gene"; 
    target = "NP_249454.1", target_type = "protein_id"; 
    target = "PA0763", target_type = "locus_tag")
    Please note that if the target matches multiple sequences, then multiple sequences will be written to outname.

    Buffers are defined in terms of bp.

    If using translate==True, then please use buffer_upstream==0, as translation starts at the start location of the sequence.
    '''
    assembly_accession = asm if asm is not None else Path(gb_file).parent.name # even if it's not actually an assembly accession, want the dir name so we can use the fasta's seq id to find the original gb_file later
    # in the context of Nextflow this might be unreliable, but it can serve as a backup
    # it turns out that it's unnecessarily complicated to get the assembly accession from a GenBank file

    if overwrite:
        with open(outname, "w"):
            pass

    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if target is None or target_type is None or feature.qualifiers.get(target_type) == [target]:
                    # Get location, handling strand direction
                    loc = feature.location
                    # Extract start/end with flanking buffer
                    start = max(0, loc.start - buffer_upstream)
                    end = min(len(record.seq), loc.end + buffer_downstream)
                    # Slice and extract
                    flanked_seq = record.seq[start:end]
                    if translate:
                        flanked_seq = flanked_seq.translate(to_stop=True)

                    locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                    protein_name = feature.qualifiers.get("product", ["N/A"])[0]
                    cds_seq = feature.extract(record.seq)

                    # Create a new SeqRecord
                    new_record = SeqRecord(
                        flanked_seq,
                        id=f"{assembly_accession}-{protein_id}-{locus_tag} | {protein_name} |",
                        description=f"buffer: upstream +{loc.start-start} bp, downstream +{end-loc.end} bp"
                    )

                    # Save to FASTA
                    with open(outname, "a") as output_handle:
                        SeqIO.write(new_record, output_handle, "fasta")
