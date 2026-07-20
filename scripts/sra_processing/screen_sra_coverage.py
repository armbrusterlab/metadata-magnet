from pathlib import Path
import pandas as pd
import json
import argparse

# example run:
# python "/home/kcw2/metadata-magnet/scripts/sra_processing/screen_sra_coverage.py" --taxdir "/home/kcw2/data/testing/testdir1/taxonomy_analysis/" --genome_length 6.26e6 --taxid 286

# # the function inputs
# taxdir = "/home/kcw2/data/testing/testdir1/taxonomy_analysis/"
# genome_length = 6.26e6 # scientific notation is saved as a float I think
# coverage_threshold = 30 # when parsing this, set the type as int or float
# taxid = 286 # int or string?

def join_data(taxdir, genome_length, taxid, coverage_threshold=30):
    esearch_file = f"{taxdir}/../metadata_esearch.csv"
    esearch = pd.read_csv(esearch_file, sep=",")

    pysradb_file = f"{taxdir}/../metadata_pysradb.tsv"
    pysradb = pd.read_csv(pysradb_file, sep="\t")

    coverage_list = []
    coverage_pass = []

    pathlist = Path(taxdir).rglob('*.json')
    for path in pathlist:
        run_id=path.stem # filename without extension- this is the run ID
        read_type = pysradb[pysradb.run_accession == run_id].reset_index(drop=True).library_layout[0]
        spot_factor = 1 + (read_type == "PAIRED") # in paired end runs, a single spot has two reads
        avg_length = esearch[esearch.Run == run_id].reset_index(drop=True).avgLength[0]
        # print(run_id, read_type, spot_factor, avg_length)

        # parse the json
        with open(path, 'r') as file:
            data = json.load(file)
            match = next((item for item in data[0]["tax_table"] if item["tax_id"] == taxid), None) # there should only be one match
            num_spots = 0 if match == None else match["total_count"]
        # print(num_spots)

        # calculate the coverage
        coverage = (spot_factor * num_spots * avg_length) / genome_length
        # print(coverage, coverage > coverage_threshold)

        # record pass/fail
        coverage_list.append(coverage)
        coverage_pass.append(coverage > coverage_threshold)

    esearch[f"coverage_taxid_{taxid}"] = coverage_list
    esearch["coverage_pass"] = coverage_pass
    df = pd.DataFrame({'run_id': esearch.Run, 'coverage': coverage_list, 'coverage_pass': coverage_pass})
    passed = df[df.coverage_pass == True].run_id
    df.head()

    df.to_csv(f"{taxdir}/../coverage_taxid_{taxid}.tsv", sep='\t', index=False) # may just name this coverage.tsv, though it's possible to use wildcards to find this tsv
    
    print("passed length:", len(passed))

    print(taxdir, genome_length, taxid, coverage_threshold)
    passed.to_csv(f"{taxdir}/../coverage_pass.txt", sep='\t', index=False, header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to filter SRA runs by coverage for a taxid of choice.")

    # positional arguments (required)
    parser.add_argument("-t", "--taxdir", type=str, help='Directory containing SRA taxonomy analyses.')
    parser.add_argument("-L", "--genome_length", type=float, help="Length of reference genome in bp.")
    parser.add_argument("-i", "--taxid", type=int, help="Taxonomic id to reference when calculating coverage.") # must be int or else it won't match in the JSON
    parser.add_argument("-c", "--coverage_threshold", type=float, default=30, help="Coverage required for the taxonomic ID.")

    args = parser.parse_args()
    join_data(args.taxdir, args.genome_length, args.taxid, args.coverage_threshold)