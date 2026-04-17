#!/usr/bin/env python3

import json
import re
import sys

def parse_orthoDB(fasta, outname):
    '''
    Input: FASTA obtained from OrthoDB website
    (Headers include JSON after sequence title, and there are empty lines between entries)
    Output: .blast file mimicking the format expected by pipeline: columns for genomeid, proteinid, sequence, evalue, title, organism
    '''
    with open(fasta, 'r') as f:
        text=f.read().split(">")

    p = r'[\w:]*\s' # intended to match the sequence name in the header, including the space between it and the JSON

    with open(outname, 'w') as out:
        for i in range(len(text)):
            if text[i] != "":
                # parse the entry
                entry = (text[i].replace(re.search(p, text[i]).group(), "")).split('\n') # remove fasta title
                json_string, seq = entry[0], entry[1]
                data = json.loads(json_string)

                # write the info to the output file
                out.write(f"0\t{data['pub_gene_id']}\t{seq}\tNA\t{data['og_name']}\t{data['organism_name']}\n")

    print(f"Finished writing to {outname}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} from_orthoDB.fasta output.blast")
        sys.exit(1)

    fasta = sys.argv[1]
    outname = sys.argv[2]

    parse_orthoDB(fasta, outname)