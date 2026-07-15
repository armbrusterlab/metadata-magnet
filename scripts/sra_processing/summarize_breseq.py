#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd
import re

def write_summary(breseq_dirs, outname, n=1, filter_intergenic = False, filter_synonymous = False):
    outname_path = Path(outname)
    outname_path.parent.mkdir(exist_ok=True, parents=True)

    df_all = pd.DataFrame() # it's generally bad practice to grow a df, but this will give flexibility wrt columns of the output df
    with open(breseq_dirs) as f:
        for line in f:
            dir = line.strip()
            print(f"Processing {dir}")
            summary = f"{dir}/index.html"
            tab = pd.read_html(summary, extract_links="body") # need link locations from the evidence column, but this turns all elements into tuples
            df=tab[1]

            df.columns = [tup[0] for tup in df.iloc[-1]] # rename df columns so they're referenceable
            df["evidence_file"] = df["evidence"].apply(lambda x: x[1] if isinstance(x, tuple) else None) # get link locations

            # clean up columns now- no more tuples
            for col in df.columns:
                if col != "evidence_file":
                    df[col] = df[col].apply(lambda x: x[0] if isinstance(x, tuple) else x)

            name = "/".join(Path(dir).parts[-1 * (n+1):-1])
            df["source"] = name
                
            if list(df.iloc[-2])[0] != "Predicted mutations":
                print(f"Error: predicted mutations table was not retrieved for {dir}")
            else:
                df = df.iloc[:-2] # last two rows don't actually have data, so crop them.
                df_all = pd.concat([df_all, df], axis=0, join='outer', ignore_index=True)
    # now that all dfs have been joined, filter if requested.
    if filter_intergenic:
        print("Filtering out intergenic variants...")
        df_all = df_all[df_all['annotation'].str.contains('intergenic', na=False) == False]

    if filter_synonymous:
        print("Filtering out synonymous mutants...")
        df_all = df_all[df_all['annotation'].apply(find_synonymous) == False]

    df_all.to_csv(outname, sep="\t", index=False, header=True)
                
    print("Done!")

def find_synonymous(s):
    pattern = r"^([A-Z])[0-9]+([A-Z])?" # capturing groups: the two nucleotides flanking some number
    found = re.match(pattern, s)
    if found:
        return found.group(1) == found.group(2)
    else: # e.g. if it's an intergenic variant or blank
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to summarize breseq outputs.")

    # positional arguments (required)
    parser.add_argument("breseq_dirs", type=str, help='A file containing one breseq output dir per line. Each output dir must contain a breseq output file "index.html".')
    parser.add_argument("outname", type=str, help="The name of the file to write the results to.")

    # named (optional) arguments
    parser.add_argument("-n", "--num", type=int, default=1, help="Number of levels of dir names to preserve in the output (source column).")
    parser.add_argument("-i", "--filter_intergenic", action="store_true", help="Filter out intergenic variants.")
    parser.add_argument("-s", "--filter_synonymous", action="store_true", help="Filter out synonymous mutations.")

    args = parser.parse_args()
    write_summary(args.breseq_dirs, args.outname, args.num, args.filter_intergenic, args.filter_synonymous)