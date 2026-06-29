#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd

def write_summary(breseq_dirs, outname, n):
    outname_path = Path(outname)
    outname_path.parent.mkdir(exist_ok=True, parents=True)

    df_all = pd.DataFrame() # it's generally bad practice to grow a df, but this will give flexibility wrt columns of the output df
    with open(breseq_dirs) as f:
        for line in f:
            dir = line.strip()
            print(f"Processing {dir}")
            summary = f"{dir}/index.html"
            tab = pd.read_html(summary)
            df=tab[1]

            # column 0 of this table featured links to read alignments;
            # don't need these, so overwrite with breseq outdir name indicating where the variant came from
            if "Predicted mutations" in list(df.iloc[-2]): # penultimate row should just be "Predicted mutations" for each column
                name = "/".join(Path(dir).parts[-1 * (n+1):-1])
                df.columns = df.iloc[-1] # very last row contains the column names
                df = df.iloc[:-2] # as described, last two rows don't actually have data
                if "evidence" in df.columns:
                    df = df.rename(columns={'evidence': 'source'}) # evidence should be the 0th column. In any case, not needed here; overwrite with name.
                df["source"] = name # if no source column, will create one.
                
                df_all = pd.concat([df_all, df], axis=0, join='outer', ignore_index=True)

    df_all.to_csv(outname, sep="\t", index=False, header=True)
                
    print("Done!")

def main(breseq_dirs, outname, n):
    write_summary(breseq_dirs, outname, n)

if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        print(f"Usage: {sys.argv[0]} breseq_outdir_list.txt output.tsv [num_dirs_in_name_column]")
        sys.exit(1)

    breseq_dirs = sys.argv[1]
    outname = sys.argv[2]
    num = int(sys.argv[3]) if len(sys.argv) == 4 else 1

    main(breseq_dirs, outname, num)