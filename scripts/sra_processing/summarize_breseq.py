#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd

def write_summary(breseq_dirs, outname, n):
    outname_path = Path(outname)
    outname_path.parent.mkdir(exist_ok=True, parents=True)
    with open(outname, "w") as out: # first overwrite the file if it exists
        out.write("source	seq_id	position	mutation	annotation	gene	description\n")

    with open(outname, "a") as out: # using write mode for this part seemed to cause issues
        with open(breseq_dirs) as f:
            for line in f:
                dir = line.strip()
                print(f"Processing {dir}")
                summary = f"{dir}/index.html"
                tab = pd.read_html(summary)
                df=tab[1]

                # column 0 of this table featured links to read alignments;
                # don't need these, so overwrite with breseq outdir name indicating where the variant came from
                name = "/".join(Path(dir).parts[-1 * (n+1):-1])
                df = df.iloc[0:-2]
                df[0] = name
                
                df.to_csv(outname, mode="a", sep="\t", index=False, header=False)
                
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