import re
from Bio import SeqIO
import os
import pandas as pd

def find_longest_match(p, s):
	'''
	Input: a compiled regex pattern p, and a string s to search.
	Output: the length of the longest match for p.
	'''
	matches = re.findall(p, s)
	if len(matches) > 0:
		longest_match_length = len(max(matches, key=len))
		return longest_match_length
	else: # no matches
		return 0

def make_repeat_pattern(motif):
	'''
	Input: a motif string
	Output: a regular expression that searches for repeat sections of the motif
	'''
	pattern = f"(?:{motif.upper()})+"
	return re.compile(pattern)

def summarize_repeats(metadata_file, motif, seq_col="sequence", outname="", remove_gaps = True):
	'''
	Inputs: a metadata in which to find repeats, a motif for repeats, 
	(optional) name of the column containing sequences,
	(optional, by default overwrite the input file) an output file name,
	(optional) whether to remove gap characters ("-") in the sequence
	Output: To outname, writes the length of the longest repeat per sequence to the "longest_repeat" column.
	If a sequence has no instances of the repeat at all, a length of 0 is recorded.
	Repeats are not case-sensitive (i.e. motif and sequences are converted to uppercase before searching)
	'''
	if outname == "":
		outname = metadata_file # overwrite the input file
	make_outdir(outname) # make outdir if necessary

	p = make_repeat_pattern(motif) # set up the regex object
	df = pd.read_csv(metadata_file, sep="\t")

	if remove_gaps:
		longest_repeat_notNormalized = [find_longest_match(p, s.upper().replace("-", "")) for s in df[seq_col]]
	else:
		longest_repeat_notNormalized = [find_longest_match(p, s.upper()) for s in df[seq_col]] # gets lengths as number of characters
	longest_repeat = [int(n / len(motif)) for n in longest_repeat_notNormalized] # normalize to get the number of times the motif appears in the repeat section

	df['longest_repeat'] = longest_repeat

	df.to_csv(outname, sep='\t', index=False)
	print(f"Done writing to {outname}")


def make_outdir(outname):
	'''
	General helper function.
	Input: name of output file.
	Output: if the directory in outname doesn't exist, create it.
	'''
	# get outdir path from outname
	directory_path = os.path.dirname(outname)

	# Check if the directory exists, and create it if it doesn't
	if not os.path.exists(directory_path):
		os.makedirs(directory_path)
		print(f"Directory '{directory_path}' created.")