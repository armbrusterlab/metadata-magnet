# Nextflow pipeline helper script
# If the blast db was
# If among the duplicated locus tags there is a title that contains that locus tag, keep that row
# Otherwise, keep all rows for that locus tag, and when left joined by the GUI, it will arbitrarily keep the first row (and the metadata table will have a few more rows than the FASTA has sequences)

library(dplyr)
library(stringr)

deduplicate <- function(metadata_file, outname) {
  df = read.csv(metadata_file, sep='\t', na.strings = character(), quote = "", colClasses = "character") # various things to prevent changes to the file beyond removal of rows
  # keep empty cells as "", not NA. quote = "" to preserve quotes in input file, colClasses = "character" so that it doesn't change how the scientific notation is written 
  
  duplicated_loci <- df |> 
    group_by(locus_tag) |> 
    summarize(n=n()) |> 
    arrange(desc(n)) |> 
    filter(n>1) |> 
    pull(locus_tag)
    
  # for each locus tag, determine whether there are any rows in which it appears in title_old
  # if so, we can drop all other rows, but if not, just keep all rows
  filtered_df <- df |> 
    filter(locus_tag %in% duplicated_loci) |> 
    mutate(match = str_detect(title_old, locus_tag))
    
  loci_to_deduplicate <- filtered_df |> 
    group_by(locus_tag) |> 
    summarize(n=n(), found = any(match)) |> 
    filter(found == TRUE) |> 
    pull(locus_tag)
    
  df_new <- df |> 
    filter(!((locus_tag %in% loci_to_deduplicate) & !str_detect(title_old, locus_tag))) # drop whichever rows we can, as explained earlier
    
  write.table(df_new, outname, sep='\t', row.names = FALSE, quote = FALSE, na = "")
  print("Completed!")
}


### Process CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 2) {
  stop('Please provide two arguments: <metadata_file> <outname>')
}

metadata_file <- args[1]
outname <- args[2]

cat("Metadata file provided:", metadata_file, "\n")
cat("Output file:", outname, "\n")

deduplicate(metadata_file, outname)