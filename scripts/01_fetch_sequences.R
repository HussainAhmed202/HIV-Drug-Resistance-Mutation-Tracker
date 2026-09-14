# 01_fetch_sequences.R
# Fetch real HIV-1 pol gene (protease + RT) sequences from NCBI GenBank

library(rentrez)

search_result <- entrez_search(
  db = "nuccore",
  term = "HIV-1[Organism] AND pol[Gene] AND protease AND reverse transcriptase AND Pakistan[Title]", 
  retmax = 40
)

# How many did we find, and what are the GI IDs?
print(search_result$count)
print(search_result$ids)


# Fetch the complete Genome(Fasta format) for the IDs we found
search_result_sequences <- entrez_fetch(
  db = "nuccore",
  id = search_result$ids,
  rettype = "fasta"
)

# Save raw FASTA to file
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
file <- file("data/raw/hiv_pol_sequences.fasta", open = "w")
writeLines(search_result_sequences, "data/raw/hiv_pol_sequences.fasta")
close(file)
cat("Saved", length(search_result$ids), "sequences to data/raw/hiv_pol_sequences.fasta\n")