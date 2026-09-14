# 01_fetch_sequences.R
# Fetch real HIV-1 pol gene (protease + RT) sequences from NCBI GenBank

library(rentrez)
library(Biostrings)
# Search NCBI nucleotide database for HIV-1 pol gene sequences
# Complete CDS not found for pol so searching for entire genome. Then pick pol gene CDS.

#  term = "HIV-1[Organism] AND pol[Gene] AND complete[SLEN]",
search_result <- entrez_search(
  db = "nuccore",
  term = "Human immunodeficiency virus 1[Organism] AND complete genome[Title] AND Pakistan[Title]", 
  retmax = 40
)

# How many did we find, and what are the GI IDs?
print(search_result$count)
print(search_result$ids)


search_result_summary <- entrez_summary(db="nuccore",id=search_result$ids)
accession_ids <- character(length(search_result$ids))
titles <- character(length(search_result$ids))

for (i in seq_along(search_result$ids)) {
  
  accession_ids[i] <- search_result_summary[[i]]$accessionversion
  titles[i] <- search_result_summary[[i]]$title
  
}

df <- data.frame(
  Accession_ID = accession_ids,
  Title = titles
)

sequences_fa_list <- vector(mode = "character", length = length(search_result$ids))
sequences_fa_headers <- vector(mode = "character", length = length(search_result$ids))
sequences_fa_code <- vector(mode = "character", length = length(search_result$ids))


# Fetch the complete Genome(Fasta format) for the IDs we found
for (i in seq_along(search_result$ids)) {
  sequences_fa_list[[i]] <- entrez_fetch(
    db = "nuccore",
    id = search_result$ids[[i]],
    rettype = "fasta"
  )
}

for (i in seq_along(sequences_fa_list)) {
  sequences_fa_headers[i] <- substr(sequences_fa_list[i],1,regexpr("\n",sequences_fa_list[i])-1)
}
# print(sequences_fa_headers)
df[["Fasta_Header"]] <- sequences_fa_headers


for (i in seq_along(sequences_fa_list)) {
  sequences_fa_code[i] <- gsub("\n", "", substr(sequences_fa_list[1],regexpr("\n",sequences_fa_list[1])+1,nchar(sequences_fa_list[1])))
}
df[["Fasta_Seq"]] <- sequences_fa_code
# print(sequences_fa_code)




# df[["Genome_FA"]] <- sequences_fa_list

# To extract the pol cDS, need to look at the start and end of the CDS found in the GenBank format
sequences_gb_list <- vector(length = length(search_result$ids))
for (i in seq_along(search_result$ids)) {
  sequences_gb_list[[i]] <- entrez_fetch(
    db = "nuccore",
    id = search_result$ids[[i]],
    rettype = "gb"
  )
}


extract_pol_indices <- function(sequences_gb_list) {

  needle <- 'gene="pol"\n     CDS             <'
  needle_len <- nchar(needle)  # 33
  
    
  results <- data.frame(
    Pol_Start_Idx = integer(length(sequences_gb_list)),
    Pol_End_Idx = integer(length(sequences_gb_list))
  )
  
  for (i in seq_along(sequences_gb_list)) {
    
    search_space_end_idx <- nchar(sequences_gb_list[[i]])
    
    idx <- regexpr(needle, sequences_gb_list[[i]])  #3616
    
    if (idx[1] == -1) {
      warning(paste("pol CDS not found for record:", i))
      next
    }
    search_space_start_idx <- idx + needle_len

    remaining_seq <- substr(sequences_gb_list[[i]], search_space_start_idx, search_space_end_idx)  

    # find the first newline character in the new search space - this is where the end of my search space

    idx <- regexpr("\n",remaining_seq)

    if (idx[1] == -1) {
      warning(paste("End index not found for record:", i))
      next
    }

    search_space_end_idx <- search_space_start_idx + idx - 2 # 9
    
    results[["Pol_Start_Idx"]][i] <- strsplit(
      substr(sequences_gb_list[[i]], search_space_start_idx, search_space_end_idx),  #1597..4608
      '\\.\\.')[[1]][1]  # 1597
    

    results[["Pol_End_Idx"]][i] <- strsplit(
      substr(sequences_gb_list[i], search_space_start_idx, search_space_end_idx),  #1597..4608
      '\\.\\.')[[1]][2]   # 4608

  }
  
  return(results)
}

df <- c(df, extract_pol_indices(sequences_gb_list))


# 1    MT222943.1 HIV-1 isolate DEURF15PK040 from Pakistan, complete genome          3649        3658
# 2    MT222942.1 HIV-1 isolate DEMA115PK021 from Pakistan, complete genome          3649        3658
# 3    KX232629.1        HIV-1 isolate PK040 from Pakistan, complete genome          3533        3542


# Save raw FASTA to file
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
writeLines(sequences_fasta, "data/raw/hiv_pol_sequences.fasta")

cat("Saved", length(search_result$ids), "sequences to data/raw/hiv_pol_sequences.fasta\n")