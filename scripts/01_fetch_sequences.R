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

  # Search Space 

  #   "LOCUS       MT222943                9002 bp    RNA     linear   VRL 13-MAY-2020\nDEFINITION  
  #   HIV-1 isolate DEURF15PK040 from Pakistan, complete genome.\nACCESSION   MT222943\nVERSION     
  #   MT222943.1\nKEYWORDS    .\nSOURCE      Human immunodeficiency virus 1 (HIV-1)\n  ORGANISM  
  #   Human immunodeficiency virus 1\n            Viruses; Riboviria; Pararnavirae; Artverviricota; 
  #   Revtraviricetes;\n            Ortervirales; Retroviridae; Orthoretrovirinae; Lentivirus;\n 
  #   Lentivirus humimdef1.\nREFERENCE   1  (bases 1 to 9002)\n  AUTHORS   Hora,B., Chen,Y., 
  #   Shah,S.A., Busch,M.P., Denny,T.N. and Gao,F.\n  TITLE     Characterization of near full-length 
  #   genome sequences for standard\n            panels of HIV-1 isolates established at the 
  #   External Quality\n            Assurance Program Oversight Laboratory (EQAPOL)\n  JOURNAL   
  #   Unpublished\nREFERENCE   2  (bases 1 to 9002)\n  AUTHORS   Hora,B., Chen,Y., Shah,S., 
  #   Busch,M., Denny,T. and Gao,F.\n  TITLE     Direct Submission\n  JOURNAL   Submitted ....
  #   .....gene="pol"\n     CDS             <.....    
  #   ....
  #   ....
  
  # Goal :: Look for the pol gene CDS start and end indexes defined in the Genbank file
  
  # Look at the highlighted portion below
                        # gene=\"pol\"\n     CDS             <1597..4608\n
  # start_idx of pol gene CDS = 1597
  # end_idx of pol gene CDS = 4608
  

      
  needle <- 'gene="pol"\n     CDS             <'
  needle_len <- nchar(needle)  # 33
  
    
  results <- data.frame(
    Pol_Start_Idx = integer(length(sequences_gb_list)),
    Pol_End_Idx = integer(length(sequences_gb_list))
  )
  
  for (i in seq_along(sequences_gb_list)) {
    
    search_space_end_idx <- nchar(sequences_gb_list[[i]])
    
    search_space_idx <- regexpr(needle, sequences_gb_list[[i]])  #3616
    
    if (idx[1] == -1) {
      warning(paste("pol CDS not found for record:", i))
      next
    }
    search_space_start_idx <- search_space_idx + needle_len

    remaining_seq <- substr(sequences_gb_list[[i]], search_space_start_idx, search_space_end_idx)  

    # find the first newline character in the new search space - this is where the end of my search space

    search_space_idx <- regexpr("\n",remaining_seq)

    if (search_space_idx[1] == -1) {
      warning(paste("End index not found for record:", i))
      next
    }

    search_space_end_idx <- search_space_start_idx + search_space_idx - 2 # 9
    
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

pol_gene <- substr(df$Fasta_Seq,df$Pol_Start_Idx,df$Pol_End_Idx)
pol_gene[["fasta_title"]] <- gsub("complete genome","pol gene complete sequence taken", df$Title)


# Save raw FASTA to file
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
file <- file("data/raw/hiv_pol_sequences_new.fasta", open = "w")

for (i in seq_len(length(df$Fasta_Seq))) {
  
  # Extract sequence
  sequence <- substr(
    df$Fasta_Seq[i],
    df$Pol_Start_Idx[i],
    df$Pol_End_Idx[i]
  )
  
  # Create FASTA header
  header <- gsub(
    "complete genome",
    "pol gene complete sequence taken",
    df$Title[i]
  )
  
  # Write header
  writeLines(
    paste0(">", header),
    file
  )
  
  # Write sequence
  writeLines(
    sequence,
    file
  )
}

close(file)

# writeLines(sequences_fasta, "data/raw/hiv_pol_sequences.fasta")
# cat("Saved", length(search_result$ids), "sequences to data/raw/hiv_pol_sequences.fasta\n")
cat("Saved", length(search_result$ids), "sequences to data/raw/hiv_pol_sequences_new.fasta\n")