# read and extend gene list

library(readr)
library(dplyr)

read_gene_list <- function(gene_csv) {
  
  use_cols <- c("geneName")
  
  gene_list <- read_csv(gene_csv, show_col_types = FALSE) %>%
    select(all_of(use_cols))
  
  # find duplicates
  duplicates <- gene_list %>%
    filter(duplicated(geneName) | duplicated(geneName, fromLast = TRUE))
  
  genes <- unique(gene_list$geneName)
  
  cat(length(genes), "genes were extracted from file.\n")
  cat(nrow(duplicates), "duplicates were found. These include:\n")
  print(unique(duplicates$geneName))
  
  return(genes)   # a character vector
}

extend_gene_set <- function(gene_set) {
  # to be implemented
}


base_dir <- dirname(getwd())
print(base_dir)
gene_csv <- file.path(base_dir, "data", "updated_gene_list.csv")
genes <- read_gene_list(gene_csv)

print(genes)