# read and extend gene list

library(readr)
library(dplyr)

# METHOD THAT READS IN CREATED CSV FILE OF KNOWN PLATELET SIGNATURE GENES
read_gene_list <- function(gene_csv) {
  
  use_cols <- c("geneName")
  
  gene_list <- read_csv(gene_csv, show_col_types = FALSE) %>%
    select(all_of(use_cols))
  
  # might have some duplicates in the list
  duplicates <- gene_list %>%
    filter(duplicated(geneName) | duplicated(geneName, fromLast = TRUE))
  
  genes <- unique(gene_list$geneName)
  
  cat(length(genes), "genes were extracted from file.\n")
  print(unique(duplicates$geneName))
  
  return(genes)
}

# METHOD THAT EXTENDS EXISTING GENE SET WITH ENRICHED GENES FROM ENRICHMENT METHODS
extend_gene_set <- function(genes=genes, gene_set) {
  extended_gene_set <- union(genes, gene_set)
  return(extended_gene_set)
}

# METHOD THAT CREATES SIGNATURE OBJECT FOR THE ENRICHMENT METHODS OUT OF THE KNOWN GENES
make_signature <- function(genes, sig_name = "Signature", method = "other") {
  # defalt methods (addmodulescore etc)
  
  if (method != "vision") {
    sig_list <- list()
    sig_list[[sig_name]] <- genes
    return(sig_list)
  }
  
  # for vision
  if (method == "vision") {

    sigData <- rep(1, length(genes))
    names(sigData) <- genes
    
    sig_obj <- VISION::createGeneSignature(
      name = sig_name,
      sigData = sigData
    )
    
    return(sig_obj)
  }
}

# RUNNING THE METHODS
# read in genes
base_dir <- getwd()
gene_csv <- file.path(base_dir, "data", "updated_gene_list.csv")
genes <- read_gene_list(gene_csv)

print(genes)

# create sig lists for enrichment methods
sig_other <- make_signature(genes, "Platelet_Signature", "other")
sig_vision <- make_signature(genes = genes, sig_name = "Platelet_Signature", method = "vision")

