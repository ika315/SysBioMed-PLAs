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
extend_gene_set <- function(
    pbmc,
    base_genes,
    score_name,
    top_n = 50,
    high_quantile = 0.9,
    min.pct = 0.05,
    test.use = "wilcox") {
  score_vec <- pbmc[[score_name, drop = TRUE]]
  
  high_cells <- names(score_vec)[
    score_vec > quantile(score_vec, high_quantile)
  ]
  
  group_col <- paste0(score_name, "_group")
  pbmc[[group_col]] <- ifelse(
    colnames(pbmc) %in% high_cells,
    "high",
    "low"
  )
  
  markers <- FindMarkers(
    pbmc,
    ident.1 = "high",
    ident.2 = "low",
    group.by = group_col,
    logfc.threshold = 0, # test all genes
    min.pct = min.pct, # genes must be expressed in at least 5% of cells in either group
    test.use = test.use
  )
  
  # outputs a dataframe with one row per gene
  # avg_log2FC -> effect size
  # p_val & p_val_adj
  # pct.1 -> % expressed in high cells
  # pct.2 -> % expressed in low cells
  
  top_genes <- rownames(
    markers[order(markers$avg_log2FC, decreasing = TRUE), ]
  )[1:top_n]
  
  
  extended <- union(base_genes, top_genes)
  
  return(list(
    extended_genes = extended,
    top_genes = top_genes,
    markers = markers
  ))
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

