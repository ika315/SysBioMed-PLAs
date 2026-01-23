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


gmt_dir_to_csv <- function(gmt_dir) {
  files <- list.files(gmt_dir, full.names = TRUE)
  
  for (f in files) {
    
    if (!grepl("\\.gmt$", f, ignore.case = TRUE)) next
    
    genes <- qusage::read.gmt(f)[[1]]
    
    write.csv(
      data.frame(geneName = genes),
      file = sub("\\.gmt$", ".csv", f),
      row.names = FALSE, quote = FALSE
    )
  }
}

#gmt_dir_to_csv("/Users/irem/SysBioMed-PLAs/data/gmt_lists")

# immune gene set
# diff of platelets & leukocyte activation
# do not use

# -----------------------------
# immune gene set
# difference of platelet & leukocyte activation
# -----------------------------

# read gene lists
#updated_genes <- read_gene_list("data/updated_gene_list.csv")

#leukocyte_activation_genes <- read_gene_list(
#  "data/gmt_lists/GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs.csv"
#)

# COMPARING THE DIFFERENT PLATELET GENES

#files <- c(
#  "data/updated_gene_list.csv",
#  "data/gmt_lists/GNATENKO_PLATELET_SIGNATURE.v2025.1.Hs.csv",
#  "data/gmt_lists/GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs.csv",
#  "data/gmt_lists/GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs.csv",
#  "data/gmt_lists/HP_ABNORMAL_PLATELET_MEMBRANE_PROTEIN_EXPRESSION.v2025.1.Hs.csv",
#  "data/gmt_lists/MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs.csv",
#  "data/gmt_lists/MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_UP.v2025.1.Hs.csv",
#  "data/gmt_lists/REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION.v2025.1.Hs.csv",
#  "data/gmt_lists/WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS.v2025.1.Hs.csv"
#)

#gene_lists <- lapply(files, read_gene_list)
#names(gene_lists) <- basename(files)

#overlap <- matrix(0, length(gene_lists), length(gene_lists))
#rownames(overlap) <- colnames(overlap) <- names(gene_lists)

#for (i in seq_along(gene_lists)) {
#  for (j in seq_along(gene_lists)) {
#    overlap[i, j] <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
#  }
#}

#diag(overlap) <- NA
#overlap

#write_csv(
#  data.frame(List = rownames(overlap), overlap),
#  "data/overlap_gene_lists.csv"
#)
#unique_genes <- names(table(unlist(gene_lists)) == 1)


