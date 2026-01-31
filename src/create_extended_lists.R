library(Seurat)
library(AUCell)
library(UCell)
library(pROC)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(mclust)


# --- KONFIGURATION ---
GT_COLUMN    <- "pla.status"
POSITIVE_VAL <- "PLA"
PATH_DATA    <- "~/SysBioMed-PLAs/data/seu_sx_integration_new.rds"

 # --- DATEN LADEN ---
pbmc <- readRDS(PATH_DATA)
new_metadata <- read.csv("~/SysBioMed-PLAs/data/external_dataset_pla-status_metatable.csv", row.names = 1)
rownames(new_metadata) <- new_metadata$barcodes_clean
common_cells <- intersect(Cells(pbmc), rownames(new_metadata))
pbmc <- subset(pbmc, cells = common_cells)
pbmc <- AddMetaData(pbmc, metadata = new_metadata[common_cells, ])

# --- Immune Config --- 
base_dir <- getwd()
IMMUNE_SIG <- "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs"
PATH_IMMUNE_SIG <- file.path(base_dir, "data", paste0(IMMUNE_SIG, ".csv"))
immune_genes <- read_gene_list(PATH_IMMUNE_SIG)
immune_genes <- immune_genes <- intersect(immune_genes, rownames(GetAssayData(pbmc, layer = "data")))

signatures <- list(
    "GNATENKO_PLATELET_SIGNATURE.v2025.1.Hs",
    "GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs",
    "HP_ABNORMAL_PLATELET_MEMBRANE_PROTEIN_EXPRESSION.v2025.1.Hs",
    "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs",
    "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_UP.v2025.1.Hs",
    "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION.v2025.1.Hs",
    "updated_gene_list",
    "WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS.v2025.1.Hs"
)

signature_names <- list(
    "GNATENKO",
    "GOBP_REG",
    "HP_ABNORMAL",
    "MANNE_DN",
    "MANNE_UP",
    "REACTOME",
    "Standard",
    "WP_PLATELET"
)


for(i in seq_along(signatures)){

    SIG_NAME <- signature_names[[i]]
    SIG_FILE_BASE <- signatures[[i]]

    # --- GENLISTE LADEN ---
    print(paste("Processing Signature:", SIG_NAME))
    base_dir <- getwd()
    source(file.path(base_dir, "src", "read_and_extend_gene_list.R"))
    PATH_SIG <- file.path(base_dir, "data", paste0(SIG_FILE_BASE, ".csv"))
    genes <- read_gene_list(PATH_SIG)
    for(METHOD_NAME in c("AUCell", "UCell", "AddModuleScore")){
        print(paste("Processing Signature:", SIG_NAME, "with Method:", METHOD_NAME))
        # --- ORDNERSTRUKTUR AUTOMATISCH ERSTELLEN ---
        dir.create("results/extended_lists_new", recursive = TRUE, showWarnings = FALSE)

    

        # --- SCORING LOGIK ---
        print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

        if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
            expression_matrix <- GetAssayData(pbmc, layer = "data")
            rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
            auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
            pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])

        } else if (METHOD_NAME == "UCell") {
            pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Orig = genes), name = NULL)
            pbmc$Raw_Score_Original <- pbmc$Platelet_Orig

        } else if (METHOD_NAME == "AddModuleScore") {
            pbmc <- AddModuleScore(pbmc, features = list(genes), name = "AMS_Orig")
            pbmc$Raw_Score_Original <- pbmc$AMS_Orig1
        }

        EXT_FILE <- paste0("results/extended_lists_new/ext_", SIG_NAME, "_", METHOD_NAME, ".csv")
        if (file.exists(EXT_FILE)) {
            message("Lade existierende Liste...")
            extended_genes <- read.csv(EXT_FILE)$geneName
        } else {
            message("Berechne neue Extension...")
            res_ext <- extend_gene_set(pbmc, base_genes = genes, score_name = "Raw_Score_Original")
            extended_genes <- res_ext$extended_genes
            write.csv(data.frame(geneName = extended_genes), EXT_FILE, row.names = FALSE)
        }
    }
}
