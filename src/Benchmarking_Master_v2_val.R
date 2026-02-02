# -------------------------------------------------------------------
# Benchmarking_Master_v2.R - FINAL VERSION
# -------------------------------------------------------------------
library(Seurat)
library(AUCell)
library(UCell)
library(pROC)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(mclust)

args <- commandArgs(trailingOnly = TRUE)
METHOD_NAME    <- if(length(args) >= 1) args[1] else "AUCell"
SIG_NAME       <- if(length(args) >= 2) args[2] else "MANNE_DN"
SIG_FILE_BASE  <- if(length(args) >= 3) args[3] else "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs"
USE_EXTENSION  <- if(length(args) >= 4) as.logical(args[4]) else TRUE
THRESH_MODE    <- if(length(args) >= 5) args[5] else "gmm_dist_dual" 

# --- KONFIGURATION ---
GT_COLUMN    <- "pla.status"
POSITIVE_VAL <- "PLA"
PATH_DATA    <- "~/SysBioMed-PLAs/data/seu_sx_integration_new.rds"

# --- ORDNERSTRUKTUR AUTOMATISCH ERSTELLEN ---
OUT_DIR <- paste0("plots/Platelet_Main_Val/", SIG_NAME, "/", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, "/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("results_val/metrics", recursive = TRUE, showWarnings = FALSE)
dir.create("results_val/celltype_data", recursive = TRUE, showWarnings = FALSE)
dir.create("results_val/extended_lists_val", recursive = TRUE, showWarnings = FALSE)

# --- DATEN LADEN ---
print(paste("Lade Daten für:", SIG_NAME, "mit", METHOD_NAME))
pbmc <- readRDS(PATH_DATA)

new_metadata <- read.csv("~/SysBioMed-PLAs/data/external_dataset_pla-status_metatable.csv", row.names = 1)
rownames(new_metadata) <- new_metadata$barcodes_clean
common_cells <- intersect(Cells(pbmc), rownames(new_metadata))
pbmc <- subset(pbmc, cells = common_cells)
pbmc <- AddMetaData(pbmc, metadata = new_metadata[common_cells, ])

# --- GENLISTE LADEN ---
base_dir <- getwd()
source(file.path(base_dir, "src", "read_and_extend_gene_list.R"))
PATH_SIG <- file.path(base_dir, "data", paste0(SIG_FILE_BASE, ".csv"))
genes <- read_gene_list(PATH_SIG)

# --- Immune Config --- 
IMMUNE_SIG <- "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs"
PATH_IMMUNE_SIG <- file.path(base_dir, "data", paste0(IMMUNE_SIG, ".csv"))
immune_genes <- read_gene_list(PATH_IMMUNE_SIG)
immune_genes <- intersect(immune_genes, rownames(pbmc))

donors <- unique(pbmc$donor)

tp <- numeric(length(donors))
fn <- numeric(length(donors))
tn <- numeric(length(donors))
fp <- numeric(length(donors))
prec <- numeric(length(donors))
rec <- numeric(length(donors))
f1 <- numeric(length(donors))

for (i in 1:length(donors)) {
    print(paste("---- Cross-Validation Fold:", i, "of", length(donors), "----"))
    test_donor_id <- donors[i]
    other_donors <- donors[-i]

    train_pbmc <- subset(pbmc, subset = donor %in% other_donors)
    test_pbmc  <- subset(pbmc, subset = donor == test_donor_id)

    # --- SCORING LOGIK ---
    print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

    if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
        train_expression_matrix <- GetAssayData(train_pbmc, layer = "data")
        train_rankings <- AUCell_buildRankings(train_expression_matrix, plotStats=FALSE)
        train_auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), train_rankings)
        train_pbmc$Raw_Score_Original <- as.numeric(getAUC(train_auc_orig)[1, ])

        train_auc_imm <- AUCell_calcAUC(list(Immune_Score = immune_genes), train_rankings)
        train_pbmc$Immune_Score <- as.numeric(getAUC(train_auc_imm)[1, ])


        expression_matrix <- GetAssayData(test_pbmc, layer = "data")
        rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
        auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
        test_pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])

        auc_imm <- AUCell_calcAUC(list(Immune_Score = immune_genes), rankings)
        test_pbmc$Immune_Score <- as.numeric(getAUC(auc_imm)[1, ])

    } else if (METHOD_NAME == "UCell") {
        train_pbmc <- AddModuleScore_UCell(train_pbmc, features = list(Platelet_Orig = genes), name = NULL)
        train_pbmc$Raw_Score_Original <- train_pbmc$Platelet_Orig
        train_pbmc <- AddModuleScore_UCell(train_pbmc, features = list(Immune_Score = immune_genes), name = NULL)
        train_pbmc$Immune_Score <- train_pbmc$Immune_Score 

        test_pbmc <- AddModuleScore_UCell(test_pbmc, features = list(Platelet_Orig = genes), name = NULL)
        test_pbmc$Raw_Score_Original <- test_pbmc$Platelet_Orig
        test_pbmc <- AddModuleScore_UCell(test_pbmc, features = list(Immune_Score = immune_genes), name = NULL)
        test_pbmc$Immune_Score <- test_pbmc$Immune_Score 

    } else if (METHOD_NAME == "AddModuleScore") {
        train_pbmc <- AddModuleScore(train_pbmc, features = list(genes), name = "AMS_Orig")
        train_pbmc$Raw_Score_Original <- train_pbmc$AMS_Orig1

        train_pbmc <- AddModuleScore(train_pbmc, features = list(immune_genes), name = "AMS_Immune")
        train_pbmc$Immune_Score <- train_pbmc$AMS_Immune1

        test_pbmc <- AddModuleScore(test_pbmc, features = list(genes), name = "AMS_Orig")
        test_pbmc$Raw_Score_Original <- test_pbmc$AMS_Orig1

        test_pbmc <- AddModuleScore(test_pbmc, features = list(immune_genes), name = "AMS_Immune")
        test_pbmc$Immune_Score <- test_pbmc$AMS_Immune1
    }

    # 2. Schritt: Gensequenz-Erweiterung (optional)
    if (USE_EXTENSION) {
        TRAIN_EXT_FILE <- paste0("results_val/extended_lists_val/ext_", SIG_NAME, "_", METHOD_NAME, "_", test_donor_id, "_test.csv")
        TEST_EXT_FILE <- paste0("results_val/extended_lists_val/ext_", SIG_NAME, "_", METHOD_NAME, "_", test_donor_id, "_val.csv")
        if (file.exists(TRAIN_EXT_FILE)) {
            message("Lade existierende Liste...")
            train_extended_genes <- read.csv(TRAIN_EXT_FILE)$geneName
            test_extended_genes <- read.csv(TEST_EXT_FILE)$geneName
        } else {
            message("Berechne neue Extension...")
            res_ext <- extend_gene_set(train_pbmc, base_genes = genes, score_name = "Raw_Score_Original")
            train_extended_genes <- res_ext$extended_genes
            write.csv(data.frame(geneName = train_extended_genes), TRAIN_EXT_FILE, row.names = FALSE)

            res_ext <- extend_gene_set(test_pbmc, base_genes = genes, score_name = "Raw_Score_Original")
            test_extended_genes <- res_ext$extended_genes
            write.csv(data.frame(geneName = test_extended_genes), TEST_EXT_FILE, row.names = FALSE)
        }
        train_final_genes <- train_extended_genes
        test_final_genes <- test_extended_genes
    } else {
        train_final_genes <- genes
        test_final_genes <- genes
    }

    # 3. Schritt: Finales Scoring
    if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
        train_pbmc$Raw_Score <- as.numeric(getAUC(AUCell_calcAUC(list(Platelet_Score = train_final_genes), train_rankings))[1, ])
        test_pbmc$Raw_Score <- as.numeric(getAUC(AUCell_calcAUC(list(Platelet_Score = test_final_genes), rankings))[1, ])
    } else if (METHOD_NAME == "UCell") {
        train_pbmc <- AddModuleScore_UCell(train_pbmc, features = list(Platelet_Score = train_final_genes), name = NULL)
        train_pbmc$Raw_Score <- train_pbmc$Platelet_Score

        test_pbmc <- AddModuleScore_UCell(test_pbmc, features = list(Platelet_Score = test_final_genes), name = NULL)
        test_pbmc$Raw_Score <- test_pbmc$Platelet_Score
    } else {
        train_pbmc <- AddModuleScore(train_pbmc, features = list(train_final_genes), name = "AMS")
        train_pbmc$Raw_Score <- train_pbmc$AMS1

        test_pbmc <- AddModuleScore(test_pbmc, features = list(test_final_genes), name = "AMS")
        test_pbmc$Raw_Score <- test_pbmc$AMS1
    }

    train_pbmc$Z_Score <- as.vector(scale(train_pbmc$Raw_Score))
    train_pbmc$Immune_Z <- as.vector(scale(train_pbmc$Immune_Score))

    test_pbmc$Z_Score <- as.vector(scale(test_pbmc$Raw_Score))
    test_pbmc$Immune_Z <- as.vector(scale(test_pbmc$Immune_Score))
    # --- THRESHOLDING & EVALUIERUNG ---
    train_pbmc$GT_Response <- ifelse(train_pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
    roc_obj <- roc(response = train_pbmc$GT_Response, predictor = train_pbmc$Z_Score, direction = "<", quiet = TRUE)

    test_pbmc$GT_Response <- ifelse(test_pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)

    THRESHOLD_I <- -Inf

    if (THRESH_MODE == "youden") {
        THRESHOLD_Z <- as.numeric(coords(roc_obj, x = "best", best.method = "youden")$threshold)
    } else if (THRESH_MODE == "percentile") {
        prob_cutoff <- 0.90 
        THRESHOLD_Z <- as.numeric(quantile(train_pbmc$Z_Score, probs = prob_cutoff))
    } else if (THRESH_MODE == "manual") {
        # MAD-Ansatz: Robust gegen Ausreißer (die PLAs)
        med <- median(train_pbmc$Z_Score)
        mad_val <- mad(train_pbmc$Z_Score)
        THRESHOLD_Z <- med + (1.5 * mad_val) 
    } else if (THRESH_MODE == "null_dist_platelet") {
        z_grid <- unique(quantile(train_pbmc$Z_Score, probs = seq(0.05, 0.95, 0.05)))
        best_f1 <- -1
        for(tz in z_grid) {
            pred <- train_pbmc$Z_Score > tz
            tp <- sum(pred & train_pbmc$GT_Response == 1); fp <- sum(pred & train_pbmc$GT_Response == 0)
            fn <- sum(!pred & train_pbmc$GT_Response == 1); prec <- tp/(tp+fp); rec <- tp/(tp+fn)
            f1 <- 2*(prec*rec)/(prec+rec)
            if(!is.na(f1) && f1 > best_f1) { best_f1 <- f1; THRESHOLD_Z <- tz }
        }
    } else if (THRESH_MODE == "null_dist_immune_dual") {
        # 2D-Optimierung (Platelet + Immune)
        z_grid <- unique(quantile(train_pbmc$Z_Score, probs = seq(0.1, 0.9, 0.1)))
        i_grid <- unique(quantile(train_pbmc$Immune_Z, probs = seq(0.1, 0.9, 0.1)))
        best_f1 <- -1
        for(tz in z_grid) {
            for(ti in i_grid) {
                pred <- (train_pbmc$Z_Score > tz) & (train_pbmc$Immune_Z > ti)
                tp <- sum(pred & train_pbmc$GT_Response == 1); fp <- sum(pred & train_pbmc$GT_Response == 0)
                fn <- sum(!pred & train_pbmc$GT_Response == 1); prec <- tp/(tp+fp); rec <- tp/(tp+fn)
                f1 <- 2*(prec*rec)/(prec+rec)
                if(!is.na(f1) && f1 > best_f1) { 
                    best_f1 <- f1; THRESHOLD_Z <- tz; THRESHOLD_I <- ti 
                }
            }
        }
    }else if (THRESH_MODE == "gmm_dist_platelet" || THRESH_MODE == "gmm_dist_dual") {
        auc_obs <- train_pbmc$Raw_Score
        auc_imm <- train_pbmc$Immune_Score

        fit_plat <- Mclust(auc_obs, G = 2)
        plat_high <- which.max(fit_plat$parameters$mean)
        train_pbmc$Platelet_High <- fit_plat$classification == plat_high
        train_pbmc$Immune_High <- TRUE

        if(THRESH_MODE == "gmm_dist_dual") {
            idx <- which(train_pbmc$Platelet_High)
            fit_imm <- Mclust(auc_imm[idx], G = 2)
            imm_high <- which.max(fit_imm$parameters$mean)

            train_pbmc$Immune_High <- FALSE
            train_pbmc$Immune_High[idx] <- fit_imm$classification == imm_high
            THRESHOLD_I <- min(train_pbmc$Immune_Z[train_pbmc$Immune_High], na.rm = TRUE)
        }

        THRESHOLD_Z <- min(train_pbmc$Z_Score[train_pbmc$Platelet_High], na.rm = TRUE)
    } else if (THRESH_MODE == "kmeans"){
        set.seed(42)
        km_data <- FetchData(train_pbmc, vars = c("Z_Score", "Immune_Z")) %>% drop_na()
        km_fit <- kmeans(km_data, centers = 4, nstart = 50, iter.max = 100)

        train_pbmc$KMeans_Cluster <- NA
        train_pbmc$KMeans_Cluster[rownames(km_data)] <- km_fit$cluster

        cluster_stats <- km_data %>pbmc%
            mutate(Cluster = km_fit$cluster) %>%
            group_by(Cluster) %>%
            summarise(
                Mean_Z = mean(Z_Score),
                Mean_Immune_Z = mean(Immune_Z),
                Score_Sum = Mean_Z + Mean_Immune_Z,
                Min_Z = min(Z_Score),        # Untergrenze für Z-Score
                Min_Immune = min(Immune_Z)   # Untergrenze für Immune Score
            )

        positive_cluster <- cluster_stats$Cluster[which.max(cluster_stats$Score_Sum)]

        train_pbmc$Platelet_High <- train_pbmc$KMeans_Cluster == positive_cluster
        train_pbmc$Immune_High <- train_pbmc$KMeans_Cluster == positive_cluster
        THRESHOLD_Z <- cluster_stats$Min_Z[cluster_stats$Cluster == positive_cluster]
        THRESHOLD_I <- cluster_stats$Min_Immune[cluster_stats$Cluster == positive_cluster]
    }


    positive_condition <- if (THRESH_MODE == "kmeans") {
        train_pbmc$KMeans_Cluster == positive_cluster

    } else if (THRESH_MODE == "gmm_dist_platelet") {
        train_pbmc$Platelet_High

    } else if (THRESH_MODE == "gmm_dist_dual") {
        test_pbmc$Z_Score > THRESHOLD_Z & test_pbmc$Immune_Z > THRESHOLD_I


    } else {
        (train_pbmc$Z_Score > THRESHOLD_Z) &
        (train_pbmc$Immune_Z > THRESHOLD_I)
    }

    test_pbmc$Prediction <- factor(ifelse(positive_condition, "Positive", "Negative"),levels = c("Negative","Positive"))

    test_pbmc$Error_Type <- case_when(
        test_pbmc$Prediction == "Positive" & test_pbmc[[GT_COLUMN]] == POSITIVE_VAL ~ "TP",
        test_pbmc$Prediction == "Positive" & test_pbmc[[GT_COLUMN]] != POSITIVE_VAL ~ "FP",
        test_pbmc$Prediction == "Negative" & test_pbmc[[GT_COLUMN]] == POSITIVE_VAL ~ "FN",
        test_pbmc$Prediction == "Negative" & test_pbmc[[GT_COLUMN]] != POSITIVE_VAL ~ "TN"
    )
    
    # Zeit und Metriken
    tp[i] <- sum(test_pbmc$Error_Type == "TP", na.rm = TRUE); fp[i] <- sum(test_pbmc$Error_Type == "FP", na.rm = TRUE)
    fn[i] <- sum(test_pbmc$Error_Type == "FN", na.rm = TRUE); tn[i] <- sum(test_pbmc$Error_Type == "TN", na.rm = TRUE)
    prec[i] <- if((tp[i] + fp[i]) > 0) tp[i] / (tp[i] + fp[i]) else 0
    rec[i]  <- if((tp[i] + fn[i]) > 0) tp[i] / (tp[i] + fn[i]) else 0
    f1[i]   <- if((prec[i] + rec[i]) > 0) 2 * prec[i] * rec[i] / (prec[i] + rec[i]) else 0

    gc()
}

# Ergebnisse speichern
metrics_df <- data.frame(
    Donor = donors,
    TP = tp,    
    FP = fp,
    FN = fn,
    TN = tn,
    Precision = prec,
    Recall = rec,
    F1_Score = f1
)

metrics_file <- paste0("results_val/metrics/metrics_", SIG_NAME, "_", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, ".csv")
write.csv(metrics_df, file = metrics_file, row.names = FALSE)

avg_metrics <- data.frame(
    Metric = c("TP", "FP", "FN", "TN", "Precision", "Recall", "F1_Score"),
    Average = c(mean(tp), mean(fp), mean(fn), mean(tn), mean(prec), mean(rec), mean(f1))
)
avg_metrics_file <- paste0("results_val/metrics/avg_metrics_", SIG_NAME, "_", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, ".csv")
write.csv(avg_metrics, file = avg_metrics_file, row.names = FALSE)


# Plot Avg. Metrics
avg_metrics_plot <- ggplot(avg_metrics, aes(x = Metric, y = Average)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    ggtitle(paste("Average Metrics -", SIG_NAME, "-", METHOD_NAME, "- Ext:", USE_EXTENSION, "- Thresh:", THRESH_MODE)) +
    ylab("Average Value") +
    xlab("Metric")
ggsave(filename = paste0(OUT_DIR, "avg_metrics_", SIG_NAME, "_", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, ".png"),
       plot = avg_metrics_plot, width = 8, height = 6)

# Plot Each F1 Score
f1_plot <- ggplot(metrics_df, aes(x = Donor, y = F1_Score)) +
    geom_bar(stat = "identity", fill = "coral") +
    theme_minimal() +
    ggtitle(paste("F1 Score per Donor -", SIG_NAME, "-", METHOD_NAME, "- Ext:", USE_EXTENSION, "- Thresh:", THRESH_MODE)) +
    ylab("F1 Score") +
    xlab("Donor")   

ggsave(filename = paste0(OUT_DIR, "f1_scores_", SIG_NAME, "_", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, ".png"),
       plot = f1_plot, width = 8, height = 6)