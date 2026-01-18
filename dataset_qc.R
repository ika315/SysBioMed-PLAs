# TODO Violin Plot - COunt Distribution per Sample pre/post filtering
# Vln Mitrochondirial content per sample pre/post filtering
# nach scDblFinder UMAP plot colored singlet/doublet
# Distribution Highly Variable Genes
# res 0.5, 1.0, 1.5, 2.0 Clustree plot (UMAP)
# Annotation by Cell Type plotten
# Standard QC Plots aus Tutorial

# Keine Thresholds, approach aus scanpy tutorial

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SoupX)
library(Matrix)
library(matrixStats)
library(clustree)
library(scry)
library(ggplot2)

save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )

    ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}


# read data
data_path <- "data/GSM5008737_RNA_3P/"
plot_dir <- "plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
counts <- Read10X(data.dir = data_path)
counts <- as(counts, "dgCMatrix")

seu <- CreateSeuratObject(
  counts = counts,  
  project = "PlateletEnrichment",
  min.cells = 3,
  min.features = 20 
)


# add metadata
metadata_path <- "data/GSE164378_sc.meta.data_3P.csv"
metadata <- read.csv(metadata_path, row.names = 1)
seu_ids <- colnames(seu)

meta_use <- metadata[seu_ids, , drop = FALSE]
seu <- AddMetaData(seu, metadata = meta_use)
seu$donor_time <- interaction(seu$donor, seu$time, drop = TRUE, sep = "_")

#set.seed(42)
#cells_test <- sample(colnames(seu), 10000)
#seu <- seu[, cells_test]

# prepare soupX
tod <- counts
#tod <- tod[Matrix::rowSums(tod) > 0, , drop = FALSE]
#tod <- tod[, Matrix::colSums(tod) > 2, drop = FALSE]  

flt <- GetAssayData(seu, layer = "counts")            
flt <- flt[rownames(tod), , drop = FALSE]    

#rm(counts)
gc()

# Violin: counts + mt pre-QC


# statistics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

hb_genes <- grep("^HB", rownames(seu), value = TRUE)
hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]
seu[["percent.hb"]] <- PercentageFeatureSet(seu, features = hb_genes)

plot_qc_prepost <- function(seu, prefix, group.by = "donor_time") {
  # Violin: absolute QC metrics per sample
  p <- VlnPlot(
    seu,
    features = c("nCount_RNA", "nFeature_RNA"),
    group.by = group.by,
    pt.size = 0,
    ncol = 2
  )

  save_plot(p, paste0("QC_", prefix, "_violin_counts_features.png"), width = 10, height = 5)

  # Violin: relative QC metrics per sample
  p <- VlnPlot(
    seu,
    features = c("percent.mt", "percent.ribo", "percent.hb"),
    group.by = group.by,
    pt.size = 0,
    ncol = 3
  )

  save_plot(p, paste0("QC_", prefix, "_violin_percent_mt_ribo_hb.png"), width = 10, height = 5)

  # Scatter: nCount vs mt
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    ggtitle(paste0(prefix, ": nCount_RNA vs percent.mt"))
  save_plot(p, paste0("QC_", prefix, "_scatter_nCount_vs_percentMT.png"))

  # Scatter: nCount vs ribo
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.ribo") +
    ggtitle(paste0(prefix, ": nCount_RNA vs percent.ribo"))
  save_plot(p, paste0("QC_", prefix, "_scatter_nCount_vs_percentRibo.png"))

  # Scatter: nCount vs hb
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.hb") +
    ggtitle(paste0(prefix, ": nCount_RNA vs percent.hb"))
  save_plot(p, paste0("QC_", prefix, "_scatter_nCount_vs_percentHB.png"))
}

plot_qc_prepost(seu, "pre_QC")

saveRDS(seu, file = "data/seu_before_filtering.rds")
seu_preQC <- readRDS("data/seu_before_filtering.rds")

# qc start

is_outlier_mad <- function(x, nmads = 5, type = c("both", "lower", "upper")){
  type <- match.arg(type)

  med <- median(x, na.rm = TRUE)
  madv <- mad(x, constant = 1, na.rm = TRUE)

  if (madv == 0 || is.na(madv)) {
    return (rep(FALSE, length(x)))
  }

  lower <- med - nmads * madv
  upper <- med + nmads * madv

  if(type == "both") {
    return (x < lower | x > upper)
  } else if(type == "lower") {
    return (x < lower)
  } else {
    return (x > upper)
  }
}

out_low_counts <- is_outlier_mad(seu$nCount_RNA, nmads = 5, type = "both")
out_low_genes <- is_outlier_mad(seu$nFeature_RNA, nmads = 5, type = "both")

out_high_mt <- is_outlier_mad(seu$percent.mt, nmads = 5, type = "upper")
out_high_ribo <- is_outlier_mad(seu$percent.ribo, nmads = 5, type = "upper")
out_high_hb <- is_outlier_mad(seu$percent.hb, nmads = 5, type = "upper")

seu$qc_outlier <- out_low_counts | out_low_genes | out_high_mt | out_high_hb | out_high_ribo

seu <- subset(seu, subset = !qc_outlier)

plot_qc_prepost(seu, "post_QC")
# cluster for soup
dims.use <- 1:20

seu <- NormalizeData(seu)

# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

sce <- as.SingleCellExperiment(seu)
hvgs <- devianceFeatureSelection(sce)
dev <- rowData(hvgs)$binomial_deviance
names(dev) <- rownames(hvgs)  # Gene als Names

## 4) Nach Deviance sortieren (absteigend) und Top 2000 Gene nehmen
hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

## 5) Diese Gen-Namen als VariableFeatures im Seurat-Objekt setzen
VariableFeatures(seu) <- hvgs_genes

seu <- ScaleData(seu, features = hvgs_genes) # vars to regress?

seu <- RunPCA(seu, features = hvgs_genes)

seu <- FindNeighbors(seu, dims = dims.use)

seu <- FindClusters(seu, resolution = 0.5)

soupx_groups <- seu$seurat_clusters
names(soupx_groups) <- colnames(seu)
soupx_groups <- soupx_groups[colnames(flt)]

saveRDS(seu, file = "data/seu_before_soup.rds")
# soup

common_cells <- intersect(colnames(flt), colnames(seu))
flt <- flt[, common_cells, drop = FALSE]
tod <- tod[rownames(flt), , drop = FALSE]  # tod mit flt syncen (gleiche gene wie flt)

sc <- SoupChannel(tod, flt, calcSoupProfile = FALSE)

gene_counts <- Matrix::rowSums(tod)
gene_counts <- gene_counts[gene_counts > 0]

soupProf <- data.frame(
  est    = as.numeric(gene_counts / sum(gene_counts)),
  counts = as.numeric(gene_counts),
  row.names = names(gene_counts)
)


sc <- setSoupProfile(sc, soupProf)

sc <- setClusters(sc, soupx_groups)
sc <- autoEstCont(sc, doPlot = TRUE)
gc()
corrected <- adjustCounts(sc, roundToInt = TRUE)

rm(tod, flt, seu, soupx_groups, keep) 
gc()

seu_sx <- CreateSeuratObject(counts = corrected)

rm(seu)
gc()


cnt <- GetAssayData(seu_sx, assay = "RNA", layer = "counts")
keep <- Matrix::rowSums(cnt > 0) >= 20
seu_sx <- seu_sx[keep, ]

# qc after soup

seu_sx[["percent.mt"]] <- PercentageFeatureSet(seu_sx, pattern = "^MT-")
seu_sx[["percent.ribo"]] <- PercentageFeatureSet(seu_sx, pattern = "^RP[SL]")
hb_genes <- grep("^HB", rownames(seu_sx), value = TRUE)
hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]
seu_sx[["percent.hb"]] <- PercentageFeatureSet(seu_sx, features = hb_genes)

out_low_counts <- is_outlier_mad(seu_sx$nCount_RNA, nmads = 5, type = "both")
out_low_genes <- is_outlier_mad(seu_sx$nFeature_RNA, nmads = 5, type = "both")

out_high_mt <- is_outlier_mad(seu_sx$percent.mt, nmads = 5, type = "upper")
out_high_ribo <- is_outlier_mad(seu_sx$percent.ribo, nmads = 5, type = "upper")
out_high_hb <- is_outlier_mad(seu_sx$percent.hb, nmads = 5, type = "upper")

seu_sx$qc_outlier <- out_low_counts | out_low_genes | out_high_mt | out_high_hb | out_high_ribo
seu_sx <- subset(seu_sx, subset = !qc_outlier)

saveRDS(seu_sx, file = "data/seu_after_soup.rds")
# clustering

seu_sx <- NormalizeData(seu_sx, normalization.method = "LogNormalize",)

# seu_sx <- FindVariableFeatures(seu_sx, selection.method = "vst", nfeatures = 2000)


sce <- as.SingleCellExperiment(seu_sx)
hvgs_sce <- devianceFeatureSelection(sce)
dev <- rowData(hvgs_sce)$binomial_deviance
names(dev) <- rownames(hvgs_sce)
hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]
VariableFeatures(seu_sx) <- hvgs_genes
gc()
seu_sx <- ScaleData(seu_sx, features = hvgs_genes) # vars.to.regress = c("nCount_RNA", "percent.mt")) # todo: andere features hinzufügen?

seu_sx <- RunPCA(seu_sx, features = hvgs_genes)

seu_sx <- FindNeighbors(seu_sx, dims = dims.use)

seu_sx <- FindClusters(seu_sx, resolution = 0.5)

# print(DimPlot(seu_sx, group.by = "seurat_clusters", label = TRUE))

# add metadata
metadata_path <- "data/GSE164378_sc.meta.data_3P.csv"
metadata <- read.csv(metadata_path, row.names = 1)
seu_ids <- colnames(seu_sx)

meta_use <- metadata[seu_ids, , drop = FALSE]
seu_sx <- AddMetaData(seu_sx, metadata = meta_use)
seu_sx$donor_time <- interaction(seu_sx$donor, seu_sx$time, drop = TRUE, sep = "_")

saveRDS(seu_sx, file = "data/seu_before_dd.rds")
# doublet detection
sce <- as.SingleCellExperiment(seu_sx)
sce <- scDblFinder(sce, clusters = "seurat_clusters", samples = "donor_time")

seu_sx$scDblFinder_class <- sce$scDblFinder.class
seu_sx$scDblFinder_score <- sce$scDblFinder.score

seu_sx <- RunUMAP(seu_sx, dims = dims.use)

p <- DimPlot(
  seu_sx,
  group.by = "scDblFinder_class",
  reduction = "umap"
) + ggtitle("scDblFinder: singlet/doublet")

save_plot(p, "UMAP_scDblFinder_class.png")

p <- FeaturePlot(
  seu_sx,
  features = "scDblFinder_score",
  reduction = "umap"
) + ggtitle("scDblFinder score")

save_plot(p, "UMAP_scDblFinder_score.png")

table(seu_sx$scDblFinder_class)
saveRDS(seu_sx, file = "data/seu_after_dd.rds")

#print(
  #DimPlot(seu_sx, group.by = "scDblFinder_class", reduction = "umap")
#)

#VlnPlot(seu_sx, "scDblFinder_score", group.by = "scDblFinder_class", pt.size = 0.1)

#FeaturePlot(seu_sx, "scDblFinder_score", reduction = "umap")

seu_sx <- subset(seu_sx, subset = scDblFinder_class == "singlet")

# final clustering

seu_sx <- NormalizeData(seu_sx, normalization.method = "LogNormalize",)

# seu_sx <- FindVariableFeatures(seu_sx, selection.method = "vst", nfeatures = 2000)


sce <- as.SingleCellExperiment(seu_sx)
gc()
hvgs_sce <- devianceFeatureSelection(sce)
gc()
dev <- rowData(hvgs_sce)$binomial_deviance
names(dev) <- rownames(hvgs_sce)

hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

dev_df <- data.frame(deviance = dev)
p <- ggplot(dev_df, aes(x = deviance)) +
  geom_histogram(bins = 50) +
  ggtitle("Distribution of binomial deviance (all genes)")

save_plot(p, "HVG_deviance_distribution.png")

VariableFeatures(seu_sx) <- hvgs_genes

seu_sx <- ScaleData(seu_sx, features = hvgs_genes) #, vars.to.regress = c("percent.mt"))

seu_sx <- RunPCA(seu_sx, features = hvgs_genes)

dev_sorted <- sort(dev, decreasing = TRUE)
df_rank <- data.frame(
  rank = seq_along(dev_sorted),
  deviance = as.numeric(dev_sorted),
 # qc start

out_low_counts <- is_outlier_mad(seu_sx$nCount_RNA, nmads = 5, type = "both")
out_low_genes <- is_outlier_mad(seu_sx$nFeature_RNA, nmads = 5, type = "both")

out_high_mt <- is_outlier_mad(seu_sx$percent.mt, nmads = 5, type = "upper")
out_high_ribo <- is_outlier_mad(seu_sx$percent.ribo, nmads = 5, type = "upper")
out_high_hb <- is_outlier_mad(seu_sx$percent.hb, nmads = 5, type = "upper")

seu_sx$qc_outlier <- out_low_counts | out_low_genes | out_high_mt | out_high_hb | out_high_ribo

seu_sx <- subset(seu_sx, subset = !qc_outlier)

# cluster for soup
dims.use <- 1:20

seu <- NormalizeData(seu)

# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

sce <- as.SingleCellExperiment(seu)
hvgs <- devianceFeatureSelection(sce)
dev <- rowData(hvgs)$binomial_deviance
names(dev) <- rownames(hvgs)  # Gene als Names

## 4) Nach Deviance sortieren (absteigend) und Top 2000 Gene nehmen
hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

## 5) Diese Gen-Namen als VariableFeatures im Seurat-Objekt setzen
VariableFeatures(seu) <- hvgs_genes

seu <- ScaleData(seu, features = hvgs_genes) # vars to regress?

seu <- RunPCA(seu, features = hvgs_genes)

seu <- FindNeighbors(seu, dims = dims.use)

seu <- FindClusters(seu, resolution = 0.5)

soupx_groups <- seu$seurat_clusters
names(soupx_groups) <- colnames(seu)
soupx_groups <- soupx_groups[colnames(flt)]

saveRDS(seu, file = "data/seu_before_soup.rds") is_hvg = seq_along(dev_sorted) <= 2000
)

p <- ggplot(df_rank, aes(x = rank, y = deviance)) +
  geom_point(aes(alpha = is_hvg), size = 0.6) +
  ggtitle("Deviance feature selection: rank vs deviance") +
  xlab("Gene rank (1 = highest deviance)") +
  ylab("Binomial deviance")

save_plot(p, "HVG_deviance_rankplot.png")

p <- ggplot(df_rank, aes(x = rank, y = deviance)) +
  geom_point(size = 0.4, alpha = 0.6) +
  geom_vline(xintercept = 2000, linetype = "dashed") +
  ggtitle("HVG selection by deviance: rank plot") +
  xlab("Gene rank (1 = highest deviance)") +
  ylab("Binomial deviance")

save_plot(p, "HVG_deviance_rankplot.png")


df_rank <- data.frame(
  rank = seq_along(dev_sorted),
  deviance = as.numeric(dev_sorted),
  HVG = seq_along(dev_sorted) <= 2000
)

top_n <- 10
top_genes <- names(dev_sorted)[1:top_n]

df_rank$gene <- names(dev_sorted)
df_rank$is_top10 <- df_rank$gene %in% top_genes

p <- ggplot(df_rank, aes(x = rank, y = deviance, color = HVG)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = 2000, linetype = "dashed") +
  geom_text_repel(
    data = subset(df_rank, is_top10),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  ggtitle("HVG selection by deviance (top 2000 highlighted)") +
  xlab("Gene rank") +
  ylab("Binomial deviance")

save_plot(p, "HVG_deviance_rankplot_highlight_top10.png")

top_n <- 20
top_genes <- head(names(sort(dev, decreasing = TRUE)), top_n)

df_top <- data.frame(
  gene = names(dev),
  deviance = as.numeric(dev),
  is_top = names(dev) %in% top_genes
)

p <- ggplot(df_top, aes(x = deviance)) +
  geom_histogram(bins = 80) +
  ggtitle("Deviance distribution (top genes flagged)")

save_plot(p, "HVG_deviance_hist_topflag.png")

saveRDS(seu_sx, file = "data/seu_after_last_pca.rds")

seu_sx <- FindNeighbors(seu_sx, dims = dims.use)

seu_sx <- FindClusters(seu_sx, resolution = c(0.5, 1.0, 1.5, 2.0)) # use leiden? currently louvain

p <- clustree(seu_sx, prefix = "RNA_snn_res.")
save_plot(p, "clustree_resolutions.png", width = 8, height = 8)

seu_sx <- RunUMAP(seu_sx, dims = dims.use)

p <- DimPlot(seu_sx, group.by = "RNA_snn_res.0.5",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 0.5")

save_plot(p, "UMAP_res_0.5.png")

p <- DimPlot(seu_sx, group.by = "RNA_snn_res.1",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 1.0")

save_plot(p, "UMAP_res_1.png")


p <- DimPlot(seu_sx, group.by = "RNA_snn_res.1.5",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 1.5")

save_plot(p, "UMAP_res_1.5.png")

p <- DimPlot(seu_sx, group.by = "RNA_snn_res.2",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 2.0")

save_plot(p, "UMAP_res_2.png")

p <- DimPlot(
  seu_sx,
  group.by = "celltype.l2",
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cell type annotation (level 2)")

save_plot(p, "UMAP_celltype_L2.png")

p <- DimPlot(
  seu_sx,
  group.by = "celltype.l1",
  reduction = "umap",
  label = TRUE
) + ggtitle("Cell classes (level 1)")

save_plot(p, "UMAP_celltype_L1.png")

p <- DimPlot(
  seu_sx,
  group.by = "celltype.l3",
  reduction = "umap",
  label = FALSE
) + ggtitle("Cell subtypes (level 3)")

save_plot(p, "UMAP_celltype_L3.png")

seu_sx[["percent.mt"]] <- PercentageFeatureSet(seu_sx, pattern = "^MT-")
p <- VlnPlot(
  seu_sx,
  features = c("nCount_RNA", "percent.mt"),
  group.by = "donor_time",
  pt.size = 0,
  ncol = 2
) + ggtitle("Post-QC: nCount_RNA & percent_final.mt per sample")

save_plot(p, "QC_post_violin_counts_mt.png", width = 10, height = 5)

p <- VlnPlot(
  seu_sx,
  features = "percent.mt",
  group.by = "donor_time",
  pt.size = 0,
  ncol = 1
) + ggtitle("Post-QC: percent_final.mt per sample")

save_plot(p, "QC_post_violin_mt.png", width = 10, height = 5)

seu_preQC[["percent.mt"]] <- PercentageFeatureSet(seu_preQC, pattern = "^MT-")
p <- VlnPlot(
  seu_preQC,
  features = "percent.mt",
  group.by = "donor_time",
  pt.size = 0,
  ncol = 1
) + ggtitle("Pre-QC: percent_final.mt per sample")

save_plot(p, "QC_pre_violin_mt_10.png", width = 10, height = 5)

seu_preQC <- subset(
  seu_preQC,
  subset = percent.mt < 10 &
           percent.hb < 5
)

p <- VlnPlot(
  seu_preQC,
  features = "nCount_RNA",
  group.by = "donor_time",
  pt.size = 0,
  ncol = 1
) + ggtitle("Pre-QC: nCount_RNA")

save_plot(p, "QC_pre_violin_count.png", width = 10, height = 5)

saveRDS(seu_sx, file = "data/seu_sx_final_new.rds")
seu_sx <- readRDS("data/seu_sx_final_new.rds")
