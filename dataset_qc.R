library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SoupX)
library(Matrix)
library(clustree)
library(scry)


# read data
data_path <- "data/GSM5008737_RNA_3P/"

counts <- Read10X(data.dir = data_path)
counts <- as(counts, "dgCMatrix")

seu <- CreateSeuratObject(
  counts = counts,  
  project = "PlateletEnrichment",
  min.cells = 3,
  min.features = 20 
)

#set.seed(42)
#cells_test <- sample(colnames(seu), 10000)
#seu <- seu[, cells_test]

# prepare soupX
tod <- counts
tod <- tod[Matrix::rowSums(tod) > 0, , drop = FALSE]
tod <- tod[, Matrix::colSums(tod) > 2, drop = FALSE]  

flt <- GetAssayData(seu, layer = "counts")            
flt <- flt[rownames(tod), , drop = FALSE]    

#rm(counts)
gc()

# statistics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

hb_genes <- grep("^HB", rownames(seu), value = TRUE)
hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]
seu[["percent.hb"]] <- PercentageFeatureSet(seu, features = hb_genes)

#VlnPlot(
##  seu,
#  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#  ncol = 3,
# pt.size = 0
#)

#VlnPlot(seu_sx, features = c("nCount_RNA","nFeature_RNA","percent.mt"), group.by = "lane", pt.size = 0)
#VlnPlot(seu_sx, features = c("nCount_RNA","nFeature_RNA","percent.mt"), group.by = "donor", pt.size = 0)


##FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
#FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.hb")
#FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.ribo")
# qc start
seu <- subset(
  seu,
  subset = percent.mt < 20 &
           percent.hb < 5
)

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

# TODO annotate platelet clusters

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

#rm(tod, flt, seu, soupx_groups, keep) 
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

seu_sx <- subset(
  seu_sx,
  subset = percent.mt < 20 &
           percent.hb < 5
)

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

# seu_sx <- RunUMAP(seu_sx, dims = dims.use)

# print(DimPlot(seu_sx, group.by = "seurat_clusters", label = TRUE))

# add metadata
metadata_path <- "data/GSE164378_sc.meta.data_3P.csv"
metadata <- read.csv(metadata_path, row.names = 1)
seu_ids <- colnames(seu_sx)

meta_use <- metadata[seu_ids, , drop = FALSE]
seu_sx$donor_time <- interaction(seu_sx$donor, seu_sx$time, drop = TRUE, sep = "_")
seu_sx <- AddMetaData(seu_sx, metadata = meta_use)

# doublet detection
sce <- as.SingleCellExperiment(seu_sx)
sce <- scDblFinder(sce, clusters = "seurat_clusters", samples = "donor_time")

seu_sx$scDblFinder_class <- sce$scDblFinder.class
seu_sx$scDblFinder_score <- sce$scDblFinder.score

table(seu_sx$scDblFinder_class)

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

hvgs_sce <- devianceFeatureSelection(sce)
gc()
dev <- rowData(hvgs_sce)$binomial_deviance
names(dev) <- rownames(hvgs_sce)

hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

VariableFeatures(seu_sx) <- hvgs_genes

seu_sx <- ScaleData(seu_sx, features = hvgs_genes) #, vars.to.regress = c("percent.mt")) # todo: andere features hinzufügen?

seu_sx <- RunPCA(seu_sx, features = hvgs_genes)

seu_sx <- FindNeighbors(seu_sx, dims = dims.use)

seu_sx <- FindClusters(seu_sx, resolution = 0.5)

seu_sx <- RunUMAP(seu_sx, dims = dims.use)

print(DimPlot(seu_sx, group.by = "seurat_clusters", label = TRUE))

saveRDS(seu_sx, file = "data/seu_sx_final.rds")
# seu_sx <- readRDS("data/seu_sx_final.rds")
