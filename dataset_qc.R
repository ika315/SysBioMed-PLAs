library(Seurat)
library(DoubletFinder)
# read data
data_path <- "data/GSM5008737_RNA_3P/"
counts <- Read10X(data.dir = data_path)

seu <- CreateSeuratObject(
  counts = counts,
  project = "PlateletEnrichment",
  min.cells = 3,
  min.features = 50 
)

rm(counts)

# statistics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")lot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"),
  ncol = 5,
  pt.size = 0
)

FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.hb")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.ribo")

seu <- subset(
  seu,
  subset =  nFeature_RNA > 100 &
    percent.mt < 20   
)

seu <- NormalizeData(seu)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu <- ScaleData(seu)

seu <- RunPCA(seu)

dims.use <- 1:20

seu <- FindNeighbors(seu, dims = dims.use)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = dims.use)

DimPlot(seu, group.by = "seurat_clusters", label = TRUE)

# doublet detection
sweep.res.list <- paramSweep(seu, PCs = dims.use, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_optimal <- bcmvn$pK[which.max(bcmvn$BCmetric)]

doublet.rate <- 0.075
nExp_poi <- round(doublet.rate * ncol(seu))

df_seu <- seu

# Remove unnecessary assays, scale.data and dimreductions
df_seu <- DietSeurat(
  df_seu,
  counts    = TRUE,
  data      = TRUE,      
  scale.data = FALSE,
  assays    = "RNA",    
  dimreducs = "pca"     
)


DimPlot(seu, group.by = grep("DF.classifications", colnames(seu@meta.data), value = TRUE), 
        reduction = "umap")

