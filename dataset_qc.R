library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SoupX)
library(Matrix)
library(clustree)
# read data
data_path <- "/nfs/home/students/f.mathis/SysBioMed-PLAs/data/GSM5008737_RNA_3P/"
counts <- Read10X(data.dir = data_path)
counts <- as(counts, "dgCMatrix")

seu <- CreateSeuratObject(
  counts = counts,  
  project = "PlateletEnrichment",
  min.cells = 3,
  min.features = 50 
)

tod <- counts
flt <- seu@assays$RNA@counts

rm(counts)
gc()

# statistics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^HB(?!P)", regex = TRUE)

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.hb")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.ribo")

# qc start
seu <- subset(
  seu,
  subset = percent.mt < 20 &
           percent.ribo < 40 &
           percent.hb < 5
)

dims.use <- 1:20

seu <- NormalizeData(seu)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu <- ScaleData(seu)

seu <- RunPCA(seu)

seu <- FindNeighbors(seu, dims = dims.use)

seu <- FindClusters(seu, resolution = 0.5)

soupx_groups <- seu$seurat_clusters

sc <- SoupChannel(tod, flt, calcSoupProfile = TRUE)

sc <- setClusters(sc, soupx_groups)
sc <- autoEstCont(sc, doPlot = TRUE)

corrected <- adjustCounts(sc, roundToInt = TRUE)

seu_sx <- CreateSeuratObject(counts = corrected)

keep <- Matrix::rowSums(seu_sx@assays$RNA@counts > 0) >= 20
seu_sx <- seu_sx[keep, ]

rm(tod, flt, seu, soupx_groups, keep) 
gc()

seu_sx[["percent.mt"]] <- PercentageFeatureSet(seu_sx, pattern = "^MT-")
seu_sx[["percent.ribo"]] <- PercentageFeatureSet(seu_sx, pattern = "^RP[SL]")
seu_sx[["percent.hb"]] <- PercentageFeatureSet(seu_sx, pattern = "^HB(?!P)", regex = TRUE)

seu_sx <- subset(
  seu_sx,
  subset = percent.mt < 20 &
           percent.ribo < 40 &
           percent.hb < 5
)

seu_sx <- NormalizeData(seu_sx)

seu_sx <- FindVariableFeatures(seu_sx, selection.method = "vst", nfeatures = 2000)

seu_sx <- ScaleData(seu_sx)

seu_sx <- RunPCA(seu_sx)

seu_sx <- FindNeighbors(seu_sx, dims = dims.use)

seu_sx <- FindClusters(seu_sx, resolution = 1.0)

seu_sx <- RunUMAP(seu_sx, dims = dims.use)

DimPlot(seu_sx, group.by = "seurat_clusters", label = TRUE)

# doublet detection
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce)

seu$scDblFinder_class <- sce$scDblFinder.class
seu$scDblFinder_score <- sce$scDblFinder.score

table(seu$scDblFinder_class)

print(
  DimPlot(seu, group.by = "scDblFinder_class", reduction = "umap")
)

VlnPlot(seu, "scDblFinder_score", group.by = "scDblFinder_class", pt.size = 0.1)

FeaturePlot(seu, "scDblFinder_score", reduction = "umap")

pbmc_singlets <- subset(pbmc_seu, subset = scDblFinder_class == "singlet")