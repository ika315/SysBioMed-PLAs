library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SoupX)
library(Matrix)
library(clustree)
library(scry) # TODO implementieren statt normaler feature selection (findvariablefeatures)

# read data
data_path <- "/nfs/home/students/f.mathis/SysBioMed-PLAs/data/GSM5008737_RNA_3P/"
counts <- Read10X(data.dir = data_path)
counts <- as(counts, "dgCMatrix")

seu <- CreateSeuratObject(
  counts = counts,  
  project = "PlateletEnrichment",
  min.cells = 3,
  min.features = 20 
)

tod <- counts
tod <- tod[Matrix::rowSums(tod) > 0, , drop = FALSE]
tod <- tod[, Matrix::colSums(tod) > 2, drop = FALSE]  

flt <- GetAssayData(seu, slot = "counts")            
flt <- flt[rownames(tod), , drop = FALSE]    

rm(counts)
gc()

# statistics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^HB(?!P)", regex = TRUE)  # fix

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
           percent.hb < 5
)

# cluster for soup
dims.use <- 1:20

seu <- NormalizeData(seu)

# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

hvgs <- devianceFeatureSelection(seu)
VariableFeatures(seu) <- hvgs[1:2000]


seu <- ScaleData(seu) # vars to regress?

seu <- RunPCA(seu)

seu <- FindNeighbors(seu, dims = dims.use)

seu <- FindClusters(seu, resolution = 0.5)

soupx_groups <- seu$seurat_clusters

# TODO annotate platelet clusters

# soup

flt <- flt[Matrix::rowSums(flt) > 0, , drop = FALSE]
tod <- tod[rownames(flt), , drop = FALSE]

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

corrected <- adjustCounts(sc, roundToInt = TRUE)

rm(tod, flt, seu, soupx_groups, keep) 

seu_sx <- CreateSeuratObject(counts = corrected)

keep <- Matrix::rowSums(seu_sx@assays$RNA@counts > 0) >= 20
seu_sx <- seu_sx[keep, ]


gc()

# qc after soup

seu_sx[["percent.mt"]] <- PercentageFeatureSet(seu_sx, pattern = "^MT-")
seu_sx[["percent.ribo"]] <- PercentageFeatureSet(seu_sx, pattern = "^RP[SL]")
seu_sx[["percent.hb"]] <- PercentageFeatureSet(seu_sx, pattern = "^HB(?!P)", regex = TRUE)

seu_sx <- subset(
  seu_sx,
  subset = percent.mt < 20 &
           percent.hb < 5
)

# clustering

seu_sx <- NormalizeData(seu_sx, normalization.method = "LogNormalize",)

# seu_sx <- FindVariableFeatures(seu_sx, selection.method = "vst", nfeatures = 2000)

hvgs <- devianceFeatureSelection(seu_sx)
VariableFeatures(seu_sx) <- hvgs[1:2000]

seu_sx <- ScaleData(seu_sx, vars.to.regress = c("nCount_RNA", "percent.mt")) # todo: andere features hinzufügen?

seu_sx <- RunPCA(seu_sx)

seu_sx <- FindNeighbors(seu_sx, dims = dims.use)

seu_sx <- FindClusters(seu_sx, resolution = 0.5)

# seu_sx <- RunUMAP(seu_sx, dims = dims.use)

# print(DimPlot(seu_sx, group.by = "seurat_clusters", label = TRUE))

# doublet detection
sce <- as.SingleCellExperiment(seu_sx)
sce <- scDblFinder(sce)

seu_sx$scDblFinder_class <- sce$scDblFinder.class
seu_sx$scDblFinder_score <- sce$scDblFinder.score

table(seu_sx$scDblFinder_class)

print(
  DimPlot(seu_sx, group.by = "scDblFinder_class", reduction = "umap")
)

VlnPlot(seu_sx, "scDblFinder_score", group.by = "scDblFinder_class", pt.size = 0.1)

FeaturePlot(seu_sx, "scDblFinder_score", reduction = "umap")

seu_sx <- subset(seu_sx, subset = scDblFinder_class == "singlet")

# final clustering

seu_sx <- NormalizeData(seu_sx, normalization.method = "LogNormalize",)

# seu_sx <- FindVariableFeatures(seu_sx, selection.method = "vst", nfeatures = 2000)

hvgs <- devianceFeatureSelection(seu_sx)
VariableFeatures(seu_sx) <- hvgs[1:2000]

seu_sx <- ScaleData(seu_sx) #, vars.to.regress = c("percent.mt")) # todo: andere features hinzufügen?

seu_sx <- RunPCA(seu_sx)

seu_sx <- FindNeighbors(seu_sx, dims = dims.use)

seu_sx <- FindClusters(seu_sx, resolution = 0.5)

seu_sx <- RunUMAP(seu_sx, dims = dims.use)
saveRDS(seu_sx, file = "data/seu_sx_final.rds")
print(DimPlot(seu_sx, group.by = "seurat_clusters", label = TRUE))

# seu_sx <- readRDS("data/seu_sx_final.rds")
