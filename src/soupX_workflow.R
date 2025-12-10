library(SoupX)

raw_dir <- "./data/GSM5008737_RNA_3P"
filtered_dir <- "./data/filtered_GSM5008737_RNA_3P"

toc <- Seurat::Read10X(filtered_dir)

## 2) SoupChannel nur aus dieser Matrix aufbauen
##    (entspricht dem offiziellen „only counts”-Beispiel in der Doku)
sc <- SoupChannel(toc, toc, calcSoupProfile = FALSE)

## 3) Soup-Profil manuell aus der Gesamt-Expression schätzen
rowSums <- Matrix::rowSums
soupProf <- data.frame(
  row.names = rownames(toc),
  est    = rowSums(toc) / sum(toc),
  counts = rowSums(toc)
)
sc <- setSoupProfile(sc, soupProf)

## 4) Optional: wenn du Clusterlabels hast, hier einlesen und setzen:
## clusters <- read.table("soupx_clusters.tsv", header = TRUE, row.names = 1, sep = "\t")$cluster
## sc <- setClusters(sc, clusters)

## 5) Kontaminationsfraktion automatisch schätzen
sc <- autoEstCont(sc, doPlot = FALSE)

## 6) Korrigierte Counts berechnen (Ganzzahlen für Scanpy)
out <- adjustCounts(sc, roundToInt = TRUE)

## 7) Als 10X-ähnliches Format wieder rausschreiben
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

Matrix::writeMM(out, file.path(output_dir, "matrix.mtx"))
write.table(
  rownames(out),
  file.path(output_dir, "features.tsv"),
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)
write.table(
  colnames(out),
  file.path(output_dir, "barcodes.tsv"),
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)

cat("SoupX correction finished.\n")