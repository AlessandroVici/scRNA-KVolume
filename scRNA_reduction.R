
library(Seurat)
library(ggplot2)
library(pheatmap)
library(factoextra)
library(patchwork)



# Import scRNA-Seq object in R through Seurat package  
# Based on the data stored in GEO (Gene Expression Omnibus, GEO ACCESSION: GSE261385) used in the current work 
# we first download the folder containing the information for each sample used for scRNA-Seq experiment
# The folder is located at the end of the following page (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261385)

# Setting commands
gc()
options(Seurat.object.assay.version = "v5")

# IMPORTANT: Set as working directory the folder where the imported data have been located

# Obtain the single directory in the folder 
dirs <- list.dirs(path = '.', recursive = FALSE, full.names = FALSE)

# Create Seurat object in R for every directory individuated
for(x in dirs){
  name <- x
  
  # Create count matrix (genes x cells)
  cts <- ReadMtx(mtx = paste0('',x,'/matrix.mtx.gz'),
                 features = paste0('',x,'/features.tsv.gz'),
                 cells = paste0('',x,'/barcodes.tsv.gz'))

  # IMPORTANT FILTER: this command is used to eliminate all the genes with zero counts in the current cts matrix
  cts <- cts[rowSums(cts != 0) > 0, ]
  
  # Creating Seurat object with filtered cts matrix
  assign(name, CreateSeuratObject(counts = cts))
}




# Manipulation of Seurat object
# For this work, we focused only on the untreated sample (corresponding to untreated Seurat object in R)

# We will follow the standard Seurat pipeline for data's quality control 
# see https://github.com/quadbio/scRNAseq_analysis_vignette/blob/master/Tutorial.pdf for more informations

untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^mt[-\\.]")
View(untreated@meta.data)

VlnPlot(untreated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(untreated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


plot1 <- FeatureScatter(untreated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(untreated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# First reduction applied for the number of cells, based on information contained in the metadata 
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 1)


# Normalisation phase 

untreated <- NormalizeData(untreated)

# Scaling phase

all_genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all_genes)


# Identify highly-variable features (in our case genes)
# To strongly reduce the variability and potential source of errors we decided to select
# only the top-100 most variable genes

untreated <- FindVariableFeatures(untreated, selection.method =  'vst', nfeatures = 100)
var_genes <- VariableFeatures(untreated)

untreated <- subset(untreated, features = var_genes)

# After quality control phase, we extract from Seurat object "untreated" the dataframe
# containing the normalized/scaled reduced genes expression matrix

untreated_scaled <- data.frame(untreated@assays[["RNA"]]@layers[["scale.data"]])
rownames(untreated_scaled) <- rownames(untreated)
colnames(untreated_scaled) <- colnames(untreated)

# For further analysis we inverted the order of row and columns

untreated_scaled <- data.frame(t(untreated_scaled))
summary(untreated_scaled)

# Graphical correlation for the 100 genes
corr <- round(cor(untreated_scaled),2)
corr <- data.frame(corr)

pheatmap(
  corr,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  labels_row = c(rep(" ",100)),
  labels_col = c(rep(" ",100))
)


# Apply PCA setting 50 as the maximum number of components for the analysis

npcs <- 50 
pca <- prcomp(untreated_scaled,center = T,scale. = T,rank. = npcs)

summary(pca)

# Percentage of Variance explained 

fviz_eig(pca, addlabels = TRUE,ncp=10, main = " ")
