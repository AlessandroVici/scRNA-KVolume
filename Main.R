
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--- LIST OF PACKAGE USED IN THE CODE ---#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#                                        

library(readxl)
library("writexl")
library(gplots)
library(ggplot2)
library(cluster)
library(NetworkToolbox)
library(plotfunctions)
library(RColorBrewer)
library(scales)
library(tidyverse)
library(MASS)
library(readr)
library(clusterProfiler)
library(GOfuncR)
library(dplyr)
library(mclm)
library(plotly)
library(pheatmap)
library(umap)
library(msigdbr)
library(igraph)
library(stringr)
library(extrafont)
library(RColorBrewer)
library(circlize)
library(cluster)
library(factoextra)
library(Seurat)
library(ggplot2)
library(Matrix)
library(org.Mm.eg.db)
library(Mus.musculus)
library(ggpubr)
library(geometry)
library(rgl)
library(plyr)
library(magick)
library(skmeans)
library(Kmedians)
library(mclust)
library(factoextra)
library(sp)
library(clue)
library(celldex)
library(SingleR)
library(installr)
#library(scrapper)
#library(patchwork)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--- LIST OF FUNCTION USED IN THE CODE ---#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  

# K-Area Clustering (2D-Space)
k_center_2d <- function(data, k, seed = 42) {
  set.seed(seed)
  centers <- data[sample(1:nrow(data), 1), , drop = FALSE]
  
  for (i in 2:k) {
    dists <- apply(data, 1, function(p)
      min(apply(centers, 1, function(c) sqrt(sum((p - c)^2)))))
    next_center <- data[which.max(dists), , drop = FALSE]
    centers <- rbind(centers, next_center)
  }
  
  assignments <- apply(data, 1, function(p)
    which.min(apply(centers, 1, function(c) sqrt(sum((p - c)^2)))))
  
  return(list(centers = centers, assignments = assignments))
}

# K-Center clustering function (3D-Space)
k_center <- function(data, k, seed = 42) {
  set.seed(seed)
  centers <- data[sample(1:nrow(data), 1), , drop = FALSE]
  
  for (i in 2:k) {
    dists <- apply(data, 1, function(p)
      min(apply(centers, 1, function(c) sqrt(sum((p - c)^2)))))
    next_center <- data[which.max(dists), , drop = FALSE]
    centers <- rbind(centers, next_center)
  }
  
  assignments <- apply(data, 1, function(p)
    which.min(apply(centers, 1, function(c) sqrt(sum((p - c)^2)))))
  
  return(list(centers = centers, assignments = assignments))
}


# Graphical representation of Convex-Hull
plotbag_allpoints <- function(x, cluster_id, col = "black") {
  x2 <- x[x[, ncol(x)] == cluster_id, ]
  
  if (nrow(x2) >= 4) {
    ids <- t(convhulln(as.matrix(x2[, 1:3]), options = "Tv"))
    triangles3d(x2[ids, 1], x2[ids, 2], x2[ids, 3],
                col = col, alpha = 0.5, shininess = 80)
    spheres3d(x2[ids, 1], x2[ids, 2], x2[ids, 3], r = 0.1, color = col) 
  }
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#---- LIST OF PARAMETERS FOR 3D PLOTS ----#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  


zoom <- 0.9  

aspectr <- c(1, 1, 1)

windowRect <- c(50, 50, 1000, 800)

cex <- 1.2  

pointsalpha <- 0.5  

userMatrix <- rotationMatrix(pi/3, 1, 0, 0) %*% 
  rotationMatrix(pi/6, 0, 1, 0) %*% 
  rotationMatrix(0, 0, 0, 1) %*%
  rotationMatrix(-pi/10, 0, 1, 0)  

# P.S. If you should encounter problems with 3D visualization through the use of plot3D function,
# change the inside command "type = s" with a different type (for example " type = p").
# This issue probably depends on the version of R used to execute the code.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------- SIMULATED CLUSTERING CODE -------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 


#Set seed to reproduce results
set.seed(123)

# Create three almost-spherical clusters
cluster1 <- cbind(runif(10, min = 0, max = 0.5), runif(10, min = -0.2, max = 0.7))
cluster2 <- cbind(runif(10, min = 1, max = 2), runif(10, min = 1, max = 1.5))
cluster3 <- cbind(runif(10, min = 2.5, max = 3.5), runif(10, min = 0, max = 1))

# Create a collinear cluster
cluster4 <- cbind(c(1.3, 1.5, 1.7, 1.9, 2.1), c(0.1, 0.2, 0.4, 0.6, 0.7))

# Dataset containing the four clusters
data <- rbind(cluster1, cluster2, cluster3, cluster4)

# df is the dataframe that will contain different assignments for each clustering type
df <- as.data.frame(data)
df$Cluster <- factor(c(rep("1", 10), rep("2", 10), rep("3", 10), rep("4", 5))) # original assignment


# Set different colors for each cluster
colors <- c("red", "skyblue1", "green3", "violet")
names(colors) <- levels(df$Cluster)


par(mgp = c(2, 0.4, 0))  #optimal graphical settings


plot(df$V1, df$V2, pch = 21, bg = colors[df$Cluster], col = colors[df$Cluster], cex = 1.2,
     xlab = "V1", ylab = "V2", main = "Projected Data into 2D Space",xaxt = "n", yaxt = "n")

# Custom x-axis
axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
     tck = -0.02, mgp = c(0, 0.4, 0))

# Custom y-axis
axis(2, at = c(0, 0.5, 1, 1.5),
     labels = c(0, 0.5, 1, 1.5),
     tck = -0.02, mgp = c(0, 0.4, 0))

legend("topleft", legend = paste("Cluster", levels(df$Cluster)),
       pt.bg = colors, pch = 21, pt.cex = 1.2, col = colors)

#lims <- par("usr")

# Add letter at the bottom of the plot
#text(x = lims[1] + 1.8, y = lims[3] - 0.5, labels = "(a)", xpd = NA, cex = 1, font = 2, family = "Palatino Linotype")




#~~~ Application of K-Means clustering for K=4 clusters ~~~

km <- kmeans(df[, 1:2], centers = 4)
df$km_cluster <- factor(km$cluster) #K-Means assignment

# Indexes calculation

# Adjusted Rand Index (ARI)
ari_kmeans <- adjustedRandIndex(df$Cluster, df$km_cluster)

# Silhouette Score
sil <- silhouette(as.numeric(df$km_cluster), dist(df[,1:2]))
silhouette_kmeans <- mean(sil[,3])


cat("Adjusted Rand Index (K-Means):", round(ari_kmeans, 4), "\n")
cat("Silhouette Score (K-Means):", round(silhouette_kmeans, 4), "\n")

# 2D Lacunarity Index Calculation Procedure
total_area <- convhulln(df[,1:2], options = "FA")$vol
cat("Total area of the dataset:", round(total_area, 4), "\n")

# Area calculation for each cluster
cluster_areas <- numeric()
for (k in levels(df$km_cluster)) {
  group_data <- df[df$km_cluster == k, 1:2]
  if (nrow(group_data) >= 3) {
    area <- convhulln(as.matrix(group_data), options = "FA")$vol
  } else {
    area <- 0
  }
  cluster_areas <- c(cluster_areas, area)
  cat("Cluster", k, "Area:", round(area, 4), "\n")
}

# Lacunarity score
sum_cluster_areas <- sum(cluster_areas)
lacunarita <- (total_area - sum_cluster_areas) / total_area
cat("Sum clusters'areas:", round(sum_cluster_areas, 4), "\n")
cat("Lacunarity", round(lacunarita, 4), "\n")



# Set different colors for each cluster
colors <- c("red", "skyblue1", "green3", "violet")
names(colors) <- levels(df$km_cluster)


par(mgp = c(2, 0.4, 0))


plot(df$V1, df$V2, pch = 21, bg = colors[df$km_cluster], col = colors[df$km_cluster], cex = 1.2,
     xlab = "V1", ylab = "V2", main = "K-Means Clustering with Convex Hulls",xaxt = "n", yaxt = "n")

# Custom x-axis
axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
     tck = -0.02, mgp = c(0, 0.4, 0))

# Custom y-axis
axis(2, at = c(0, 0.5, 1, 1.5),
     labels = c(0, 0.5, 1, 1.5),
     tck = -0.02, mgp = c(0, 0.4, 0))

# Add centroids for each cluster
points(km$centers[,1], km$centers[,2], pch = 3, cex = 2, lwd = 2)

# Add Convex-Hulls to delimitate clusters
for (i in levels(df$km_cluster)) {
  cluster_points <- df[df$km_cluster == i, 1:2]
  if (nrow(cluster_points) >= 3) {
    hull_indices <- chull(cluster_points)
    polygon(cluster_points[hull_indices, ], border = colors[i], lwd = 2)
  }
}

legend("topleft", legend = paste("Cluster", levels(df$km_cluster)),
       pt.bg = colors, pch = 21, pt.cex = 1.2, col = colors)

#lims <- par("usr")

# Add letter at the bottom of the plot
#text(x = lims[1] + 1.8, y = lims[3] - 0.5, labels = "(b)", xpd = NA, cex = 1, font = 2, family = "Palatino Linotype")




#~~~~ Application of K-Area clustering for K=4 clusters ~~~


# Clustering
k <- 4
clustering <- k_center_2d(df[,1:2], k)
centers <- clustering$centers

# Create contingency matrix
contingency <- table(df$Cluster, clustering$assignments)

# Apply optimal assignment algorithm 
mapping <- solve_LSAP(contingency, maximum = TRUE)

# New alligned labels 
aligned_assignments <- mapping[clustering$assignments]

# Original labels vs alligned labels
table(Original = df$Cluster, Aligned = aligned_assignments)

df$Cluster_K_vol <- factor(aligned_assignments)


# Adjusted Rand Index (ARI)
ari_kmeans1 <- adjustedRandIndex(df$Cluster, aligned_assignments)

# Silhouette Score
sil1 <- silhouette(aligned_assignments, dist(df[,1:2]))
silhouette_kmeans1 <- mean(sil1[,3])

cat("Adjusted Rand Index (K-Area):", round(ari_kmeans1, 4), "\n")
cat("Silhouette Score (K-Area):", round(silhouette_kmeans1, 4), "\n")


total_area <- convhulln(df[,1:2], options = "FA")$vol
cat("Total area of the dataset:", round(total_area, 4), "\n")


cluster_areas <- numeric()
for (k in levels(df$Cluster_K_vol)) {
  group_data <- df[df$Cluster_K_vol == k, 1:2]
  if (nrow(group_data) >= 3) {
    area <- convhulln(as.matrix(group_data), options = "FA")$vol
  } else {
    area <- 0
  }
  cluster_areas <- c(cluster_areas, area)
  cat("Cluster", k, " Area:", round(area, 4), "\n")
}

# Lacunarity score
sum_cluster_areas <- sum(cluster_areas)
lacunarita <- (total_area - sum_cluster_areas) / total_area
cat("Sum clusters'areas:", round(sum_cluster_areas, 4), "\n")
cat("Lacunarity:", round(lacunarita, 4), "\n")




# Set different colors for each cluster
colors <- c("red", "skyblue1", "green3", "violet")
names(colors) <- levels(df$Cluster_K_vol)


par(mgp = c(2, 0.4, 0))

plot(df$V1, df$V2, pch = 21, bg = colors[df$Cluster_K_vol], col = colors[df$Cluster_K_vol], cex = 1.2,
     xlab = "V1", ylab = "V2", main = "K-Area Clustering with Convex Hulls",xaxt = "n", yaxt = "n")


axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
     tck = -0.02, mgp = c(0, 0.4, 0))


axis(2, at = c(0, 0.5, 1, 1.5),
     labels = c(0, 0.5, 1, 1.5),
     tck = -0.02, mgp = c(0, 0.4, 0))


for (cl in levels(df$Cluster_K_vol)) {
  pts <- df[df$Cluster_K_vol == cl, c("V1", "V2")]
  if (nrow(pts) >= 3) {
    hull <- chull(pts)
    polygon(pts[hull, ], border = colors[cl], lwd = 2)
  }
}


points(centers[,1], centers[,2], pch = 3, cex = 2, lwd = 2)


legend("topleft", legend = paste("Cluster", levels(df$Cluster_K_vol)),
       pt.bg = colors, pch = 21, pt.cex = 1.2, col = colors)

#lims <- par("usr")

# Add letter at the bottom of the plot
#text(x = lims[1] + 1.8, y = lims[3] - 0.5, labels = "(c)", xpd = NA, cex = 1, font = 2, family = "Palatino Linotype")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------- REAL SCRNA DATASET CLUSTERING CODE -------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#~~~ Creation of Seurat object ~~~                            

# Import scRNA-Seq object in R through Seurat package  
# Based on the data stored in GEO (Gene Expression Omnibus, GEO ACCESSION: GSE261385) used in the current work 
# we first download the folder containing the information for each sample used for scRNA-Seq experiment
# The folder is located at the end of the following page (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261385)

# Setting commands
gc()
options(Seurat.object.assay.version = "v5")

# IMPORTANT: Set as working directory the folder where the imported data have been located
# For our case, the data are stored in this working directory:

setwd("C:/Users/vici_alessandro/OneDrive - Istituto Superiore di Sanità/ISS/Lavori/Lavoro_AG_K_volume/GSE261385_processed")

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



#~~~ Manipulation of Seurat object ~~~

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

#PCA taking the first 20 PCs
scRNA_all <- RunPCA(untreated, features = var_genes,npcs=20)

# Standard devations of PCs
std_devs <- scRNA_all@reductions[["pca"]]@stdev

# Variance of each components
variances <- std_devs^2

# % of explained Variance
explained_var <- variances / sum(variances) * 100

explained_var/100

pc_labels <- paste0("PC", 1:length(explained_var))


barplot(explained_var[1:10],
        names.arg = pc_labels[1:10],
        col = "skyblue",
        border = "black",
        main = " ",
        xlab = " ",
        ylab = "Percentage of explained variance",
        ylim = c(0, 30),las=2)







#We extract the scores for the cells

pca_scores <- data.frame(scRNA_all@reductions[["pca"]]@cell.embeddings)


df <- pca_scores 
df <- df[, c(1, 2, 3)]  #first 3 PCs
data <- as.matrix(df)


#A first plot of cells' score coordinates
open3d(windowRect = windowRect)

par3d(zoom = zoom, userMatrix = userMatrix)

plot3d(data[, 1], data[, 2], data[, 3],
       aspect = aspectr,
       col = "brown", size = 1, type='s',
       main = " ",
       xlab = "PC1",
       ylab = "PC2",
       zlab = "PC3"
)




#~~~ Clustering Phase procedure ~~~ 
   
#Let's start with K-volume clustering
#k_vol_df will be the dataset used for further analysis with K-volume procedure

k_vol_df <- df

# Total volume of the points contained in the dataset
total_volume <- convhulln(data, options = "FA")$vol
cat("Total volume of the dataset:", total_volume, "\n")

# Vector created to store the value of volume for each k clustering solution proposed
cluster_volume_sums <- numeric()

# We will test different k values (from 3 to 10)
for (k in 3:10) {
  clustering <- k_center(data, k)
  centers <- clustering$centers
  assignments <- clustering$assignments
  df$Cluster <- factor(assignments)
  
  cat("\n--- k =", k, "---\n")
  cluster_volumes <- numeric(k)
  
  for (i in 1:k) {
    group_data <- df[df$Cluster == i, 1:3]
    
    if (nrow(group_data) >= 4) {
      vol <- convhulln(as.matrix(group_data), options = "FA")$vol
    } else {
      vol <- 0
    }
    
    cluster_volumes[i] <- vol
    cat("Volume cluster", i, ":", vol, "\n")
  }
  
  cluster_volume_sums[k] <- sum(cluster_volumes)
  cat("Sum of total volumes for k =", k, ":", cluster_volume_sums[k], "\n")
  cat("Silhouette score...", "\n")
  
  sil <- silhouette(as.numeric(df$Cluster), dist(df[,1:3]))
  silhouette_means <- mean(sil[,3])
  cat("Silhouette Score for k =", k, ":", silhouette_means, "\n")
}


cluster_volume_sums <- na.omit(cluster_volume_sums) #remove NA values
lacunarity_k <- (total_volume - cluster_volume_sums) / total_volume
cat("Lacunarity for k =", k, ":", lacunarity_k, "\n")


par(mgp = c(2, 0.4, 0))

plot(c(3:10), lacunarity_k, type = "b",
     xlab = "Number of clusters (K)",
     ylab = "Lacunarity index",
     main = " ",
     xaxt = "n", yaxt = "n",
     panel.first = rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
                        col = "white", border = NA),
     pch = 16, col = "grey30")

axis(1, at = c(3, 5, 7, 9), labels = c(3, 5, 7, 9),
     tck = -0.02, mgp = c(0, 0.4, 0))

axis(2, at = c(0.90, 0.92, 0.94, 0.96, 0.98),
     labels = c("0.90", "0.92", "0.94", "0.96", "0.98"),
     tck = -0.02, mgp = c(0, 0.4, 0))

#Optimal value
points(5, lacunarity_k[3], pch = 1, col = "red", cex = 3, lwd = 2)


#We represent the 3D plot for the optimal k value selected previously

k <- 5
clustering <- k_center(data, k)
centers <- clustering$centers
assignments <- clustering$assignments
k_vol_df$Cluster <- factor(assignments)


cluster <- k_vol_df$Cluster  

palette <- c("red", "lightblue1", "lightgreen", "orange", "purple1")
colors <- c()

# Code for the correspondence between colors and clusters

for(i in 1:length(cluster)){
  colors[i] <- ifelse(cluster[i] == 1,"red",
                      ifelse(cluster[i] == 2,"lightblue1",
                             ifelse(cluster[i] == 3,"lightgreen",
                                    ifelse(cluster[i] == 4,"orange",
                                           ifelse(cluster[i] == 5,"purple1")))))
}


colors_transparent <- adjustcolor(colors, alpha.f = pointsalpha) #graphical optimization

open3d(windowRect = windowRect)

par3d(zoom = zoom, userMatrix = userMatrix)

plot3d(k_vol_df[, 1], k_vol_df[, 2], k_vol_df[, 3],
       aspect = aspectr,
       col = colors_transparent, size = 1, type='s',
       main = " ",
       xlab = "PC1",
       ylab = "PC2",
       zlab = "PC3"
)

legend3d(x=.30, y=.75, legend = paste("Cluster", as.integer(levels(k_vol_df$Cluster))),
         col = palette, pch = 16, cex = 1.2, inset = c(0.02))








#K-means clustering
#k_means_df will be the dataset used for further analysis with K-means procedure

k_means_df <- df

total_volume <- convhulln(data, options = "FA")$vol
cat("Total volume:", total_volume, "\n")


cluster_volume_sums <- numeric()


for (k in 3:10) {
  set.seed(1234)
  clustering <- kmeans(data, k)
  centers <- clustering$centers
  assignments <- clustering$cluster
  df$Cluster <- factor(assignments)
  
  cat("\n--- k =", k, "---\n")
  cluster_volumes <- numeric(k)
  
  for (i in 1:k) {
    group_data <- df[df$Cluster == i, 1:3]
    
    if (nrow(group_data) >= 4) {
      vol <- convhulln(as.matrix(group_data), options = "FA")$vol
    } else {
      vol <- 0
    }
    
    cluster_volumes[i] <- vol
    cat("Cluster volume", i, ":", vol, "\n")
  }
  
  cluster_volume_sums[k] <- sum(cluster_volumes)
  cat("Sum of total volumes for k =", k, ":", cluster_volume_sums[k], "\n")
  cat("Silhouette score...", "\n")

  sil <- silhouette(as.numeric(df$Cluster), dist(df[,1:3]))
  silhouette_means <- mean(sil[,3])
  cat("Silhouette Score for k =", k, ":", silhouette_means, "\n")
}


cluster_volume_sums <- na.omit(cluster_volume_sums)
lacunarity_k <- (total_volume - cluster_volume_sums) / total_volume
cat("Lacunarity per k =", k, ":", lacunarity_k, "\n")


par(mgp = c(2, 0.4, 0))

plot(c(3:10), lacunarity_k, type = "b",
     xlab = "Number of clusters (K)",
     ylab = "Lacunarity index",
     main = " ",
     xaxt = "n", yaxt = "n",
     panel.first = rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
                        col = "white", border = NA),
     pch = 16, col = "grey30")


axis(1, at = c(3, 5, 7, 9), labels = c(3, 5, 7, 9),
     tck = -0.02, mgp = c(0, 0.4, 0))


axis(2, at = c(0.90, 0.92, 0.94, 0.96, 0.98),
     labels = c("0.90", "0.92", "0.94", "0.96", "0.98"),
     tck = -0.02, mgp = c(0, 0.4, 0))


points(5, lacunarity_k[3], pch = 1, col = "red", cex = 3, lwd = 2)



# Optimal k value is k=5
k <- 5
clustering <- kmeans(data, k)
centers <- clustering$centers
assignments <- clustering$cluster
k_means_df$Cluster <- factor(assignments)


cluster <- k_means_df$Cluster  


palette <- c("red", "lightblue1", "lightgreen", "orange", "purple1")
colors <- c()


for(i in 1:length(cluster)){
  colors[i] <- ifelse(cluster[i] == 1,"red",
                      ifelse(cluster[i] == 2,"lightblue1",
                             ifelse(cluster[i] == 3,"lightgreen",
                                    ifelse(cluster[i] == 4,"orange",
                                           ifelse(cluster[i] == 5,"purple1")))))
}

colors_transparent <- adjustcolor(colors, alpha.f = pointsalpha)


open3d(windowRect = windowRect)

par3d(zoom = zoom, userMatrix = userMatrix)


plot3d(k_means_df[, 1], k_means_df[, 2], k_means_df[, 3],
       aspect = aspectr,
       col = colors_transparent, size = 1, type='s',
       main = " ",
       xlab = "PC1",
       ylab = "PC2",
       zlab = "PC3"
)

legend3d(x=.30, y=.75, legend = paste("Cluster", as.integer(levels(k_means_df$Cluster))),
         col = palette, pch = 16, cex = 1.2, inset = c(0.02))




#Louvain clustering
#Louvain_df will be the dataset used for further analysis with Louvain procedure


# Clustering ============================
scRNA_all <- FindNeighbors(scRNA_all, dims = 1:3)
scRNA_all <- FindClusters(scRNA_all, resolution = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5))

metadata <- data.frame(scRNA_all@meta.data)
metadata$Cells <- rownames(metadata)

metadata <- metadata[,-c(1:4,13)]
pca_scores$Cells <- rownames(pca_scores)
metadata <- merge(metadata,pca_scores[,c(1:3,21)],by="Cells")
metadata <- metadata[,c(1,10:12,2:9)]



total_volume <- convhulln(data, options = "FA")$vol
cat("Total volume of the dataset:", total_volume, "\n")

cluster_volume_sums <- numeric()


for (k in 1:8) {
  set.seed(1234)
  
  n_group <- length(levels(metadata[,k+4]))
  # Calcolo volumi dei cluster
  cat("\n--- k =", n_group, "per ", colnames(metadata)[k+4],"---\n")
  cluster_volumes <- numeric(n_group)
  
  for (i in 0:n_group-1) {
    group_data <- metadata[metadata[,k+4] == i, 2:4]
    
    if (nrow(group_data) >= 4) {
      vol <- convhulln(as.matrix(group_data), options = "FA")$vol
    } else {
      vol <- 0
    }
    
    cluster_volumes[i+1] <- vol
    cat("Volume cluster", i, ":", vol, "\n")
  }
  
  cluster_volume_sums[k] <- sum(cluster_volumes)
  cat("Sum of clusters' volume for k =", k, ":", cluster_volume_sums[k], "\n")
  cat("Silhouette score...", "\n")
  
  sil <- silhouette(as.numeric(metadata[,k+4]), dist(metadata[,2:4]))
  silhouette_means <- mean(sil[,3])
  cat("Silhouette Score for k =", n_group, ":", silhouette_means, "\n")
}


cluster_volume_sums <- na.omit(cluster_volume_sums)
lacunarity_k <- (total_volume - cluster_volume_sums) / total_volume
cat("Lacunarity for k =", k, ":", lacunarity_k, "\n")



assignments <- as.numeric(as.character(metadata[,7]))
names(assignments) <- metadata$Cells

#Update of the original dataset after the clustering optimal solution 
Louvain_df <- cbind(metadata[,c(2:4)],assignments)

cluster <- Louvain_df$assignments  

palette <- c("yellow1", "red", "lightblue1", "lightgreen", "orange", "purple1")
colors <- c()

for(i in 1:length(cluster)){
  colors[i] <- ifelse(cluster[i] == 0, "yellow1",ifelse(cluster[i] == 1,"red",
                      ifelse(cluster[i] == 2,"lightblue1",
                             ifelse(cluster[i] == 3,"lightgreen",
                                    ifelse(cluster[i] == 4,"orange",
                                           ifelse(cluster[i] == 5,"purple1"))))))
}

colors_transparent <- adjustcolor(colors, alpha.f = pointsalpha)


open3d(windowRect = windowRect)

par3d(zoom = zoom, userMatrix = userMatrix)

plot3d(Louvain_df[, 1], Louvain_df[, 2], Louvain_df[, 3],
       aspect = aspectr,
       col = colors_transparent, size = 1, type='s',
       main = " ",
       xlab = "PC1",
       ylab = "PC2",
       zlab = "PC3"
)

legend3d(x=.30, y=.75, legend = paste("Cluster", as.integer(levels((as.factor((Louvain_df$assignments)))))),
         col = palette, pch = 16, cex = 1.2, inset = c(0.02))









#~~~ Minor PC Evaluation ~~~ 


#We assess the statistical significance of the minor components (distance from the standard normal in terms of mean 
#zero and standard deviation one). In practice, for very small clusters, we apply Fisher’s exact test to evaluate 
#the null hypothesis that 5% of observations fall outside the two-standard-error interval on either side, based on 
#clusters obtained using the major components. Since the minor components were NOT used for clustering, this is 
#effectively an 'external validation' test of the relevance of the minor components.

#We use Fisher’s exact test if the cluster has fewer than 30 observations; otherwise, we use the Chi-square test.
#When dealing with a small cluster and wanting to verify whether its observations show anomalous behavior along a 
#specific principal component—such as the sixth one—Fisher’s exact test can be applied.

#To do this, we start by defining what is meant by “anomaly”: in this context, a value is considered anomalous 
#if it lies outside the interval between two standard deviations above or below the mean, as expected from a 
#standard normal distribution.

#Once this criterion is established, we count the observations in the cluster that fall outside this interval 
#and those that fall within it. These numbers are then compared to what would be expected if the data truly 
#followed a standard normal distribution, in which only a small percentage of observations should lie beyond 
#two standard deviations.

#At this point, a contingency table is constructed, where one row represents the observed data in the cluster 
#and the other row represents the expected values according to the theoretical distribution. The first column 
#of the table contains the number of observations outside the interval, while the second column contains those within it.

#Fisher’s exact test is then applied to this table to calculate the probability of observing such a marked deviation—or 
#an even more extreme one—purely by chance. If this probability, i.e., the p-value, is sufficiently low, we can conclude 
#that the cluster exhibits anomalous behavior along that principal component, suggesting that it may contain 
#relevant information not captured by the major components.



# We create a single dataset containing the first 3 PCs and one column with the clustering results
# for each method used (thus, 3 additional columns)

colnames(k_vol_df)[4] <- "Cluster_Kvol"
colnames(k_means_df)[4] <- "Cluster_Kmeans"
colnames(Louvain_df)[4] <- "Cluster_Louvain"

total <- cbind(k_vol_df,k_means_df$Cluster_Kmeans,Louvain_df$Cluster_Louvain)
colnames(total)[c(5,6)] <- c("Cluster_Kmeans","Cluster_Louvain")

# We use the various datasets obtained one at a time

total$Cells <- rownames(total)

# We remove the first 3 PCs which we already have
pca_scores <- pca_scores[,-c(1,2,3)]

total <- merge(total,pca_scores,by="Cells")
rownames(total) <- total$Cells
total <- total[,-c(1,15:24)]
total <- total[,c(4:6,1:3,7:13)]

total[,c(4:13)] <- round(total[,c(4:13)],3)

total$Cluster_Louvain <- as.factor(total$Cluster_Louvain)

data_list <- list()
methods <- c("K-volume","K-means","Louvain")

for(m in 1:3){
  # Main loop that iterates over each method used (3 in this case)
  
  cat("Current method...", methods[m], "\n")
  cat("\n")
  cat("\n")
  
  n_clusters <- as.numeric(table(total[,m]))
  
  # # Sort the dataset by Cluster (K-means for example)
  # total <- total %>%
  #   arrange(total[,m])
  
  # We create the matrix that will contain the p-values and will be inserted
  # into the "data_list" list
  
  matrix_pval <- data.frame(matrix(NA,nrow=length(n_clusters), ncol=8))
  colnames(matrix_pval) <- c("PC4","PC5","PC6","PC7","PC8","PC9","PC10","Num")
  rownames(matrix_pval) <- as.character(levels(total[,m]))
  
  # Now we move on to the actual algorithm
  
  for(i in 1:length(n_clusters)){
    
    current_cluster_label <- levels(total[,m])[i]
    cat("Current cluster: ", current_cluster_label, "that contains ", n_clusters[i], "cells", "\n")
    
    matrix_pval[i,8] <- n_clusters[i]
    
    # We create and fill the 2x2 table to be used for the subsequent statistical tests
    
    for(j in 7:13){
      
      # j handles the various PCs
      
      cat("Cluster ", current_cluster_label, "related to ", colnames(total)[j], "\n")  
      
      matrix_test <- data.frame(matrix(NA,nrow=2,ncol=2))
      rownames(matrix_test) <- c("Current","Expected")
      colnames(matrix_test) <- c("Out ±2 SD","In ±2 SD")
      
      # Filling the second row (expected counts)...
      
      matrix_test[2,1] <- round((n_clusters[i]/100)*5,0)
      matrix_test[2,2] <- round((n_clusters[i]/100)*95,0)
      
      # To fill the first row instead...
      
      current <- total[total[,m] == current_cluster_label,]
      
      matrix_test[1,1] <- sum(current[,j] >= 2 | current[,j] <= -2)
      matrix_test[1,2] <- n_clusters[i] - matrix_test[1,1]
      
      if(n_clusters[i] <= 30){
        
        # Perform the test
        risultato <- fisher.test(matrix_test)
        
        # Extract the p-value
        p <- risultato$p.value
        
        matrix_pval[i,j-6] <- p
        
        cat("Fisher’s exact test between ", current_cluster_label, "and ", colnames(total)[j], "has p-value: ", p , "\n") 
      }
      else{
        # Perform the test
        risultato <- chisq.test(matrix_test,simulate.p.value = T)
        
        # Extract the p-value
        p <- risultato$p.value
        
        matrix_pval[i,j-6] <- p
        
        cat("Chi-square test between ", current_cluster_label, "and ", colnames(total)[j], "has p-value: ", p , "\n") 
      }
    }
  }
  
  data_list[[m]] <- matrix_pval
  
}

round(data_list[[1]],3)
round(data_list[[2]],3)
round(data_list[[3]],3)

# For K-Volume, cluster C4 stands out with respect to PC7
# For K-means, cluster C4 stands out with PC7
# For Louvain, cluster C3 stands out with PC7


a <- total[total$Cluster_Kvol == 4,c(1,4:13)]
b <- total[total$Cluster_Kmeans == 4,c(1,4:13)]
c <- total[total$Cluster_Louvain == 3,c(1,4:13)]





#~~~ Biological characterization ~~~ 

# Use of SingleR package
set.seed(42)

surveyReferences()


# Get reference dataset
#This celldex dataset has label.main, which has rougher or more general cell type labels, and label fine. 
#We’ll go with label.main for general purpose

ref <- celldex::ImmGenData()
unique(ref$label.main)


# SingleR can use raw or normalised counts from our dataset

norm_counts <- LayerData(scRNA_all, assay = "RNA", layer = 'data') 


# SingleR can take many arguments, but the main ones are:
  
# test, our normalised/raw counts.
# ref, our reference counts, normalised (you can give it just the matrix of log transformed expression values, 
# but since I have the summarizedexperiment format already, I’m just going to give it this).
# labels, the labels of each of the cells in our reference dataset.
# de.method, the differential expression method it uses to predict the cell types. 
# For single-cell data, the authors recommend the Wilcoxon ranked test.


# Run SingleR
ct_ann <- SingleR(test = norm_counts, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.main,
                  de.method = 'wilcox')

# The function returns a dataframe. We have our cell IDs, the predicted label, and then some statistics on that prediction.

# What about delta and pruned labels?
  
# Well, by itself, the SingleR algorithm will always assign a label to every cell. 
# But what if the cell’s true label isn’t in the reference dataset? It will still assign it a label. 
# For this and other reasons, SingleR can assign incorrect labels.

# The developers tried to mitigate this  by removing poor-quality predictions with “low” scores. 
# They basically compute a “delta” value for each cell, which is the gap, or the difference between the 
# score for the assigned label and the median score across all labels. If the delta is small, 
# this indicates that the cell matches all labels with the same confidence, so the assigned label is 
# not very meaningful.

# This way, SingleR can discard cells with low delta values caused by (i) ambiguous assignments with 
# closely related reference labels and (ii) incorrect assignments that match poorly to all 
# reference labels – so in the pruned_labels column you will find ‘cleaner’ or more reliable labels.


unique(ct_ann$pruned.labels)
table(ct_ann$pruned.labels)
table(ct_ann$labels)
summary(is.na(ct_ann$pruned.labels))


# Inspect quality of the predictions
plotScoreHeatmap(ct_ann)

# Check the type of cells contained in the clusters previously individuated 

labels <- ct_ann$pruned.labels
cells <- rownames(ct_ann)
result <- data.frame(cbind(cells,labels)) 

result[result$cells %in% rownames(a),] 
result[result$cells %in% rownames(b),] 
result[result$cells %in% rownames(c),] 

#It seems that each cluster related to PC7 are composed of endothelial cells


scRNA_all <- AddMetaData(scRNA_all, ct_ann$pruned.labels, col.name = 'SingleR_HCA')
meta <- scRNA_all@meta.data
meta$Cells <- rownames(meta)
total$Cells <- rownames(total)

meta <- merge(meta,total[,c(14,1:3)],by="Cells")



#Let's continue the biological characterization usign Gene Ontology
#First we select the most important genes for PC7

VizDimLoadings(scRNA_all, dims = 7, reduction = "pca",nfeatures = 50)

DimHeatmap(scRNA_all, dims = 7, nfeatures=100, cells = 321, balanced = TRUE)

#To work on the genes we have to take the loadings of the PCA
loadings <- round(data.frame(scRNA_all@reductions[["pca"]]@feature.loadings),3)

#Genes of interest
loadings$Gene <- rownames(loadings)

#We select only the genes that have correlation values above 0.1 and below -0.1
genes <- as.character(loadings[loadings$PC_7 <= -0.1 | loadings$PC_7 >= 0.1,21])


#Gene Ontology procedure
Input_Genes <- genes
Input_Genes <- data.frame(Input_Genes)


Input_Genes$Candidates <- 1
Go_Enrich_Out <- go_enrich(Input_Genes,organismDb = 'Mus.musculus')
Results <- Go_Enrich_Out$results

Over_Representation <- Results[Results$raw_p_overrep<=0.05,]
Under_Representation <- Results[Results$raw_p_underrep<=0.05,]

Genes <- Go_Enrich_Out$genes
Candidate_Gene <- Genes[Genes$Candidates==1,]
Gene_all_GO <- get_anno_categories(Candidate_Gene[,1],database = 'Mus.musculus')

Out_Over_Representation <- merge(Over_Representation,Gene_all_GO, by.x = "node_id",by.y = "go_id", all.x = TRUE, all.y = FALSE)
Out_Under_Representation <- merge(Under_Representation,Gene_all_GO, by.x = "node_id",by.y = "go_id", all.x = TRUE, all.y = FALSE)


Results_Over_Representation <- Out_Over_Representation %>% group_by(node_id) %>% mutate(gene= paste(gene, collapse=",")) %>% unique %>% na.omit
Results_Under_Representation <- Out_Under_Representation %>% group_by(node_id) %>% mutate(gene= paste(gene, collapse=",")) %>% unique %>% na.omit
Results_Over_Representation[ ,c(9,10)] <- NULL 
Results_Under_Representation[ ,c(9,10)] <- NULL


Results_Over_Representation$ontology <- as.factor(Results_Over_Representation$ontology)
table(Results_Over_Representation$ontology)
barplot(table(Results_Over_Representation$ontology))

Results_Over_Representation <- Results_Over_Representation %>%
  arrange(ontology, raw_p_overrep)


# We show the top 10 most significant pathways for each category (Biological Process, Molecular Function, Cellular Component)
top_pathways <- Results_Over_Representation %>%
  group_by(ontology) %>%
  top_n(-10, raw_p_overrep)

top_pathways$log10_p <- -log10(top_pathways$raw_p_overrep)

top_pathways <- data.frame(top_pathways)


# df1 will be the reduced dataset with Gene Ontology results
df1 <- data.frame(
  Term = top_pathways$node_name,
  Value = top_pathways$log10_p,
  Category = top_pathways$ontology)

df1$Category <- factor(df1$Category, levels = unique(df1$Category))
df1[4,1] <- c("positive regulation of chemokine CXCL2 production")
df1$Term <- factor(df1$Term, levels = df1$Term)
df1$Category <- factor(df1$Category, levels = c("biological_process", "cellular_component", "molecular_function"))



setwd("C:/Users/vici_alessandro/OneDrive - Istituto Superiore di Sanità/ISS/Lavori/Lavoro_AG_K_volume/Figures")
tiff('Figure7c.tiff', units="in", width=9.5, height=10.5, res=300, compression = "lzw", bg = "white")  # Sfondo esterno colorato

#Bar plot for Gene Ontology 
ggplot(df1, aes(x = Value, y = Term)) + 
  geom_bar(stat='identity', aes(fill=Category), width=.5, show.legend = TRUE) + 
  labs(title= " ") + 
  labs(y = " ", x = expression(-log[10]~"(" * italic(p) * "-value" * ")")) + 
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA), 
    panel.border = element_rect(color = "black", fill = NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),  
    axis.text.y = element_text(size = 16, color = "black"),
    legend.position = c(0.7, 0.50),  
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10)) +  
  scale_x_continuous(limits = c(0, 9)) +  
  scale_fill_manual(values = c("biological_process" = "darkgreen", "cellular_component" = "blue", "molecular_function" = "red"),  
                    labels = c("Biological Process", "Cellular Component", "Molecular Function"))  


dev.off()








#~~~ UMAP approach using Seurat ~~~

# scRNAall is the seurat object that we have used before to analyze single cell RNA-Seq dataset.
# The most commonly used non-linear dimension 
# reduction methods in scRNA-seq data analysis are t-distributed Stochastic Neighbor Embedding (t-SNE) and Uniform Manifold 
# Approximation and Projection (UMAP). Both methods try to place every sample in a low-dimensional space (2D/3D), 
# so that distances or neighborhood relationships between different samples (here cells) in the original space are 
# largely retained in the low-dimensional space. The detailed mathematically descriptions of the two methods are out 
# of the scope of this tutorial, but for those who are interested in, you may check this video for tSNE, and this 
# blog of Nikolay Oskolkov for UMAP. There are also more methods to create other low-dimensional embeddings for 
# visualization, including but not limiting to SPRING, PHATE. Now let's focus on tSNE and UMAP which Seurat has included. 
# The top PCs in the PCA analysis are used as the input to create a tSNE and UMAP embedding of the data.

scRNA_all <- FindNeighbors(scRNA_all, dims = 1:3)
scRNA_all <- FindClusters(scRNA_all, resolution = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5))
DimPlot(scRNA_all, group.by = 'RNA_snn_res.0.3', label = TRUE)
Idents(scRNA_all) <- 'RNA_snn_res.0.3' # set identity of clusters





scRNA_all <- RunUMAP(scRNA_all, dims = 1:3,n.components = 3)  #based on 3 first PC (as PCA before)
DimPlot(scRNA_all, reduction = 'umap')


umap <- data.frame(scRNA_all@reductions[["umap"]]@cell.embeddings)
umap$Cells <- rownames(umap)
umap <- merge(umap,meta[,c(1,15:18)],by="Cells")


open3d(windowRect = windowRect)

par3d(zoom = zoom, userMatrix = userMatrix)

plot3d(umap[, 2], umap[, 3], umap[, 4],
       aspect = aspectr,
       col = "brown", size = 1, type='s',
       main = " ",
       xlab = "UMAP_1",
       ylab = "UMAP_2",
       zlab = "UMAP_3"
)




# Now we assign different colors to the clusters 
cluster <- as.numeric(meta$RNA_snn_res.0.3)  # oppure un vettore separato


palette <- c("red", "lightblue1", "lightgreen", "orange", "purple1")
colors <- c()

# Assignment phase

for(i in 1:length(cluster)){
  colors[i] <- ifelse(cluster[i] == 1,"red",
                      ifelse(cluster[i] == 2,"lightblue1",
                             ifelse(cluster[i] == 3,"lightgreen",
                                    ifelse(cluster[i] == 4,"orange",
                                           ifelse(cluster[i] == 5,"purple1")))))
}

colors_transparent <- adjustcolor(colors, alpha.f = pointsalpha)


open3d(windowRect = windowRect)

par3d(zoom = zoom, userMatrix = userMatrix)


plot3d(umap[, 2], umap[, 3], umap[, 4],
       aspect = aspectr,
       col = colors_transparent, size = 1, type='s',
       main = " ",
       xlab = "UMAP_1",
       ylab = "UMAP_2",
       zlab = "UMAP_3"
)

legend3d(x=.30, y=.75, legend = c("Neutrophils","Cluster 1","Cluster 2","	Endothelial cells","Basophils"),
         col = palette, pch = 16, cex = 1.2, inset = c(0.02))



 
 
