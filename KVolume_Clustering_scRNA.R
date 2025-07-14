
library(rgl)
library(geometry)




# We use the pca object obtained previously 
# Extraction of the PCA scores

scores_pca <- data.frame(pca$x)

# For our purpose, we select the first three PC in terms of scores

df <- scores_pca 
df <- df[, c(1, 2, 3)]  #prime 3 PCs
data <- as.matrix(df)


# K-Volume Clustering functions

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

#This function allows to create 3D-Convex Hull for the clusters

plotbag_allpoints <- function(x, cluster_id, col = "black") {
  x2 <- x[x[, ncol(x)] == cluster_id, ]
  
  if (nrow(x2) >= 4) {
    ids <- t(convhulln(as.matrix(x2[, 1:3]), options = "Tv"))
    triangles3d(x2[ids, 1], x2[ids, 2], x2[ids, 3],
                col = col, alpha = 0.5, shininess = 80)
    spheres3d(x2[ids, 1], x2[ids, 2], x2[ids, 3], r = 0.1, color = col) 
  }
}




# For better visualization we used some graphical parameters

zoom <- 0.8
aspectr <- c(1, 1, 1)
windowRect <- c(100, 100, 800, 600)

cex <- 1
pointsalpha <- 1
userMatrix <- matrix(c(0.80, -0.60, 0.022, 0,
                       0.23,  0.34, 0.91,  0,
                       -0.55, -0.72, 0.41,  0,
                       0,     0,    0,     1), ncol = 4, byrow = TRUE)


# Lacunarity Index calculation

# Let's start with the total Convex-Hull calculation for all the points (cells)
total_volume <- convhulln(data, options = "FA")$vol

# We repeat K-Volume adopting different possible clustering solutions 

K_low <- 3
K_max <- 10

# We initialize an empty vector that will contain the Convex-Hull value
# for various K solutions

cluster_volume_sums <- numeric()


for (k in K_low:K_max) {
  clustering <- k_center(data, k)
  centers <- clustering$centers
  assignments <- clustering$assignments
  df$Cluster <- factor(assignments)
  
  cat("\n--- k =", k, "---\n")
  cluster_volumes <- numeric(k) # numeric vector that reports volume of the single clusters inside the current K solution
  
  for (i in 1:k) {
    group_data <- df[df$Cluster == i, 1:3]
    
    # If the cluster has 3 or less units (cells), the Convex-Hull volume won't be calculated
    if (nrow(group_data) >= 4) {
      vol <- convhulln(as.matrix(group_data), options = "FA")$vol
    } else {
      vol <- 0
    }
    
    cluster_volumes[i] <- vol
    cat("Cluster volume", i, ":", vol, "\n")
  }
  
  cluster_volume_sums[k] <- sum(cluster_volumes)
  cat("Sum of cluster volumes for k =", k, ":", cluster_volume_sums[k], "\n")
}


cluster_volume_sums <- na.omit(cluster_volume_sums)

#Lacunarity Index for each K tested
lacunarity_k <- (total_volume - cluster_volume_sums) / total_volume
cat("Lacunarity Index for k =", k, ":", lacunarity_k, "\n")


# The following plot is used to show graphically the optimal K-solution 
# based on Lacunarity Index 

par(mgp = c(2, 0.4, 0))

plot(c(3:10), lacunarity_k, type = "b",
     xlab = "Number of clusters (K)",
     ylab = "Lacunarity index",
     main = " ",
     xaxt = "n", yaxt = "n",
     panel.first = rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
                        col = "white", border = NA),
     pch = 16, col = "grey30")

# Custom x-axis
axis(1, at = c(3, 5, 7, 9), labels = c(3, 5, 7, 9),
     tck = -0.02, mgp = c(0, 0.4, 0))

# Custom y-axis
axis(2, at = c(0.90, 0.92, 0.94, 0.96, 0.98),
     labels = c("0.90", "0.92", "0.94", "0.96", "0.98"),
     tck = -0.02, mgp = c(0, 0.4, 0))

points(5, lacunarity_k[3], pch = 1, col = "red", cex = 3, lwd = 2)





# Graphical representation of K-Volume clustering for optimal value of K

k <- 5
clustering <- k_center(data, k)
centers <- clustering$centers
assignments <- clustering$assignments
df$Cluster <- factor(assignments)

# Set colors for clusters
cols <- rainbow(k)

# Open 3D window using graphical parameters
open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect, antialias = 8)


cluster_ids <- sort(unique(assignments)) 
for (i in seq_along(cluster_ids)) {
  plotbag_allpoints(df, cluster_id = cluster_ids[i], col = cols[i])
}

axes3d(color = "black", drawfront = TRUE, box = TRUE, alpha = 1, labels = FALSE, tick = TRUE)

title3d(main = " ",
        xlab = "PC1" ,
        ylab = "PC2",
        zlab = "PC3",
        color = "black")

#aspect3d(aspectr)

# Centroids
texts3d(centers[, 1], centers[, 2], centers[, 3],
        texts = "+", col = "black", cex = 1.5, adj = c(0.5, 0.5))

legend3d("topright", legend = paste("Cluster", cluster_ids),
         col = cols, pch = 16, cex = 1.2, inset = c(0.02))

