

library(ggplot2)
library(cluster)
library(grDevices)

#Set seed to reproduce results
set.seed(123)

# Create three almost-spherical clusters
cluster1 <- cbind(runif(10, min = 0, max = 0.5), runif(10, min = -0.2, max = 0.7))
cluster2 <- cbind(runif(10, min = 1, max = 2), runif(10, min = 1, max = 1.5))
cluster3 <- cbind(runif(10, min = 2.5, max = 3.5), runif(10, min = 0, max = 1))

# Create a collinear cluster
cluster4 <- cbind(c(1.3, 1.5, 1.7, 1.9, 2.1), c(0.1, 0.2, 0.4, 0.6, 0.7))

# Combina tutto
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





# Application of K-Means clustering for K=4 clusters 

km <- kmeans(df[, 1:2], centers = 4)
df$km_cluster <- factor(km$cluster) #K-Means assignment

# Set different colors for each cluster
colors <- c("red", "skyblue1", "green3", "violet")
names(colors) <- levels(df$km_cluster)



par(mgp = c(2, 0.4, 0)
    

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




# Application of K-Area clustering for K=4 clusters 

# First we create a function for the K-Area Clustering (2D-Space)

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



# Then we apply the previous function setting K=4 

k <- 4
clustering <- k_center_2d(df[,1:2], k)

# We store the results in different R object
centers <- clustering$centers # centroids
assignments <- clustering$assignments # label assignment
df$Cluster_K_vol <- factor(assignments) #K-Area assignment

# Set different colors for each cluster
colors <- c("red", "skyblue1", "green3", "violet")
names(colors) <- levels(df$Cluster_K_vol)



par(mgp = c(2, 0.4, 0))

plot(df$V1, df$V2, pch = 21, bg = colors[df$Cluster_K_vol], col = colors[df$Cluster_K_vol], cex = 1.2,
     xlab = "V1", ylab = "V2", main = "K-Area Clustering with Convex Hulls",xaxt = "n", yaxt = "n")

# Custom x-axis
axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
     tck = -0.02, mgp = c(0, 0.4, 0))

# Custom y-axis
axis(2, at = c(0, 0.5, 1, 1.5),
     labels = c(0, 0.5, 1, 1.5),
     tck = -0.02, mgp = c(0, 0.4, 0))


# Add Convex-Hulls to delimitate clusters
for (cl in levels(df$Cluster_K_vol)) {
  pts <- df[df$Cluster_K_vol == cl, c("V1", "V2")]
  if (nrow(pts) >= 3) {
    hull <- chull(pts)
    polygon(pts[hull, ], border = colors[cl], lwd = 2)
  }
}

# Centroids
points(centers[,1], centers[,2], pch = 3, cex = 2, lwd = 2)

legend("topleft", legend = paste("Cluster", levels(df$Cluster_K_vol)),
       pt.bg = colors, pch = 21, pt.cex = 1.2, col = colors)

