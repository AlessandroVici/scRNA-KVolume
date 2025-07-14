





# Create a dataframe object that contains, for each cells, the scores related to
# the minor PCs (i.e. PCs not used in the clustering) and specifies the cluster to which each cell belongs

# df is the dataframe that stored the result of K-Volume clustering

df$Cells <- rownames(df)

# Remove the PCs used for clustering and all PCs over PC10
scores_pca <- scores_pca[,-c(1,2,3,11:50)] 

# The resulting dataframe will be:

clustering_minor_PC <- merge(df,scores_pca,by="Cells")

# Some improvements...

rownames(clustering_minor_PC) <- clustering_minor_PC$Cells
clustering_minor_PC <- clustering_minor_PC[,c(5,1:4,6:12)]
clustering_minor_PC <- clustering_minor_PC[,-c(2:5)]

# We check for number of cells inside the clusters...

table(clustering_minor_PC$Cluster)

# As cluster 5 has only one cell, we decided to remove for further analysis

clustering_minor_PC <- clustering_minor_PC[!c(clustering_minor_PC$Cluster == 5),]

# We have to reorder the levels 

clustering_minor_PC$Cluster <- as.character(clustering_minor_PC$Cluster)
clustering_minor_PC$Cluster <- as.factor(clustering_minor_PC$Cluster)

# Change the names of the levels
levels(clustering_minor_PC$Cluster) <- c("C1", "C2", "C3", "C4")

# We create also a vector that contains the number of cells inside each cluster

Cluster_num <- as.numeric(table(clustering_minor_PC$Cluster))



#------Reasoning behind the minor component signal/significance validation algorithm------#


# The statistical significance of the minor components is calculated (distance from the standard normal in terms of
# mean zero and standard deviation 1, therefore in practice, for very small clusters, Fisher's exact test on the null hypothesis
# of 5% of observations outside the range of two standard errors to the left and right,
# on the clusters obtained with the major components. Since the minor components were NOT used for
# clustering, this is effectively an 'external validation' test of the significance of the minor components.

# We use the Fisher's exact test if the cluster has fewer than 30 observations, otherwise the chi-square test.

# When you have a small cluster and want to check whether its observations show anomalous behavior
# along a given principal component, such as the sixth, you can apply the Fisher's exact test.
# To do this, we start by defining what is meant by "anomaly": in this context, an anomalous value is considered to be
# that is found at outside the range of two standard deviations above or below the mean, as expected
# from a standard normal distribution.

# Once this criterion is established, the observations in the cluster that fall outside this range are counted
# and those that do fall within it. These numbers are then compared to what would be expected if the data
# actually followed a standard normal distribution, in which only a small percentage of
# observations should fall outside two standard deviations.

# At this point, a two-way table is constructed, in which one row represents the observed data
# in the cluster and the other row represents the expected values according to the theoretical distribution. The first column
# of the table contains the number of observations outside the range, while the second column contains those within.

# The Fisher exact test is then applied to this table to calculate the probability of observing such a
# deviation, or even more extreme, simply by chance. If this probability, i.e., the p-value,
# is sufficiently low, We can conclude that the cluster exhibits anomalous behavior along that
# principal component, suggesting that it may contain relevant information not
# captured by the larger principal components.



# The following code creates an empty matrix that will store the p-value obtained from the statistical test
# The matrix has rows equal to the number of clusters
# and columns equal to the number of minor PCs (from PC4 to PC10)

matrix_pval <- data.frame(matrix(NA,nrow=4, ncol=7))
colnames(matrix_pval) <- c("PC4","PC5","PC6","PC7","PC8","PC9","PC10")
rownames(matrix_pval) <- c("C1","C2","C3","C4")


# For each cluster...

for(i in 1:length(Cluster_num)){
  
  current_cluster_label <- levels(clustering_minor_PC$Cluster)[i]
  cat("Current cluster: ", levels(clustering_minor_PC$Cluster)[i], "that contains ", Cluster_num[i], "cells", "\n")
  
  # We test all the minor PCs inside the dataframe
  # As the first column is for the clustering labels...
  
  for(j in 2:ncol(clustering_minor_PC)){

    cat("Cluster ", levels(clustering_minor_PC$Cluster)[i], "related to ", colnames(clustering_minor_PC[j]), "\n")  
    
    matrix_test <- data.frame(matrix(NA,nrow=2,ncol=2))
    rownames(matrix_test) <- c("Current","Expected")
    colnames(matrix_test) <- c("Out ±2 SD","In ±2 SD")
    
    # We fill the second row of the 2x2 matrix
    matrix_test[2,1] <- round((Cluster_num[i]/100)*5,0)
    matrix_test[2,2] <- round((Cluster_num[i]/100)*95,0)
    
    # We fill the first row of the 2x2 matrix
    current <- clustering_minor_PC[clustering_minor_PC$Cluster == current_cluster_label,]
    
    matrix_test[1,1] <- sum(current[,j] >= 2 | current[,j] <= -2)
    matrix_test[1,2] <- Cluster_num[i] - matrix_test[1,1]
    
    if(Cluster_num[i] <= 30){
      
      # Exact Fisher test
      risultato <- fisher.test(matrix_test)
      
      # P-value
      p <- risultato$p.value
      
      matrix_pval[i,j-1] <- p
      
      cat("Exact Fisher test between ", levels(clustering_minor_PC$Cluster)[i], "and ", colnames(clustering_minor_PC[j]), "has p-value: ", p , "\n") 
    }
    else{
      # Chi-Square Test
      risultato <- chisq.test(matrix_test)
      
      # P-value
      p <- risultato$p.value
      
      matrix_pval[i,j-1] <- p
      
      cat("Chi-Square Test between ", levels(clustering_minor_PC$Cluster)[i], "and ", colnames(clustering_minor_PC[j]), "has p-value: ", p , "\n") 
    }
  }
}



