

library(clusterProfiler)
library(GOfuncR)
library(msigdbr)
library(ggplot2)
library(org.Mm.eg.db)
library(Mus.musculus)
library(ggpubr)



# We first take a look at the loadings for the first ten PCs

loadings <- data.frame(pca$rotation[,1:10])


# The next steps are just for improving the visualization of the barplot
# for the loadings. 

loadings$Name <- rownames(loadings)  # create new column for car names
loadings <- loadings[,c(11,7)]
loadings$PC7 <- round(loadings$PC7, 3)  

loadings <- loadings[order(loadings$PC7,decreasing = T), ]  # sort
loadings$Name <- factor(loadings$Name, levels = loadings$Name)  # convert to factor to retain sorted order in plot.

# Create a sort of scale to better individuate which genes are the most important for PC7
loadings$type <- ifelse(loadings$PC7 >= 0.1, "top",
                        ifelse(loadings$PC7 < 0.1 & loadings$PC7 >= 0,"above",
                               ifelse(loadings$PC7 <= -0.1, "down","below")))  


# Horizontal Barplot for loadings

theme_set(theme_bw())  
ggplot(loadings, aes(x=Name, y=PC7, label=PC7)) + 
  geom_bar(stat='identity', aes(fill=type), width=.5, show.legend = FALSE) +
  scale_fill_manual(values = c("top"="darkgreen", "above"="green3", "below"="salmon", "down"="red3")) + 
  labs(title= "PC7 Loadings/Genes barplot") + 
  labs(x= " ") +
  coord_flip() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Sfondo interno bianco
        plot.background = element_rect(fill = "white", color = NA), # Sfondo esterno con colore desiderato
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())



# For Gene Ontology analysis the criteria to select the list of genes 
# will be the following: 
# We will consider as input for GO analysis the genes that have loadings above 0.1
# and below -0.1


genes <- as.character(loadings[loadings$PC7 <= -0.1 | loadings$PC7 >= 0.1,1])

Input_Genes <- genes
Input_Genes <- data.frame(Input_Genes)


Input_Genes$Candidati <- 1

#Let's start with GO analysis...

Go_Enrich_Out <- go_enrich(Input_Genes,organismDb = 'Mus.musculus')
Results <- Go_Enrich_Out$results #we store here the results of GO analysis

Over_Representation<-Results[Results$raw_p_overrep<=0.05,]
Under_Representation<- Results[Results$raw_p_underrep<=0.05,]

Genes<- Go_Enrich_Out$genes
Candidate_Gene <-Genes[Genes$Candidati==1,]

Gene_all_GO <- get_anno_categories(Candidate_Gene[,1],database = 'Mus.musculus')

Out_Over_Representation<-merge(Over_Representation,Gene_all_GO, by.x = "node_id",by.y = "go_id", all.x = TRUE, all.y = FALSE)

Out_Under_Representation<-merge(Under_Representation,Gene_all_GO, by.x = "node_id",by.y = "go_id", all.x = TRUE, all.y = FALSE)


Results_Over_Representation<- Out_Over_Representation %>% group_by(node_id) %>% mutate(gene= paste(gene, collapse=",")) %>% unique %>% na.omit
Results_Under_Representation<- Out_Under_Representation %>% group_by(node_id) %>% mutate(gene= paste(gene, collapse=",")) %>% unique %>% na.omit
Results_Over_Representation[ ,c(9,10)] <- NULL 
Results_Under_Representation[ ,c(9,10)] <- NULL

# GO presents a intrisic classification into three different categories:
# 1) Biological Classification
# 2) Molecular Function
# 3) Cellular Component

# Each up-regulated pathway belong to a specific category (this is stored in the column "ontology")

Results_Over_Representation$ontology <- as.factor(Results_Over_Representation$ontology)
table(Results_Over_Representation$ontology)
barplot(table(Results_Over_Representation$ontology)) #Just for visualization...

# Pathways order with respect ontology and significance value
Results_Over_Representation <- Results_Over_Representation %>%
  arrange(ontology, raw_p_overrep)


# For each category, we select the top ten over-represented pathways and transform the 
# "significance value" into a more common measure that involves p-value (-log10 p-value)

first_10_values <- Results_Over_Representation %>%
  group_by(ontology) %>%
  top_n(-10, raw_p_overrep)

first_10_values$log10_p <- -log10(first_10_values$raw_p_overrep)

first_10_values <- data.frame(first_10_values)


df <- data.frame(
  Term = first_10_values$node_name,
  Value = first_10_values$log10_p,
  Category = first_10_values$ontology)

# Unique level for each category --- just check...
df$Category <- factor(df$Category, levels = unique(df$Category))

df[4,1] <- c("positive regulation of chemokine CXCL2 production")


df$Term <- factor(df$Term, levels = df$Term)
df$Category <- factor(df$Category, levels = c("biological_process", "cellular_component", "molecular_function"))

#Visualization...

ggplot(df, aes(x = Value, y = Term)) + 
  geom_bar(stat='identity', aes(fill=Category), width=.5, show.legend = TRUE) + 
  labs(title= " ") + 
  labs(y = " ", x = expression('-log' [10] * '(p-value' * ')')) + 
  theme_minimal() +  # Tema senza griglia
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # Sfondo interno bianco
    plot.background = element_rect(fill = "white", color = NA), # Sfondo esterno con colore desiderato
    panel.border = element_rect(color = "black", fill = NA),  # Aggiungi un box attorno
    panel.grid.major = element_blank(),  # Rimuovi la griglia principale
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),  # Aumenta la grandezza delle etichette asse x
    axis.text.y = element_text(size = 16, color = "black"), # Aumenta la grandezza del testo delle etichette asse y
    legend.position = c(0.7, 0.50),  # Posiziona la legenda in alto a destra
    legend.title = element_text(size = 12),  # Aumenta la grandezza del titolo della legenda
    legend.text = element_text(size = 10)) +  # Aumenta la grandezza del testo della legenda
  scale_x_continuous(limits = c(0, 9)) +  # Estendi l'asse x fino a 9
  scale_fill_manual(values = c("biological_process" = "darkgreen", "cellular_component" = "blue", "molecular_function" = "red"),  # Colori personalizzati per le categorie
                    labels = c("Biological Process", "Cellular Component", "Molecular Function"))  # Etichette personalizzate per la legenda





