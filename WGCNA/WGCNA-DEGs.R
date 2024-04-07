rm(list = ls())

#Set the working directory
setwd("D:/Academic2024/MSc_Bioinformatics/Modules/Module-10-Computational Biology for Complex Systems-DL BIO - 37415/Data/Data2")


#########################################
#Required packages
suppressMessages({
  library(dplyr)
  library(BioNERO)
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(EnhancedVolcano)
  library(textshaping)
  library(WGCNA)
  library(CorLevelPlot)
  library(gridExtra)
  library(pheatmap)
  #library(limma)
  library(tidyverse)
  #Go enrichment analysis
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

d1 <- read.csv("DESeq2/deseq_res_G1_sig_genes_1Xwithsmybols_notfilteredwithpadj.csv",row.names=1,header=TRUE)



dds <- readRDS(file = "DESeq2/dds_object.rds")
#save the vst data
vsdata <- vst(dds,blind=FALSE)
vsdata_df <- as.data.frame(assay(vsdata))



rownames <- (rownames(d1))

# vsdata_df'de combined_rownames ile eşleşen satırlar bulunur
matched_rows <- vsdata_df[rownames(vsdata_df) %in% rownames, ]

# Eşleşen satırları yeni bir dataframe'e aktar
data <- matched_rows


#gruplari birlestir

#group the sites (group the replicates together)
# Define the groups and corresponding column names
samples <- c("Control", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12")
sample_columns <- list(
  Control = c("CK07", "CK08", "CK09", "CK10", "CK11", "CK12")
)

# Generate column names for D01 to D12 based on observed pattern
for (sample in samples[-1]) { # Skip "Control" since it's already defined
  sample_columns[[sample]] <- paste0(sample, "A", 1:6)
}

# Calculate the mean for each group and combine into a new dataframe
data2 <- data.frame(row.names = rownames(data))
for (sample in samples) {
  data2[[sample]] <- rowMeans(data[,sample_columns[[sample]]], na.rm = TRUE)
}

data <- data2

##############################
#Create a metadata
meta <- data.frame(
  Description = c("Control", "Group2", "Group2", "Group1", "Group1", "Group2", "Group2", "Group1", "Group2", "Group1", "Group2", "Group1", "Group2")
)

# Row names (satır isimleri) ayarlama
rownames(meta) <- c("Control", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12")

#QUALITY CONTROL


#Heatmap of sample correlations
p <- plot_heatmap(data, type = "samplecor", show_rownames = FALSE)
p

#PCA plot
plot_PCA(data,meta,PCs = c(1,2),size = 2)


exclude_samples <- c("D06", "D10","D11","D12")
data_f <- data[,!(colnames(data) %in% exclude_samples)]
meta_f <- meta %>%
  filter(!(row.names(.) %in% exclude_samples))

#after removing the samples Heatmap of sample correlations
p <- plot_heatmap(data_f, type = "samplecor", show_rownames = FALSE)
p

##after removing the samples-PCA plot
plot_PCA(data_f,meta_f,PCs = c(1,2),size = 2)

#Set the soft threshold
sft <- SFT_fit(data_f, net_type = "signed", cor_method = "pearson", )
#check the plot
sft$plot
power <- sft$power

'Network'
net <- exp2gcn(
  data_f, net_type = "signed",module_merging_threshold = 0.8,  SFTpower = power, 
  cor_method = "pearson"
)

#check the module stability
module_stability(data_f, net, nRuns = 5)


# Dendro and colors
plot_dendro_and_colors(net)

# Eigengene networks plot (similarity between modules)
plot_eigengene_network(net)

#how many genes in modules
plot_ngenes_per_module(net)

#Module trait correlation plot (which module which group)
MEtrait <- module_trait_cor(exp = data_f, metadata=meta_f, MEs = net$MEs)
head(MEtrait)
plot_module_trait_cor(MEtrait)

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "brown",
  bg_line="mean"
)

#get hub genes
hubs <- get_hubs_gcn(data_f, net)
head(hubs)

############################################
#network construction -group1 specific

#For darkslateblue
edges <- get_edge_list(net, module="darkslateblue")
edges_filtered <- get_edge_list(net, module = "darkslateblue", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
darkslateblue_genes <-  subset(module_genes, Modules == "darkslateblue")


#ENSEMBL mapping
df.top <- darkslateblue_genes 
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_darkslateblue <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_darkslateblue, showCategory = 20)
####################################
#For tan4
edges <- get_edge_list(net, module="tan4")
edges_filtered <- get_edge_list(net, module = "tan4", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
tan4_genes <-  subset(module_genes, Modules == "tan4")


#ENSEMBL mapping
df.top <- tan4_genes 
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_tan4 <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_tan4, showCategory = 20)

####################################
#For green
edges <- get_edge_list(net, module="green")
edges_filtered <- get_edge_list(net, module = "green", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
green_genes <-  subset(module_genes, Modules == "green")


#ENSEMBL mapping
df.top <- green_genes 
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_green <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_green, showCategory = 20)



####################################
#For coral4
edges <- get_edge_list(net, module="coral4")
edges_filtered <- get_edge_list(net, module = "coral4", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
coral4_genes <-  subset(module_genes, Modules == "coral4")


#ENSEMBL mapping
df.top <- coral4_genes 
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_coral4 <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_coral4, showCategory = 10)
#####################################################################################
#For blue2
edges <- get_edge_list(net, module="blue2")
edges_filtered <- get_edge_list(net, module = "blue2", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
blue2_genes <-  subset(module_genes, Modules == "blue2")


#ENSEMBL mapping
df.top <- blue2_genes 
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_blue2 <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_blue2, showCategory = 10)


#####################################################################################
#For midnightblue
edges <- get_edge_list(net, module="midnightblue")
edges_filtered <- get_edge_list(net, module = "midnightblue", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
midnightblue_genes <-  subset(module_genes, Modules == "midnightblue")


#ENSEMBL mapping
df.top <- midnightblue_genes
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_midnightblue <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_midnightblue, showCategory = 10)

#####################################################################################
#For bisque4
edges <- get_edge_list(net, module="bisque4")
edges_filtered <- get_edge_list(net, module = "bisque4", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
bisque4_genes <-  subset(module_genes, Modules == "bisque4")


#ENSEMBL mapping
df.top <- bisque4_genes
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_bisque4 <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_bisque4, showCategory = 10)


#####################################################################################
#For purple
edges <- get_edge_list(net, module="purple")
edges_filtered <- get_edge_list(net, module = "purple", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
purple_genes <-  subset(module_genes, Modules == "purple")


#ENSEMBL mapping
df.top <- purple_genes
ids <- df.top$Genes

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, sep = '\t', stringsAsFactors = F)
colnames(id.maps) <- c("Daphnia_gene", "ENSEMBL_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, Daphnia_gene %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$ENSEMBL_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "Genes", by.y = "Daphnia_gene", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("ENSEMBL_ids", "Daphnia_gene"))]

GO_top <- df.top$ENSEMBL
#GO lightsteelblue1module
GO_purple <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_purple, showCategory = 10)
################################################################################









