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


#############################################
#Create a metadata
meta <- data.frame(
  Description = c("Control", "Group2", "Group2", "Group1", "Group1", "Group2", "Group2", "Group1", "Group2", "Group1", "Group2", "Group1", "Group2")
)

# Row names (satır isimleri) ayarlama
rownames(meta) <- c("Control", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12")

#########################
#reopen the dds file from saved
dds <- readRDS(file = "DESeq2/dds_object.rds")
#save the vst data
vsdata <- vst(dds,blind=FALSE)
vsdata_df <- as.data.frame(assay(vsdata))

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
data <- data.frame(row.names = rownames(vsdata_df))
for (sample in samples) {
  data[[sample]] <- rowMeans(vsdata_df[,sample_columns[[sample]]], na.rm = TRUE)
}

#QUALITY CONTROL


#Heatmap of sample correlations
p <- plot_heatmap(data, type = "samplecor", show_rownames = FALSE)
p

exclude_samples <- c("D06", "D10","D11","D12")
data_f <- data[,!(colnames(data) %in% exclude_samples)]
meta_f <- meta %>%
  filter(!(row.names(.) %in% exclude_samples))

#after removing the samplesHeatmap of sample correlations
p <- plot_heatmap(data_f, type = "samplecor", show_rownames = FALSE)
p

#PCA plot
#plot_PCA(data,meta,PCs = c(1,2),size = 2)

data_f <- filter_by_variance(data_f, n = 2000)

#Set the soft threshold
sft <- SFT_fit(data_f, net_type = "signed hybrid", cor_method = "pearson")
#check the plot
sft$plot
power <- sft$power

'Network'
net <- exp2gcn(
  data_f, net_type = "signed hybrid",module_merging_threshold = 0.80,  SFTpower = power, 
  cor_method = "pearson"
)
#check the module stability
module_stability(data_f, net, nRuns = 5)

# save() fonksiyonunu kullanarak nesneyi kaydetmek
#save(net, file = "net.RData")
#net <- load("net.RData")

#check the list of names in the net
names(net)

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

#network construction

#For brown
edges <- get_edge_list(net, module="brown")
edges_filtered <- get_edge_list(net, module = "brown", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for brown

module_genes <- as.data.frame(net$genes_and_modules)
brown_genes <-  subset(module_genes, Modules == "brown")


#ENSEMBL mapping
df.top <- brown_genes 
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
GO_brown <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_brown, showCategory = 20)

########################

#For salmon

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "salmon",
  bg_line="mean"
)


edges <- get_edge_list(net, module="salmon")
edges_filtered <- get_edge_list(net, module = "salmon", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for salmon

module_genes <- as.data.frame(net$genes_and_modules)
salmon_genes <-  subset(module_genes, Modules == "salmon")


#ENSEMBL mapping
df.top <- salmon_genes 
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
GO_salmon <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_salmon, showCategory = 10)

##########################################################

#For darkred

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "darkred",
  bg_line="mean"
)


edges <- get_edge_list(net, module="darkred")
edges_filtered <- get_edge_list(net, module = "darkred", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for darkred

module_genes <- as.data.frame(net$genes_and_modules)
darkred_genes <-  subset(module_genes, Modules == "darkred")


#ENSEMBL mapping
df.top <- darkred_genes 
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
#GO darkred
GO_darkred <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_darkred, showCategory = 10)

##########################################################

#For red

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "red",
  bg_line="mean"
)


edges <- get_edge_list(net, module="red")
edges_filtered <- get_edge_list(net, module = "red", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for darkred

module_genes <- as.data.frame(net$genes_and_modules)
red_genes <-  subset(module_genes, Modules == "red")


#ENSEMBL mapping
df.top <- red_genes 
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
#GO darkred
GO_red <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_red, showCategory = 10)

##########################################################

#For pink

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "pink",
  bg_line="mean"
)


edges <- get_edge_list(net, module="pink")
edges_filtered <- get_edge_list(net, module = "pink", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for darkred

module_genes <- as.data.frame(net$genes_and_modules)
pink_genes <-  subset(module_genes, Modules == "pink")


#ENSEMBL mapping
df.top <- pink_genes 
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
#GO darkred
GO_pink <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_pink, showCategory = 10)

#For pink

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "pink",
  bg_line="mean"
)


edges <- get_edge_list(net, module="pink")
edges_filtered <- get_edge_list(net, module = "pink", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for darkred

module_genes <- as.data.frame(net$genes_and_modules)
pink_genes <-  subset(module_genes, Modules == "pink")


#ENSEMBL mapping
df.top <- pink_genes 
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
#GO darkred
GO_pink <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_pink, showCategory = 10)
############################################

#For black

#plot the expression profile in the module
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "black",
  bg_line="mean"
)


edges <- get_edge_list(net, module="black")
edges_filtered <- get_edge_list(net, module = "black", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for darkred

module_genes <- as.data.frame(net$genes_and_modules)
black_genes <-  subset(module_genes, Modules == "black")


#ENSEMBL mapping
df.top <- black_genes 
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
#GO darkred
GO_black <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_black, showCategory = 10)

###################################################


#For Grey
plot_expression_profile(
  exp = data_f, 
  metadata=meta_f,
  net = net, 
  plot_module = TRUE, 
  modulename = "grey",
  bg_line="mean"
)



edges <- get_edge_list(net, module="grey")
edges_filtered <- get_edge_list(net, module = "grey", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for lightsteelblue1

module_genes <- as.data.frame(net$genes_and_modules)
grey_genes  <- subset(module_genes, Modules == "grey")


#ENSEMBL mapping
df.top <- grey_genes 
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
GO_grey <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_grey, showCategory = 10)

########################








