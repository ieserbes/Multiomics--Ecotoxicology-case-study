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

#Set the soft threshold
sft <- SFT_fit(data_f, net_type = "signed hybrid", cor_method = "pearson")
#check the plot
sft$plot
power <- sft$power

'Network'
net <- exp2gcn(
  data_f, net_type = "signed hybrid",module_merging_threshold = 0.60,  SFTpower = power, 
  cor_method = "pearson"
)
#check the module stability
module_stability(data_f, net, nRuns = 5)

# save() fonksiyonunu kullanarak nesneyi kaydetmek
save(net, file = "net.RData")
net <- load("net.RData")

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
  modulename = "antiquewhite",
  bg_line="mean"
)

#get hub genes
hubs <- get_hubs_gcn(data_f, net)
head(hubs)

#network construction

#For lightsteelblue1
edges <- get_edge_list(net, module="lightsteelblue1")
edges_filtered <- get_edge_list(net, module = "lightsteelblue1", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for lightsteelblue1

module_genes <- as.data.frame(net$genes_and_modules)
lightsteelblue1_genes <-  subset(module_genes, Modules == "lightsteelblue1")


#ENSEMBL mapping
df.top <- lightsteelblue1_genes 
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
GO_lightsteelblue1 <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

########################

#For Grey
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
barplot(GO_grey, showCategory = 20)

########################

#For deeppink2
edges <- get_edge_list(net, module="deeppink2")
edges_filtered <- get_edge_list(net, module = "deeppink2", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for deeppink2

module_genes <- as.data.frame(net$genes_and_modules)
pink_genes  <- subset(module_genes, Modules == "deeppink2")


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
#GO lightsteelblue1module
GO_pink <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_pink, showCategory = 10)


########################

#For antiquewhite
edges <- get_edge_list(net, module="antiquewhite")
edges_filtered <- get_edge_list(net, module = "antiquewhite", filter = TRUE)

plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

#extract the genes for antiquewhite

module_genes <- as.data.frame(net$genes_and_modules)
antiquewhite_genes  <- subset(module_genes, Modules == "antiquewhite")


#ENSEMBL mapping
df.top <- antiquewhite_genes 
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
GO_antiquewhite <- enrichGO(gene = GO_top, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# plot the results
barplot(GO_antiquewhite, showCategory = 10)




















###################################################
rm(mapping_list, mappings, merged_data,net_edges)
net_edges <- edges_filtered

# load mapping files

#Mapping the ENSEMBL

mappings <- read.delim("annotations/rna_dma_to_hsa_mappings.tsv", header = FALSE, sep = "\t")

# Eşleme Oluşturma
mappings$V2 <- strsplit(as.character(mappings$V2), ";")  # ENSEMBL ID'leri ayırma

# Eşleme listesini genişletme
mapping_list <- setNames(lapply(mappings$V2, function(x) ifelse(is.na(x), NA, x)), mappings$V1)

# Eşleme Uygulama - fromNode ve toNode için
net_edges$Var1 <- sapply(net_edges$Var1, function(x) paste(unique(unlist(mapping_list[x])), collapse = ";"))
net_edges$Var2 <- sapply(net_edges$Var2, function(x) paste(unique(unlist(mapping_list[x])), collapse = ";"))

# Function to select the first ENSEMBL ID where multiple IDs are present
select_first_ensembl <- function(ids) {
  sapply(strsplit(as.character(ids), ";"), function(x) ifelse(length(x) > 0, x[1], NA))
}

# Apply the function to both 'fromNode' and 'toNode' columns

net_edges$Var1 <- select_first_ensembl(net_edges$Var1)
net_edges$Var2 <- select_first_ensembl(net_edges$Var2)

#Öncelikle, fromNode veya toNode'da NA olmayan satırları filtrele
net_edges <- net_edges %>%
  filter(!is.na(Var1) & !is.na(Var2))

# Gerekli kütüphaneyi yükle
library(org.Hs.eg.db)

# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
net_edges$Var1 <- mapIds(org.Hs.eg.db,
                             keys=net_edges$Var1,
                             keytype="ENSEMBL",
                             column="SYMBOL",
                             multiVals="first")

# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
net_edges$Var2 <- mapIds(org.Hs.eg.db,
                         keys=net_edges$Var2,
                         keytype="ENSEMBL",
                         column="SYMBOL",
                         multiVals="first")

plot_gcn(
  edgelist_gcn = net_edges, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)




#############################################################


