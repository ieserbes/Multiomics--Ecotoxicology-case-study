###title: "DESeq2-Group1-Group2-10X"###
#clean the memory
rm(list = ls())

#Set the working directory
setwd("D:/Academic2024/MSc_Bioinformatics/Modules/Module-10-Computational Biology for Complex Systems-DL BIO - 37415/Data/Data2")


#########################################
#Required packages
suppressMessages({
  library(dplyr)
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(EnhancedVolcano)
  library(textshaping)
  library(WGCNA)
  library(pheatmap)
  #library(limma)
  library(tidyverse)
  #Go enrichment analysis
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

#############################################
#Load the data
##load the meta data

meta <- read.csv('sample_sheet.csv', row.names = 1, header = TRUE, sep = ',')

# Function for loading the data
load.dataset <- function(data.file, data.sep = ',', meta) {
  data <- read.table(data.file, header = TRUE, sep = data.sep, row.names = 1, check.names = FALSE)
  data <- as.matrix(data)
  
  # Meta veri ve ana veri arasındaki eşleşmeyi bul
  mids <- match(row.names(meta), colnames(data))
  data <- t(data[,mids[!is.na(mids)]])
  matched_meta <- meta[!is.na(mids),]
  
  stopifnot(row.names(matched_meta) == row.names(data))
  
  return(list(meta.data = matched_meta, data.matrix = data))
}

# Load the RNAseq data and connect with metadata
dataset <- list()
dataset$RNA_seq <- load.dataset(data.file = 'rna_raw_counts.csv', meta = meta)
#Assign the RNA counts in data 
data <- dataset$RNA_seq$data.matrix

#Filter the dataset as I want to use only 1X data 

patterns_to_exclude <- c("B") #Here I can exclude more things if I want

# Filter rows from the data matrix 
filtered_data <- data[!grepl(paste(patterns_to_exclude, collapse = "|"), rownames(data)), ]

# Transpose the filtered data
count_data <- t(filtered_data)

# round the count data 
count_data <-  round(count_data)

#assign the meta_rna
meta_rna <- meta
#Assign the groups
Group1 <- c("D03", "D04", "D07" , "D09","D11")
Group2 <- c("D01","D02", "D05", "D06", "D08","D10", "D12")

# Reorganise the Description by assigning the groups
row_names <- row.names(meta_rna)

meta_rna <- meta_rna %>%
  mutate(Description = case_when(
    Site %in% Group1 ~ "Group1",
    Site %in% Group2 ~ "Group2",
    TRUE ~ Description # Eğer Site Group1 veya Group2 ile eşleşmiyorsa mevcut Description değerini koru
  ))

# Rewrite the rownames
rownames(meta_rna) <- row_names

#Filter the metadata
filtered_meta_rna <- meta_rna[!grepl(paste(patterns_to_exclude, collapse = "|"), rownames(meta_rna)), ]
metadata <- filtered_meta_rna

##############################################################

#Some Checks before DEseq2

#1.check the number of columns in count_data same in number of rows in metadata
ncol(count_data) == nrow(metadata)

#2.check their names in the same order
all((colnames(count_data))==rownames(metadata))

#3.Box plot for raw counts- for outlier detection 
# Sıfır değerlerle karşılaşmamak için 1 ekleyerek log2 dönüşümü yap
log2_count_data <- log2(count_data + 1)

# log2 dönüştürülmüş count_data matrisini uzun formata dönüştür
log2_count_data_long <- melt(log2_count_data)

# Sütun adlarını düzenle
names(log2_count_data_long) <- c("Gene", "Sample", "Log2Count")

# ggplot ile boxplot çiz
ggplot(log2_count_data_long, aes(x = Sample, y = Log2Count,fill = Sample)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + # X ekseni etiketlerini eğik yap
  labs(title = "Sample-wise Log2 Transformed Gene Raw Count Distribution", x = "Sample", y = "Log2(Count + 1)")
# replikalar arasi cok fazla varyans degisimi var. Medianlari farkli. Belki DEseq2 nin normalizasyonu bunu cozer.

############################################
#4.Filtreleme
#RNAseq data cok fazla zero-count (0) sahip. Bu sikintili sonuc dogurabiliyor. Ozellikle p-value distributionda bu
#sorun rahatlikla gozlemleniyor. Dogru bir istatistiksel sonuc icin bu countlari eliyorum. 

#filtreleme oncesi gen sayisina bak
dim(count_data)
#26840 gen

#Minimum threshold for count

threshold_per_sample <- 2
total_threshold <- threshold_per_sample * ncol(count_data) # Örnek sayısı ile çarp

# Toplam sayım değeri hesaplama ve filtreleme
total_count_filter <- rowSums(count_data) >= total_threshold

# İlk filtre uygulama
filtered_count_data <- count_data[total_count_filter,]

# İkinci koşul: Örneklerin belirli bir yüzdesinden fazlasında sıfır olmayan genleri koru
filter_threshold <- 0.8 # Sıfır sayımlarının eşik değeri
sufficient_expression_filter <- rowSums(filtered_count_data == 0) / ncol(filtered_count_data) <= filter_threshold

# İkinci filtreyi uygula
count_data <- filtered_count_data[sufficient_expression_filter,]

#Filtreleme sonrasi gen sayisina bak
dim(count_data)
#17943


###############################################################
#DESeq2 analizi

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ Description)

#Set the factor level
dds$Description <- relevel(dds$Description,ref="Control")

#Perform DESeq
dds <- DESeq(dds)


#Save the dds file
#saveRDS(dds, file = "dds_object.rds")
#reopen the dds file from saved
#dds <- readRDS(file = "dds_object.rds")

##QC check the results
#check the outliers by using WGCNA package function called goodsamplesgenes

#save the vst data
vsdata <- vst(dds,blind=FALSE)

#Perform QC check the results
#PCA plot
plotPCA(vsdata,intgroup="Description")
#groups were not separated
##Save the plot
# PNG dosyası olarak kaydetmek için bir grafik cihazı aç
png("PCA_1X_G1_G2_DeSeq2_QC.png", width = 400, height = 400)

# plotDispEsts ile dispersion estimates grafiğini çiz
plotPCA(vsdata,intgroup="Description")

# Grafik cihazını kapat ve grafiği kaydet
dev.off()

vsdata_df <- as.data.frame(assay(vsdata))

#check the samples and genes that are outliers with vst counts
gsg <- goodSamplesGenes(t(vsdata_df))
summary(gsg)
gsg$allOK
#It seems like I don't have any outlier and samples.
#check any outliers
table(gsg$goodGenes)
table(gsg$goodSamples)


library(reshape2)
#Box plot icin uzun format
#Gen adlarini sütuna dönüstür
vsdata_df_long <- vsdata_df
vsdata_df_long$row_names <- rownames(vsdata_df_long)

vsdata_df_long <- melt(vsdata_df_long, variable.name = "Sample", value.name = "VstCount", rownames = "Gene")

# Sütun adlarını düzenle
names(vsdata_df_long) <- c("Gene", "Sample", "VstCount")

# ggplot ile boxplot çiz
ggplot(vsdata_df_long, aes(x = Sample, y = VstCount,fill = Sample)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + # X ekseni etiketlerini eğik yap
  labs(x = "Sample", y = "vst counts")

#check the contrast
resultsNames(dds)

#Histogram of p values for all tests. The area shaded in blue indicates the subset of those that pass the filtering, the area in khaki those that do not pass:
res <- results(dds, alpha=0.05)
use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


#dispersion plot
#Plot Dispersion estimate
# PNG dosyası olarak kaydetmek için bir grafik cihazı aç
png("DispEsts_plot_DESeq2.png", width = 800, height = 600)

# plotDispEsts ile dispersion estimates grafiğini çiz
plotDispEsts(dds)

# Grafik cihazını kapat ve grafiği kaydet
dev.off()


#MA-plot 
png("MA-plot for Group1_Group2_Control_DESeq2.png", width = 600, height = 400)
plotMA(res, alpha=0.05, ylim=c(-2,2))
dev.off()

"In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the 
mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored blue if the 
adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles 
pointing either up or down."


# Set up the PNG device with the specified parameters
#png("cooks_distances_boxplot_for_outliers.png", width = 15, height = 6, units = "in", res = 300)

# Create the boxplot
#boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)


# Close the PNG device
#dev.off()
#########################################
# Plot counts for the gene with the smallest adjusted p-value

res <- results(dds)
d <- plotCounts(dds, gene = which.min(res$padj), intgroup = "Description", returnData = TRUE)

# Create the ggplot object
p <- ggplot(d, aes(x = Description, y = count, color = Description)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) + 
  scale_y_log10(breaks = c(25, 100, 400)) +
  labs(title = paste("Gene with lowest padj value:", rownames(res)[which.min(res$padj)]))  # Add title with gene name

# Save the plot as a file (e.g., PNG format)
ggsave("gene_count_for_lowest_padj_value_plot.png", p, width = 6, height = 4, units = "in")

print(p)

###############################################
#Heatmap for each samples

# Assuming dds is your DESeqDataSet object
ntd <- normTransform(dds)

# Select the top 20 rows based on row means
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]

# Set up the PNG device with the specified parameters
png("normalized_counts_heatmap.png", width = 10, height = 8, units = "in", res = 300)

# Create the heatmap plot
heatmap <- pheatmap(assay(ntd)[select, ], cluster_rows = FALSE, show_rownames = FALSE,
                    cluster_cols = FALSE, annotation_col = df)

# Close the PNG device
dev.off()
# Print the heatmap plot to console
print(heatmap)

################################################
#Another heatmap plot with group of sites

library("pheatmap")

# Assuming dds is your DESeqDataSet object
ntd <- normTransform(dds)

# Group samples by site
sample_groups <- colData(dds)$Site

# Calculate average normalized counts for each group
group_means <- aggregate(t(assay(ntd)) ~ sample_groups, FUN = mean)

group <- t(group_means)
colnames(group) <- group[1, ]  # Extract first row and assign to column names
group <- group[-1, ]           # Remove the first row (now headers) from the data

write.csv(group,file="off.csv")
off <- read.csv("off.csv", row.names=1,header=TRUE)
type(off)

# Data frame
Group1 <- c("D03", "D04", "D07" , "D09","D11")
Group2 <- c("D01","D02", "D05", "D06", "D08","D10", "D12")
samples <- c("Control", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12")

# Her örneği ilgili grupla eşleştirin
group_labels <- sapply(samples, function(x) {
  if (x %in% Group1) {
    "Group1"
  } else if (x %in% Group2) {
    "Group2"
  } else {
    "Control"
  }
})

# Veri çerçevesini oluşturun
data_frame <- data.frame(
  Sample = samples,
  Description = group_labels,
  row.names = samples
)
#Remove samples column
data_frame$Sample <- NULL

# Top 20 geni seçin
select <- head(order(rowMeans(off), decreasing = TRUE), 20)

# Extract the sample descriptions from colData
df <- as.data.frame(colData(dds)[, "Description", drop = FALSE])

# Heatmap'i oluşturun (burada satır isimlerini göstermek istemiyorsanız show_rownames = FALSE yapın)

heatmap <- pheatmap(off[select, ], cluster_rows = FALSE, show_rownames = TRUE, 
                    cluster_cols = FALSE, annotation_col = data_frame)

# PNG dosyasına çizimi kaydet
png("normalized_counts_heatmap_grouped.png", width = 10, height = 8, units = "in", res = 300)

# Isı haritasını çizdir
print(heatmap)

# PNG çizimini kapat
dev.off()
##################################
#A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples.
library("RColorBrewer")

# Calculate sample distances
sampleDists <- dist(t(assay(vsdata)))

# Create a matrix from sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata$Site, vsdata$type, sep = "-")
colnames(sampleDistMatrix) <- paste(vsdata$Site, vsdata$type, sep = "-")

# Define colors for heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Set up the PNG device with the specified parameters

# Create the heatmap plot
heatmap2 <- pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
png("heatmap_for site correlations.png", width = 12, height = 10, units = "in", res = 300)



# Print the heatmap plot to console
print(heatmap2)

# Close the PNG device
dev.off()
#################################################
#RESULTS
#Get the group1 results 
res_G1 <- results(dds, contrast=c("Description","Group1","Control"))
res_G1 <- na.omit(res_G1)

#filter the baseMean above 50 (lower is noise)
res_G1 <- res_G1[res_G1$baseMean>50,]
#convert to df
res_G1_df <- as.data.frame(res_G1)


###########################
#Before volcano plot
#ENSEMBL mapping
df.top <- res_G1_df 

ids <- row.names(df.top)

# load mapping files
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
# Convert rownames to a new column for both dataframes
df.top$Daphnia_gene <- rownames(df.top)
# Convert rownames to a column
id.maps.matched$Daphnia_gene <- rownames(id.maps.matched)

# Split the ENSEMBL IDs and take only the first one
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$V2, split = ';'), `[`, 1)

# Now, merge the ID map with your DESeq2 results without expanding
merged_data <- merge(df.top, id.maps.matched, by.x = "Daphnia_gene", by.y = "row.names", all.x = TRUE)

# Set the 'Daphnia_gene' column as rownames
row.names(merged_data) <- merged_data$Daphnia_gene

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene <- NULL

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene.y <- NULL

# Also remove the 'V2' column as it's no longer needed
merged_data$V2 <- NULL


library(org.Hs.eg.db)
# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
merged_data$Symbol <- mapIds(org.Hs.eg.db,
                             keys=merged_data$ENSEMBL,
                             keytype="ENSEMBL",
                             column="SYMBOL",
                             multiVals="first")

#for saving the data frame (bu charactere donusuyor sonrasinda kullanma)
merged_data_df <- merged_data
merged_data_df$Symbol <- sapply(merged_data_df$Symbol, function(x) paste(x, collapse = ", "))
merged_data_df <- as.data.frame(merged_data_df)
write.csv(merged_data_df,file="deseq_res_G1_sig_genes_1Xwithsmybols_notfilteredwithpadj.csv",row.names = TRUE)

##################
#assign the merged_data to 
res_G1_df <- merged_data
#volcano plot for Group1
library(EnhancedVolcano)
library(textshaping)

# Grafiği çıktı olarak ayarla
png("volcano_plot_Group1_1X.png", width = 10, height = 6, units = "in", res = 300)

# Labels
labels <- res_G1_df$Symbol

# EnhancedVolcano fonksiyonunu çağırma
EnhancedVolcano(res_G1_df, x = "log2FoldChange", y = "padj", lab = labels,
                pCutoff = 1e-2, FCcutoff = 0.8, labSize = 3.0, pointSize = 1.0, xlim = c(-5, 5))

# Çıktıyı kapat
dev.off()
###################
#Filter based on padj<0.05

# First filter the p_Value is lower than 0.05
res_G1_sig <- results(dds, alpha=0.05,contrast=c("Description", "Group1", "Control"))

#check the summary
summary(res_G1_sig)
#122 upregulated
#268 downregulated

#drop the NAs
res_G1_sig <- na.omit(res_G1_sig)
#filter the significant genes based on (padj<0.05)
res_G1_sig <- res_G1_sig[res_G1_sig$padj<0.05,]
#convert to a df
res_G1_sig <- as.data.frame(res_G1_sig)

#Mapping again:

#Before volcano plot
#ENSEMBL mapping
df.top <- res_G1_sig 

ids <- row.names(df.top)

# load mapping files
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
# Convert rownames to a new column for both dataframes
df.top$Daphnia_gene <- rownames(df.top)
# Convert rownames to a column
id.maps.matched$Daphnia_gene <- rownames(id.maps.matched)

# Split the ENSEMBL IDs and take only the first one
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$V2, split = ';'), `[`, 1)

# Now, merge the ID map with your DESeq2 results without expanding
merged_data <- merge(df.top, id.maps.matched, by.x = "Daphnia_gene", by.y = "row.names", all.x = TRUE)

# Set the 'Daphnia_gene' column as rownames
row.names(merged_data) <- merged_data$Daphnia_gene

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene <- NULL

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene.y <- NULL

# Also remove the 'V2' column as it's no longer needed
merged_data$V2 <- NULL


library(org.Hs.eg.db)
# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
merged_data$Symbol <- mapIds(org.Hs.eg.db,
                             keys=merged_data$ENSEMBL,
                             keytype="ENSEMBL",
                             column="SYMBOL",
                             multiVals="first")
#Assign to res_G1_sig
res_G1_sig <- merged_data
#for saving the data frame (bu charactere donusuyor sonrasinda kullanma)
merged_data_df <- merged_data
merged_data_df$Symbol <- sapply(merged_data_df$Symbol, function(x) paste(x, collapse = ", "))
merged_data_df <- as.data.frame(merged_data_df)
write.csv(merged_data_df,file="deseq_res_G1_sig_genes_1X_withsymbolwithfilterpadj_005.csv",row.names = TRUE)

######################################################
#Reopen again
res_G1_sig <- read.csv("deseq_res_G1_sig_genes_1X_withsymbolwithfilterpadj_005.csv",row.names=1,header=TRUE)

##Go enrichment analysis
#subset upregulated genes
genes_up <-  res_G1_sig[res_G1_sig$log2FoldChange > 1, "ENSEMBL"]

#GO upregulated
GO_up <- enrichGO(gene = genes_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#plot the results

# plot the results
barplot(GO_up, showCategory = 20)

# Grafiği PNG dosyasına kaydetme
png("GO_results_up_regulated_Group1_20.png", res = 250, width = 1400, height = 3000)
barplot(GO_up, showCategory = 20)
dev.off()

##Go enrichment analysis
#subset upregulated genes
genes_down <-  res_G1_sig[res_G1_sig$log2FoldChange < -1, "ENSEMBL"]

#GO down
GO_down <- enrichGO(gene = genes_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#plot the results

# plot the results
barplot(GO_down, showCategory = 20)

# Grafiği PNG dosyasına kaydetme
png("GO_results_down_regulated_Group1_20.png", res = 250, width = 1400, height = 3000)
barplot(GO_down, showCategory = 20)
dev.off()

genes <- res_G1_sig[abs(res_G1_sig$log2FoldChange)>1,"ENSEMBL"]
GO <- enrichGO(gene = genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#plot the results

# plot the results
barplot(GO, showCategory = 20)

# Grafiği PNG dosyasına kaydetme
png("GO_results_top_genes_Group1_20.png", res = 250, width = 1400, height = 3000)
barplot(GO_down, showCategory = 20)
dev.off()

########################################
#Heatmap
#filter the genes have only logfold change>1 absolute


heat_top <- res_G1_sig[ (abs(res_G1_sig$log2FoldChange) > 1),]
heat_top <- na.omit(heat_top) #(this NA related to symbols)
#save the duplicated rows
duplicated_rows <- duplicated(heat_top$Symbol)
# remove one of the duplicated rows (this duplication related to symbols)
heat_top <- heat_top[!duplicated_rows, ]
#set the fold change with decreasing order
heat_top <- heat_top[order(heat_top$log2FoldChange, decreasing = TRUE),]

#For heatmap we need z_score normalized, we need log2FoldChange 
#and baseMean.We will do heatmap for top25 and bottom25 significantly differentiated genes


#get the normalised data, and scale them
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object

mat<-assay(rlog_out)[rownames(heat_top), rownames(metadata)] #sig genes x samples

# Gruplar için ortalama değerlerin hesaplanması

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
mean_df <- data.frame(row.names = rownames(mat))
for (sample in samples) {
  mean_df[[sample]] <- rowMeans(mat[,sample_columns[[sample]]], na.rm = TRUE)
}

# Display the first few rows of the new dataframe
head(mean_df)

base_mean <- rowMeans(mean_df)
mat.scaled <- t(apply(mean_df, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mean_df)

##########
#select only subset of genes
num_keep <- 10
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

##########
l2_val <- as.matrix(heat_top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(heat_top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = T, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h

png("heatmap_gene_expression_forG1genesbut all groups.png", res = 300, width = 2500, height = 4000)
print(h)
dev.off()

#################

#only for group1 sites heatmap 

heat_top <- res_G1_sig[ (abs(res_G1_sig$log2FoldChange) > 1),]
heat_top <- na.omit(heat_top) #(this NA related to symbols)
#save the duplicated rows
duplicated_rows <- duplicated(heat_top$Symbol)
# remove one of the duplicated rows (this duplication related to symbols)
heat_top <- heat_top[!duplicated_rows, ]
#set the fold change with decreasing order
heat_top <- heat_top[order(heat_top$log2FoldChange, decreasing = TRUE),]

#For heatmap we need z_score normalized, we need log2FoldChange 
#and baseMean.We will do heatmap for top25 and bottom25 significantly differentiated genes


#get the normalised data, and scale them
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object

mat<-assay(rlog_out)[rownames(heat_top), rownames(metadata)] #sig genes x samples

# Gruplar için ortalama değerlerin hesaplanması

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
mean_df <- data.frame(row.names = rownames(mat))
for (sample in samples) {
  mean_df[[sample]] <- rowMeans(mat[,sample_columns[[sample]]], na.rm = TRUE)
}
#exclude group2
Group2 <- c("D01","D02", "D05", "D06", "D08","D10", "D12")

# Exclude the specified columns from mean_df
mean_df_filtered <- dplyr::select(mean_df, -dplyr::one_of(Group2))

# Display the first few rows of the new dataframe
View(mean_df_filtered)

mean_df <- mean_df_filtered

base_mean <- rowMeans(mean_df)
mat.scaled <- t(apply(mean_df, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mean_df)

##########
#select only subset of genes
num_keep <- 10
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

##########
l2_val <- as.matrix(heat_top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(heat_top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = T, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h

png("heatmap_gene_expression_forG1genesonlyG1.png", res = 300, width = 2500, height = 4000)
print(h)
dev.off()

###########################################
#for all genes (bu data frame i WGCNA icin kullan)
results <- results(dds)
results <- as.data.frame(results)
#filtreleme
results_top <- results[(results$padj<0.05)&(results$baseMean>50)&(abs(results$log2FoldChange)>1),]

#ENSEMBL mapping
df.top <- results_top

ids <- row.names(df.top)

# load mapping files
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
# Convert rownames to a new column for both dataframes
df.top$Daphnia_gene <- rownames(df.top)
# Convert rownames to a column
id.maps.matched$Daphnia_gene <- rownames(id.maps.matched)

# Split the ENSEMBL IDs and take only the first one
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$V2, split = ';'), `[`, 1)

# Now, merge the ID map with your DESeq2 results without expanding
merged_data <- merge(df.top, id.maps.matched, by.x = "Daphnia_gene", by.y = "row.names", all.x = TRUE)

# Set the 'Daphnia_gene' column as rownames
row.names(merged_data) <- merged_data$Daphnia_gene

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene <- NULL

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene.y <- NULL

# Also remove the 'V2' column as it's no longer needed
merged_data$V2 <- NULL


library(org.Hs.eg.db)
# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
merged_data$Symbol <- mapIds(org.Hs.eg.db,
                             keys=merged_data$ENSEMBL,
                             keytype="ENSEMBL",
                             column="SYMBOL",
                             multiVals="first")

results_top <- merged_data

results_top <- na.omit(results_top) #(this NA related to symbols)
#save the duplicated rows
duplicated_rows <- duplicated(results_top$Symbol)
# remove one of the duplicated rows (this duplication related to symbols)
results_top <- results_top[!duplicated_rows, ]
#set the fold change with decreasing order
results_top <- results_top[order(results_top$log2FoldChange, decreasing = TRUE),]

heat_top <- results_top

#get the normalised data, and scale them
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object

mat<-assay(rlog_out)[rownames(heat_top), rownames(metadata)] #sig genes x samples

# Gruplar için ortalama değerlerin hesaplanması

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
mean_df <- data.frame(row.names = rownames(mat))
for (sample in samples) {
  mean_df[[sample]] <- rowMeans(mat[,sample_columns[[sample]]], na.rm = TRUE)
}


base_mean <- rowMeans(mean_df)
mat.scaled <- t(apply(mean_df, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mean_df)

##########
#select only subset of genes
num_keep <- 10
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

##########
l2_val <- as.matrix(heat_top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(heat_top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = T, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h

png("heatmap_gene_expression_forallgroups.png", res = 300, width = 3000, height = 4000)
print(h)
dev.off()

######################################
#Results for Group2
#Get the group2 results 
res_G2 <- results(dds, contrast=c("Description","Group2","Control"))
res_G2 <- na.omit(res_G2)

#filter the baseMean above 50 (lower is noise)
res_G2 <- res_G2[res_G2$baseMean>50,]
#convert to df
res_G2_df <- as.data.frame(res_G2)


###########################
#Before volcano plot
#ENSEMBL mapping
df.top <- res_G2_df 

ids <- row.names(df.top)

# load mapping files
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
# Convert rownames to a new column for both dataframes
df.top$Daphnia_gene <- rownames(df.top)
# Convert rownames to a column
id.maps.matched$Daphnia_gene <- rownames(id.maps.matched)

# Split the ENSEMBL IDs and take only the first one
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$V2, split = ';'), `[`, 1)

# Now, merge the ID map with your DESeq2 results without expanding
merged_data <- merge(df.top, id.maps.matched, by.x = "Daphnia_gene", by.y = "row.names", all.x = TRUE)

# Set the 'Daphnia_gene' column as rownames
row.names(merged_data) <- merged_data$Daphnia_gene

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene <- NULL

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene.y <- NULL

# Also remove the 'V2' column as it's no longer needed
merged_data$V2 <- NULL


library(org.Hs.eg.db)
# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
merged_data$Symbol <- mapIds(org.Hs.eg.db,
                             keys=merged_data$ENSEMBL,
                             keytype="ENSEMBL",
                             column="SYMBOL",
                             multiVals="first")

#for saving the data frame (bu charactere donusuyor sonrasinda kullanma)
merged_data_df <- merged_data
merged_data_df$Symbol <- sapply(merged_data_df$Symbol, function(x) paste(x, collapse = ", "))
merged_data_df <- as.data.frame(merged_data_df)
write.csv(merged_data_df,file="deseq_res_G2_sig_genes_1Xwithsmybols_notfilteredwithpadj.csv",row.names = TRUE)

##################
#assign the merged_data to 
res_G2_df <- merged_data
#volcano plot for Group1
library(EnhancedVolcano)
library(textshaping)

# Grafiği çıktı olarak ayarla
png("volcano_plot_Group2_1X.png", width = 10, height = 6, units = "in", res = 300)

# Labels
labels <- res_G2_df$Symbol

# EnhancedVolcano fonksiyonunu çağırma
EnhancedVolcano(res_G2_df, x = "log2FoldChange", y = "padj", lab = labels,
                pCutoff = 1e-2, FCcutoff = 0.8, labSize = 3.0, pointSize = 1.0, xlim = c(-5, 5))

# Çıktıyı kapat
dev.off()
###################
#Filter based on padj<0.05

# First filter the p_Value is lower than 0.05
res_G2_sig <- results(dds, alpha=0.05,contrast=c("Description", "Group2", "Control"))

#check the summary
summary(res_G2_sig)
#121 upregulated
#334 downregulated

#drop the NAs
res_G2_sig <- na.omit(res_G2_sig)
#filter the significant genes based on (padj<0.05)
res_G2_sig <- res_G2_sig[res_G2_sig$padj<0.05,]
#convert to a df
res_G2_sig <- as.data.frame(res_G2_sig)

#Mapping again:

#Before volcano plot
#ENSEMBL mapping
df.top <- res_G2_sig 

ids <- row.names(df.top)

# load mapping files
id.maps <- read.csv("annotations/rna_dma_to_hsa_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
# Convert rownames to a new column for both dataframes
df.top$Daphnia_gene <- rownames(df.top)
# Convert rownames to a column
id.maps.matched$Daphnia_gene <- rownames(id.maps.matched)

# Split the ENSEMBL IDs and take only the first one
id.maps.matched$ENSEMBL <- sapply(strsplit(id.maps.matched$V2, split = ';'), `[`, 1)

# Now, merge the ID map with your DESeq2 results without expanding
merged_data <- merge(df.top, id.maps.matched, by.x = "Daphnia_gene", by.y = "row.names", all.x = TRUE)

# Set the 'Daphnia_gene' column as rownames
row.names(merged_data) <- merged_data$Daphnia_gene

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene <- NULL

# Remove the 'Daphnia_gene' column as it is now the rownames
merged_data$Daphnia_gene.y <- NULL

# Also remove the 'V2' column as it's no longer needed
merged_data$V2 <- NULL


library(org.Hs.eg.db)
# fromNode için ENSEMBL ID'leri gen sembollerine dönüştür
merged_data$Symbol <- mapIds(org.Hs.eg.db,
                             keys=merged_data$ENSEMBL,
                             keytype="ENSEMBL",
                             column="SYMBOL",
                             multiVals="first")
#Assign to res_G1_sig
res_G2_sig <- merged_data
#for saving the data frame (bu charactere donusuyor sonrasinda kullanma)
merged_data_df <- merged_data
merged_data_df$Symbol <- sapply(merged_data_df$Symbol, function(x) paste(x, collapse = ", "))
merged_data_df <- as.data.frame(merged_data_df)
write.csv(merged_data_df,file="deseq_res_G2_sig_genes_1X_withsymbolwithfilterpadj_005.csv",row.names = TRUE)

######################################################
#Reopen again
res_G2_sig <- read.csv("deseq_res_G2_sig_genes_1X_withsymbolwithfilterpadj_005.csv",row.names=1,header=TRUE)

##Go enrichment analysis
#subset upregulated genes
genes_up <-  res_G2_sig[res_G2_sig$log2FoldChange > 2, "ENSEMBL"]

#GO upregulated
GO_up <- enrichGO(gene = genes_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#plot the results

# plot the results
barplot(GO_up, showCategory = 10)

# Grafiği PNG dosyasına kaydetme
png("GO_results_up_regulated_Group2.png", res = 250, width = 1400, height = 1800)
barplot(GO_up, showCategory = 10)
dev.off()

##Go enrichment analysis
#subset upregulated genes
genes_down <-  res_G2_sig[res_G1_sig$log2FoldChange < -1, "ENSEMBL"]

#GO down
GO_down <- enrichGO(gene = genes_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#plot the results

# plot the results
barplot(GO_down, showCategory = 10)

# Grafiği PNG dosyasına kaydetme
png("GO_results_down_regulated_Group2.png", res = 250, width = 1400, height = 1800)
barplot(GO_down, showCategory = 10)
dev.off()

genes <- res_G2_sig[abs(res_G2_sig$log2FoldChange)>1,"ENSEMBL"]
GO <- enrichGO(gene = genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#plot the results

# plot the results
barplot(GO, showCategory = 10)

# Grafiği PNG dosyasına kaydetme
png("GO_results_top_genes_Group2.png", res = 250, width = 1400, height = 1800)
barplot(GO_down, showCategory = 10)
dev.off()

########################################
#Heatmap
#filter the genes have only logfold change>1 absolute

heat_top <- res_G2_sig[ (abs(res_G2_sig$log2FoldChange) > 1),]
heat_top <- na.omit(heat_top) #(this NA related to symbols)
#save the duplicated rows
duplicated_rows <- duplicated(heat_top$Symbol)
# remove one of the duplicated rows (this duplication related to symbols)
heat_top <- heat_top[!duplicated_rows, ]
#set the fold change with decreasing order
heat_top <- heat_top[order(heat_top$log2FoldChange, decreasing = TRUE),]

#For heatmap we need z_score normalized, we need log2FoldChange 
#and baseMean.We will do heatmap for top25 and bottom25 significantly differentiated genes


#get the normalised data, and scale them
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object

mat<-assay(rlog_out)[rownames(heat_top), rownames(metadata)] #sig genes x samples

# Gruplar için ortalama değerlerin hesaplanması

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
mean_df <- data.frame(row.names = rownames(mat))
for (sample in samples) {
  mean_df[[sample]] <- rowMeans(mat[,sample_columns[[sample]]], na.rm = TRUE)
}

# Display the first few rows of the new dataframe
head(mean_df)

base_mean <- rowMeans(mean_df)
mat.scaled <- t(apply(mean_df, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mean_df)

##########
#select only subset of genes
num_keep <- 10
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

##########
l2_val <- as.matrix(heat_top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(heat_top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = T, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h

png("heatmap_gene_expression_forG2genesbut all groups.png", res = 300, width = 2500, height = 4000)
print(h)
dev.off()

#################

#only for group2 sites heatmap 

heat_top <- res_G2_sig[ (abs(res_G2_sig$log2FoldChange) > 1),]
heat_top <- na.omit(heat_top) #(this NA related to symbols)
#save the duplicated rows
duplicated_rows <- duplicated(heat_top$Symbol)
# remove one of the duplicated rows (this duplication related to symbols)
heat_top <- heat_top[!duplicated_rows, ]
#set the fold change with decreasing order
heat_top <- heat_top[order(heat_top$log2FoldChange, decreasing = TRUE),]

#For heatmap we need z_score normalized, we need log2FoldChange 
#and baseMean.We will do heatmap for top25 and bottom25 significantly differentiated genes


#get the normalised data, and scale them
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object

mat<-assay(rlog_out)[rownames(heat_top), rownames(metadata)] #sig genes x samples

# Gruplar için ortalama değerlerin hesaplanması

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
mean_df <- data.frame(row.names = rownames(mat))
for (sample in samples) {
  mean_df[[sample]] <- rowMeans(mat[,sample_columns[[sample]]], na.rm = TRUE)
}
#exclude group1
Group1 <- c("D03", "D04", "D07" , "D09","D11")

# Exclude the specified columns from mean_df
mean_df_filtered <- dplyr::select(mean_df, -dplyr::one_of(Group2))

# Display the first few rows of the new dataframe
View(mean_df_filtered)

mean_df <- mean_df_filtered

base_mean <- rowMeans(mean_df)
mat.scaled <- t(apply(mean_df, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mean_df)

##########
#select only subset of genes
num_keep <- 10
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

##########
l2_val <- as.matrix(heat_top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(heat_top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = T, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = heat_top$Symbol[rows_keep], 
              cluster_rows = T, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h

png("heatmap_gene_expression_forG2genesonlyG2.png", res = 300, width = 2500, height = 4000)
print(h)
dev.off()

#########################
#pathfinder
#For Group1 significant genes
library(pathfindR)
#arrange the data frame
res_G1_sig_p <- res_G1_sig %>%
  dplyr::filter(!is.na(log2FoldChange) & !is.na(padj) & !is.na(Symbol)) %>%
  # Step 2: Remove duplicates based on Symbol
  dplyr::distinct(Symbol, .keep_all = TRUE) %>%
  # Step 3: Rearrange and rename columns
  dplyr::select(Gene_symbol = Symbol, logFC = log2FoldChange, padj)

#res_G1_sig_p <- res_G1_sig_p[ (abs(res_G1_sig_p$logFC) > 1),]

output_df_G1 <- run_pathfindR(res_G1_sig_p)

#For group2
res_G2_sig_p <- res_G2_sig %>%
  dplyr::filter(!is.na(log2FoldChange) & !is.na(padj) & !is.na(Symbol)) %>%
  # Step 2: Remove duplicates based on Symbol
  dplyr::distinct(Symbol, .keep_all = TRUE) %>%
  # Step 3: Rearrange and rename columns
  dplyr::select(Gene_symbol = Symbol, logFC = log2FoldChange, padj)

#res_G1_sig_p <- res_G1_sig_p[ (abs(res_G1_sig_p$logFC) > 1),]

output_df_G2 <- run_pathfindR(res_G2_sig_p)

