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
  library(MOFA2)
  library(MOFAdata)
  library(data.table)
})

###############################################

# Meta veriyi yükle
meta <- read.csv('sample_sheet.csv', row.names = 1, header = TRUE, sep = ',')

# Grupları tanımla
Group1 <- c("D03", "D04", "D07" , "D09","D11")
Group2 <- c("D01","D02", "D05", "D06", "D08","D10", "D12")

# Meta veride Description'ı yeniden düzenle
meta$Description <- case_when(
  meta$Site %in% Group1 ~ "Group1",
  meta$Site %in% Group2 ~ "Group2",
  TRUE ~ as.character(meta$Description) # Varsayılan olarak mevcut Description değerini koru
)

# Hariç tutulacak desenleri tanımla
patterns_to_exclude <- c("B","D06A2","D07A3")

# Meta veride ve diğer veri setlerinde "B" içeren satırları hariç tut
meta <- meta[!grepl(paste(patterns_to_exclude, collapse = "|"), row.names(meta)), ]

# Veri setlerini yükleme fonksiyonu
load.dataset <- function(data.file, data.sep = ',', meta) {
  data <- read.table(data.file, header = TRUE, sep = data.sep, row.names = 1, check.names = FALSE)
  data <- as.matrix(data)
  
  # "B" içeren satırları hariç tut
  data <- data[!grepl(paste(patterns_to_exclude, collapse = "|"), row.names(data)), ]
  
  # Meta veri ve ana veri arasındaki eşleşmeyi bul
  mids <- match(row.names(meta), colnames(data))
  data <- t(data[,mids[!is.na(mids)]])
  matched_meta <- meta[!is.na(mids),]
  
  stopifnot(row.names(matched_meta) == row.names(data))
  
  return(list(meta.data = matched_meta, data.matrix = data))
}

# Her bir veri setini yükle ve meta veriyi ilişkilendir
dataset <- list()
dataset$RNA_seq <- load.dataset(data.file = 'rna_vst_counts.csv', meta = meta)
dataset$RNA_seq$data.matrix <- t(dataset$RNA_seq$data.matrix)
dataset$polar_neg <- load.dataset(data.file = 'polar_neg_pqn_imputed_glog.csv', meta = meta)
dataset$polar_neg$data.matrix <- t(dataset$polar_neg$data.matrix)
dataset$polar_pos <- load.dataset(data.file = 'polar_pos_pqn_imputed_glog.csv', meta = meta)
dataset$polar_pos$data.matrix <- t(dataset$polar_pos$data.matrix)


################################################
#Filter the RNA_seq data by variance

dataset$RNA_seq$data.matrix <- filter_by_variance(dataset$RNA_seq$data.matrix, n=7000)

#Ortak bir veriseti olustur

## Ortak örnekleri bul ve birleştir
# Sütun isimlerine (colnames) göre kesişim bulma
shared.cols <- Reduce(intersect, list(
  colnames(dataset$RNA_seq$data.matrix),
  colnames(dataset$polar_neg$data.matrix),
  colnames(dataset$polar_pos$data.matrix)
))

# Ortak sütunlarla birleştirilmiş bir liste oluşturma
dataset$combined <- list(
  RNA_seq =   dataset$RNA_seq$data.matrix[, shared.cols],
  polar_neg = dataset$polar_neg$data.matrix[, shared.cols],
  polar_pos = dataset$polar_pos$data.matrix[, shared.cols]
)
#######################################
#Annotation RNA-seq
RNA_seq <- dataset$combined$RNA_seq

genes <- row.names(RNA_seq)
df.top <- data.frame(Genes = genes)
#ENSEMBL mapping


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


RNA_seq <- as.data.frame(RNA_seq)
RNA_seq$ENSEMBL <- df.top$ENSEMBL

# ENSEMBL sütununda NA olmayan satırları filtrele
RNA_seq <- RNA_seq[!is.na(RNA_seq$ENSEMBL), ]

# ENSEMBL sütununda tekrar eden değerleri bul (ilk tekrar hariç hepsi TRUE döner)
duplicates <- duplicated(RNA_seq$ENSEMBL)

# Tekrar eden satırları sil
RNA_seq <- RNA_seq[!duplicates, ]

# ENSEMBL sütunundaki değerleri satır isimleri olarak ata
rownames(RNA_seq) <- RNA_seq$ENSEMBL

# ENSEMBL sütununu sil
RNA_seq$ENSEMBL <- NULL
#reassign it
dataset$combined$RNA_seq <- as.matrix(RNA_seq)
################################################################################################################
#Annotation polar_neg
polar_neg <- dataset$combined$polar_neg

peak_id <- row.names(polar_neg)
df.top <- data.frame(peak_id = peak_id)
#ENSEMBL mapping


ids <- df.top$peak_id

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/polar_neg_pkl_to_kegg_annotations.tsv", header = F, sep = '\t', stringsAsFactors = F)
id.maps$V2 <- NULL
colnames(id.maps) <- c("peak", "KEGG_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, peak %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$KEGG <- sapply(strsplit(id.maps.matched$KEGG_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "peak_id", by.y = "peak", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("peak", "KEGG_ids"))]


polar_neg <- as.data.frame(polar_neg)
polar_neg$KEGG <- df.top$KEGG

# ENSEMBL sütununda NA olmayan satırları filtrele
polar_neg <- polar_neg[!is.na(polar_neg$KEGG), ]

# ENSEMBL sütununda tekrar eden değerleri bul (ilk tekrar hariç hepsi TRUE döner)
duplicates <- duplicated(polar_neg$KEGG)

# Tekrar eden satırları sil
polar_neg <- polar_neg[!duplicates, ]

# ENSEMBL sütunundaki değerleri satır isimleri olarak ata
rownames(polar_neg) <- polar_neg$KEGG

# ENSEMBL sütununu sil
polar_neg$KEGG <- NULL
#reassign it
dataset$combined$polar_neg <- as.matrix(polar_neg)
##############################################################################################
#Annotation polar_pos
polar_pos <- dataset$combined$polar_pos

peak_id <- row.names(polar_pos)
df.top <- data.frame(peak_id = peak_id)
#ENSEMBL mapping


ids <- df.top$peak_id

# Mapping dosyasını yükle
id.maps <- read.csv("annotations/polar_pos_pkl_to_kegg_annotations.tsv", header = F, sep = '\t', stringsAsFactors = F)
id.maps$V2 <- NULL
colnames(id.maps) <- c("peak", "KEGG_ids")  # Kolon isimlerini anlamlı hale getir

# Eşleme işlemi
id.maps.matched <- subset(id.maps, peak %in% ids)

# ENSEMBL ID'lerini ayrıştır ve sadece ilkini al
id.maps.matched$KEGG <- sapply(strsplit(id.maps.matched$KEGG_ids, split = ';'), `[`, 1)

# Eşleşen bilgileri 'df.top' ile birleştir
df.top <- merge(df.top, id.maps.matched, by.x = "peak_id", by.y = "peak", all.x = TRUE)

# Gereksiz sütunları kaldır
df.top <- df.top[, !(names(df.top) %in% c("peak", "KEGG_ids"))]


polar_pos <- as.data.frame(polar_pos)
polar_pos$KEGG <- df.top$KEGG

# ENSEMBL sütununda NA olmayan satırları filtrele
polar_pos <- polar_pos[!is.na(polar_pos$KEGG), ]

# ENSEMBL sütununda tekrar eden değerleri bul (ilk tekrar hariç hepsi TRUE döner)
duplicates <- duplicated(polar_pos$KEGG)

# Tekrar eden satırları sil
polar_pos <- polar_pos[!duplicates, ]

# ENSEMBL sütunundaki değerleri satır isimleri olarak ata
rownames(polar_pos) <- polar_pos$KEGG

# ENSEMBL sütununu sil
polar_pos$KEGG <- NULL
#reassign it
dataset$combined$polar_pos <- as.matrix(polar_pos)

##########################################################################
# Creating MOFA object from a list of matrices (features as rows, sample as columns)...
MOFAobject <- create_mofa(dataset$combined)

#check the MOFA object
MOFAobject

## Visualise the samples with values/ missing, across the whole cohort
plot_data_overview(MOFAobject)

## Get the options for preparing the MOFA run
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)

model_opts$num_factors <- 10

train_opts <- get_default_training_options(MOFAobject)
## Set the convergence mode to 'slow', for more sensitive MOFA run
train_opts$convergence_mode <- "slow"

train_opts

#prepare the MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

#Run the MOFA object
MOFAobject <- run_mofa(MOFAobject, outfile="MOFA2_symbols_dataset.hdf5", save_data=TRUE, use_basilisk = TRUE)

## Save the computed model for future reference
saveRDS(MOFAobject,file="MOFA2_symbols_Computed.rds")

## check if all metadata sample names are present in MOFAobject
all(sort(row.names(meta))==sort(unlist(samples_names(MOFAobject))))

#The MOFA object consists of multiple slots where relevant data and information is stored. 
slotNames(MOFAobject)

#how to reach the slots
names(MOFAobject@data)

#reach the RNA-seq data
View(MOFAobject@data$RNA_seq$group1)
dim(MOFAobject@data$RNA_seq$group1)


#create a sample kolon
meta$sample<-rownames(meta)
rownames(meta) <- NULL
meta <- dplyr::select(meta, sample, REF, Site, Description)

# Add sample metadata to the model
samples_metadata(MOFAobject) <- meta

#Correlation between factors
plot_factor_cor(MOFAobject)


#Variance decomposition by Factor
plot_variance_explained(MOFAobject, max_r2=15)

#Total variance explained per view
"A reasonable question is whether the model is providing a good fit to the data.
For this we can plot the total variance explained (using all factors). The 
resulting values will depend on the nature of the data set"


plot_variance_explained(MOFAobject, plot_total = T)[[2]]


##Association analysis

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("REF","Site","Description"), 
                                  plot="log_pval"
)

#characterization of Factor2

#Plot factor values

plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Factor2"
)
#RNA_seq

#Plot weights-RNA-seq
plot_weights(MOFAobject,
             view = "RNA_seq",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

#distribution of weights
plot_top_weights(MOFAobject,
                 view = "RNA_seq",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

#CES4A has a positive weight. 

plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Description",
            dodge = TRUE,
            add_violin = TRUE
)

#expression plot
f <- c("ENSG00000172824","ENSG00000189058","ENSG00000100652","ENSG00000101670")
plot_data_scatter(MOFAobject, 
                  view = "RNA_seq",
                  factor = 2,  
                  features = f,
                  sign = "negative",
                  color_by = "Description"
) + labs(y="RNA expression")


#polar negative

#Plot factor values

#Plot weights-polar negative
plot_weights(MOFAobject,
             view = "polar_neg",
             factor = 2,
             nfeatures = 8,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

#distribution of weights
plot_top_weights(MOFAobject,
                 view = "polar_neg",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

#plot by groups

plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Description",
            dodge = TRUE,
            add_violin = TRUE
)

#expression plot
f <- c("ENSG00000172824","ENSG00000189058","ENSG00000100652","ENSG00000101670")
plot_data_scatter(MOFAobject, 
                  view = "polar_neg",
                  factor = 2,  
                  features = 4,
                  sign = "negative",
                  color_by = "Description"
) + labs(y="Polar Negative Metabolite Concentrations")


#polar positive

#Plot factor values

#Plot weights-polar negative
plot_weights(MOFAobject,
             view = "polar_pos",
             factor = 2,
             nfeatures = 8,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

#distribution of weights
plot_top_weights(MOFAobject,
                 view = "polar_pos",
                 factor = 2,
                 nfeatures = 6,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

#plot by groups

plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Description",
            dodge = TRUE,
            add_violin = TRUE
)

#expression plot
f <- c("ENSG00000172824","ENSG00000189058","ENSG00000100652","ENSG00000101670")
plot_data_scatter(MOFAobject, 
                  view = "polar_pos",
                  factor = 2,  
                  features = 4,
                  sign = "negative",
                  color_by = "Description"
) + labs(y="Polar Positive Metabolite Concentrations")

############################################################################################

#RNA_seq-factor5

#Plot factor values

plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Description",
            dodge = TRUE,
            add_violin = TRUE
)


plot_factor(MOFAobject, 
            factors = 5, 
            color_by = "Description",
            dodge = TRUE,
            add_violin = TRUE
)


#Plot weights-RNA-seq
plot_weights(MOFAobject,
             view = "RNA_seq",
             factor = 5,
             nfeatures = 6,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

#distribution of weights
plot_top_weights(MOFAobject,
                 view = "RNA_seq",
                 factor = 5,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

#CES4A has a positive weight. 

plot_factor(MOFAobject, 
            factors = 5, 
            color_by = "Description",
            dodge = TRUE,
            add_violin = TRUE
)

#expression plot
f <- c("ENSG00000142973","ENSG00000064763","ENSG00000116833","ENSG00000168878")

plot_data_scatter(MOFAobject, 
                  view = "RNA_seq",
                  factor = 5,  
                  features = f,
                  sign = "negative",
                  color_by = "Description"
) + labs(y="RNA expression")





#Enrichment analysis

utils::data(reactomeGS)


# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "RNA_seq",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "RNA_seq",
                               sign = "negative"
)

#plot_enrichment_heatmap(res.positive)-general for all factors

plot_enrichment_heatmap(res.positive)

#plot_enrichment_heatmap(res.positive)-general for all factors
plot_enrichment_heatmap(res.negative)


#plot_enrichment-positive
plot_enrichment(res.positive, factor = 1, max.pathways = 5,alpha = 1)
plot_enrichment(res.positive, factor = 2, max.pathways = 5,alpha = 1)
plot_enrichment(res.positive, factor = 5, max.pathways = 5,alpha = 1)


#plot_enrichment-negative

plot_enrichment(res.negative, factor = 1, max.pathways = 5,alpha = 1)
plot_enrichment(res.negative, factor = 2, max.pathways = 5,alpha = 1) #not
plot_enrichment(res.negative, factor = 5, max.pathways = 5,alpha = 1)

###########################################################################
utils::data(MSigDB_v6.0_C5_human)




# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = MSigDB_v6.0_C5_human, 
                               view = "RNA_seq",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = MSigDB_v6.0_C5_human, 
                               view = "RNA_seq",
                               sign = "negative"
)

#plot_enrichment_heatmap(res.positive)-general for all factors

plot_enrichment_heatmap(res.positive)

#plot_enrichment_heatmap(res.positive)-general for all factors
plot_enrichment_heatmap(res.negative)


#plot_enrichment-positive
plot_enrichment(res.positive, factor = 1, max.pathways = 5,alpha = 1)
plot_enrichment(res.positive, factor = 2, max.pathways = 10,alpha = 1)
plot_enrichment(res.positive, factor = 5, max.pathways = 5,alpha = 1)


#check which genes driving this enrichment
plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 2, 
  max.pathways = 8
)

#subset the genes 
# 'GO_RESPONSE_TO_DRUG' satırını bul
satir_index <- which(rownames(res.positive$feature.sets) == "GO_RESPONSE_TO_DRUG")

# 'GO_RESPONSE_TO_DRUG' satırında değeri '1' olan sütun adlarını bul
drug_pos <- colnames(res.positive$feature.sets)[res.positive$feature.sets[satir_index, ] == 1]
drug_pos_2 <- drug_pos[1:4]

plot_data_scatter(MOFAobject, 
                  view = "RNA_seq",
                  factor = 2,  
                  features = drug_pos,
                  sign = "positive",
                  color_by = "Description"
) + labs(y="RNA expression")


#plot_enrichment-negative

plot_enrichment(res.negative, factor = 1, max.pathways = 5,alpha = 1)
plot_enrichment(res.negative, factor = 2, max.pathways = 5,alpha = 1) #not
plot_enrichment(res.negative, factor = 5, max.pathways = 10,alpha = 1)

plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 5, 
  max.pathways = 3
)




#subset the genes 
# 'GO_ORGANIC_HYDROXY_COMPOUND_METABOLIC_PROCESS' satırını bul



satir_index <- which(rownames(res.negative$feature.sets) == "GO_ORGANIC_HYDROXY_COMPOUND_METABOLIC_PROCESS")

# 'GO_RESPONSE_TO_DRUG' satırında değeri '1' olan sütun adlarını bul
organ_neg <- colnames(res.negative$feature.sets)[res.negative$feature.sets[satir_index, ] == 1]
organ_2 <- organ_neg[1:4]

plot_data_scatter(MOFAobject, 
                  view = "RNA_seq",
                  factor = 5,  
                  features = organ_2,
                  sign = "positive",
                  color_by = "Description"
) + labs(y="RNA expression")

