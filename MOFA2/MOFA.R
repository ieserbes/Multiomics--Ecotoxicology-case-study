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

dataset$RNA_seq$data.matrix <- filter_by_variance(dataset$RNA_seq$data.matrix, n=5000)

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


# Creating MOFA object from a list of matrices (features as rows, sample as columns)...
MOFAobject <- create_mofa(dataset$combined)

#check the MOFA object
MOFAobject

## Visualise the samples with values/ missing, across the whole cohort
plot_data_overview(MOFAobject)

## Get the options for preparing the MOFA run
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)

model_opts

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

#Association analysis

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("REF","Site","Description"), 
                                  plot="log_pval"
)

#Plot feature weights

plot_weights(MOFAobject,
             view = "RNA_seq",
             factor = 1,
             nfeatures = 5,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

#get weights
weights <- get_weights(MOFAobject, 
                       views = "all", 
                       factors = "all", 
                       as.data.frame = TRUE 
)
head(weights)


#get factors
factors <- get_factors(MOFAobject, 
                       factors = "all", 
                       as.data.frame = TRUE
)
head(factors)



