rm(list = ls())

#Set the working directory
setwd("D:/Academic2024/MSc_Bioinformatics/Modules/Module-10-Computational Biology for Complex Systems-DL BIO - 37415/Data/Data2")


#########################################
#Required packages
suppressMessages({
  library(dplyr)
  library(iRF)
  library(AUC)
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

source('scripts/fileio.R')
ds <- load.dataset(
  meta.file = 'sample_sheet.csv', meta.sep = ',',
  data.file = 'polar_neg_pqn_imputed_glog.csv', data.sep = ','
)

data <- ds$data.matrix
meta <- ds$meta.data

patterns_to_exclude <- c("B") #Here I can exclude more things if I want

# Filter rows from the data matrix 
data <- data[!grepl(paste(patterns_to_exclude, collapse = "|"), rownames(data)), ]

#assign the meta_rna
#Assign the groups
Group1 <- c("D03", "D04", "D07" , "D09","D11")
Group2 <- c("D01","D02", "D05", "D06", "D08","D10", "D12")
row_names <- row.names(meta)
meta <- meta %>%
  mutate(Description = case_when(
    Site %in% Group1 ~ "Group1",
    Site %in% Group2 ~ "Group2",
    TRUE ~ Description # Eğer Site Group1 veya Group2 ile eşleşmiyorsa mevcut Description değerini koru
  ))


# Rewrite the rownames
#rownames(data) <- row_names

#Filter the metadata
meta <- meta[!grepl(paste(patterns_to_exclude, collapse = "|"), rownames(meta)), ]

##################
#Prepare inputs
X <- data
Y <-  meta$Description
Y <- as.factor(Y)
Y <- as.factor(as.numeric(relevel(Y, 'Control'))-1)


n <- dim(X)[1]
p <- dim(X)[2]

# split training-testing set
train.id <- sample(1:n, size = 8*round(n/10))
test.id <- setdiff(1:n, train.id)

# fit iRF without iterations
fit <- iRF(x = X[train.id,], 
           y = Y[train.id], 
           xtest = X[test.id,], 
           ytest = Y[test.id],
           n.iter = 5, 
           iter.return = 1:5,
           n.core = 4
)


# plot ROC/AUC
rf <- fit$rf.list

png("ROC_Curve.png", width = 800, height = 600) # PNG olarak kaydetmek için
# veya jpeg("ROC_Curve.jpg", width = 800, height = 600) # JPEG olarak kaydetmek için
# veya pdf("ROC_Curve.pdf", width = 8, height = 6) # PDF olarak kaydetmek için
# veya svg("ROC_Curve.svg", width = 8, height = 6) # SVG olarak kaydetmek için

plot(0:1, 0:1, type = 'n', xlim = c(0, 1), ylim = c(0, 1), xlab = 'FPR', ylab = 'TPR', main = 'ROC Curve')

for (iter in 1:5){
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(rf[[iter]]$test$votes[,2], Y[test.id])
  lines(roc.info$fpr, roc.info$tpr, type = 'l', col = iter, lwd = 2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep = ''))
}

legend('bottomright', legend=paste('iter:', 1:5), col=1:5, lwd=2, bty='n')

dev.off() # Çizim cihazını kapat
##############################################



# plot feature weights
par(mfrow = c(1,5))
for (iter in 1:5){
  varImpPlot(rf[[iter]], n.var = 10, main = paste('Variable Importance (iter:', iter, ')')) 
}


######################


# fit iRF with iterations
fit <- iRF(x = X[train.id,], 
           y = Y[train.id], 
           xtest = X[test.id,], 
           ytest = Y[test.id],
           n.iter = 5, 
           n.bootstrap = 50, 
           select.iter = T,
           n.core = 4
)


# plot iteration stability scores 
toplot <- fit$interaction$stability
names(toplot) <- fit$interaction$int
toplot <- sort(toplot, decreasing = TRUE)

# Grafiği çiz
dotchart(rev(toplot[1:min(20, length(toplot))]), 
         xlab = 'Stability Score', 
         xlim = c(0, 1)
)

###################################
#Exploring the 313_0985112
vitamin_b5 <- data.frame(data[, 677])
colnames(vitamin_b5) <- "vitamin_b5"
vitamin_b5$Description <- meta$Description

ggplot(vitamin_b5, aes(x = Description, y = vitamin_b5, fill = Description)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Vitamin B5 Levels by Group",
       x = "Group",
       y = "Vitamin B5 Concentration")



# t-test
# Control ve Group2 grubunu seçmek için:
control_group <- vitamin_b5[vitamin_b5$Description == "Control", "vitamin_b5"]
group2 <- vitamin_b5[vitamin_b5$Description == "Group2", "vitamin_b5"]

# Varyans homojenliği testi
var.test(control_group, group2)

# Eşit varyanslar için t-testi
t.test(control_group, group2, var.equal = FALSE)

# Control ve Group1 grubunu seçmek için:

group1 <- vitamin_b5[vitamin_b5$Description == "Group1", "vitamin_b5"]

# Varyans homojenliği testi
var.test(control_group, group1)

# Eşit varyanslar için t-testi
t.test(control_group, group1, var.equal = FALSE)


