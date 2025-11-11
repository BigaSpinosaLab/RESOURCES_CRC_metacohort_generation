#------------------------------------------------------------------------------#                            
 
                    ###       Complete Metacohort       ###
 
#------------------------------------------------------------------------------#

# Used datasets:

# GSE14333:     Jorissen RN, Gibbs P, Christie M, Prakash S et al. 
#               Metastasis-Associated Gene Expression Changes 
#               Predict Poor Outcomes in Patients with Dukes 
#               Stage B and C Colorectal Cancer. Clin Cancer 
#               Res 2009 Dec 15;15(24):7642-7651. 
#               PMID: 19996206
#
# GSE143985:    Shinto E, Yoshida Y, Kajiwara Y, Okamoto K et al. 
#               Clinical Significance of a Gene Signature 
#               Generated from Tumor Budding Grade in Colon 
#               Cancer. Ann Surg Oncol 2020 Oct;27(10):4044-4054. 
#               PMID: 32328985
#
# GSE17536:     Smith JJ, Deane NG, Wu F, Merchant NB et al. 
# GSE17537:     Experimentally derived metastasis gene 
#               expression profile predicts recurrence 
#               and death in patients with colon cancer. 
#               Gastroenterology 2010 Mar;138(3):958-68. 
#               PMID: 19914252
#
# GSE33114:    	de Sousa E Melo F, Colak S, Buikhuisen J, Koster J et al. 
#               Methylation of cancer-stem-cell-associated Wnt target genes 
#               predicts poor prognosis in colorectal cancer patients. Cell 
#               Stem Cell 2011 Nov 4;9(5):476-85. 
#               PMID: 22056143
#
# GSE38832:     Tripathi MK, Deane NG, Zhu J, An H et al. 
#               Nuclear factor of activated T-cell activity is associated 
#               with metastatic capacity in colon cancer. 
#               Cancer Res 2014 Dec 1;74(23):6947-57. 
#               PMID: 25320007
#
# GSE39582:     Marisa L, de Reyniès A, Duval A, Selves J et al. Gene 
#               expression classification of colon cancer into molecular 
#               subtypes: characterization, validation, and prognostic value. 
#               PLoS Med 2013;10(5):e1001453. 
#               PMID: 23700391


#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
########################        Install libraries        #######################
#------------------------------------------------------------------------------#

library(hgu133plus2.db) # Annotation
library(GEOquery) # Get CEL files
library(affy) # RMA
library(jetset) # Jetset Score
library(sva) # ComBat
library(frma)
library(hgu133plus2frmavecs) 
library(ggplot2)
library(CMSclassifier)
library(org.Hs.eg.db)


#------------------------------------------------------------------------------#
########################        Obtain .CEL files        #######################
#------------------------------------------------------------------------------#


setwd("/CEL_Files")

data_names <- c("GSE14333", "GSE143985", "GSE17536", "GSE17537",
                "GSE33114", "GSE38832", "GSE39582")


# Download raw data CEL files for each project

# options(timeout = max(300, getOption("timeout")))  ## Use it in console if
# options(download.file.method.GEOquery = "wget")    ## "failed to download"

lapply(data_names, function(x){
  x = getGEOSuppFiles(x)
})


### "Pheno_data_filt.txt" file (1118 samples) --> "clinical_data_extraction.R"

pheno.data <- readRDS("pheno_data_august2025.rds")
dd <- readRDS("date_GSM.rds")

setwd("/CEL_Files")
celfiles <-list.files()

## Only Useful Data
celfiles.GSM <- gsub("_.*", "", celfiles)
celfiles.GSM <-gsub("\\..*", "", celfiles.GSM)
celfiles <- celfiles[celfiles.GSM %in% rownames(pheno.data)]

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
########################          fRMA (Dataset)         #######################
#------------------------------------------------------------------------------#

### fRMA (need dates/GSM from "date_extraction.R")
data_f.eset <- data.frame("Omit" = rep(NA, 54676))
setwd("/CEL_Files")

# Select x dataset
for (dataset in unique(dd$dataset)) {
  
  # Select celfiles that were analyzed on x dataset
  GSM <- dd[dd$dataset %in% dataset,]
  GSM <- GSM$GSM
  
  temp.celfiles.GSM <- gsub("_.*", "", celfiles)
  temp.celfiles.GSM <-gsub("\\..*", "", temp.celfiles.GSM)
  
  temp.celfiles <- celfiles[temp.celfiles.GSM %in% GSM] 
  pheno.batch <- pheno.data[rownames(pheno.data) %in% GSM,]
  rownames(pheno.batch) <- temp.celfiles
  pheno.batch$Row.names = NULL
  batch <- ReadAffy(filenames = temp.celfiles)
  
  frma <- frma(batch)
  
  data_eset <- as.data.frame(data.frame(frma))
  data_eset_probes <- colnames(data_eset)
  data_eset <- t(data_eset)
  rownames(data_eset) <- data_eset_probes
  
  colnames(data_eset) <- GSM
  data_f.eset <- cbind(data_f.eset, data_eset)
  
}

# Row 54676 is removed because does not contain expression data 
# (it is a "sample" row generated by fRMA)
# Column 1 is not necessary either, and therefore is also removed

data_eset <- data_f.eset[-54676,-1]

rownames2 <- lapply(rownames(data_eset), function(x) {
  if (startsWith(x, "X")) {
    x <- sub('.', '', x)
  }
  x
})

rownames <- rownames(data_eset)
rownames(data_eset) <- rownames2


#------------------------------------------------------------------------------#
###########################   Jetset Scoring   #################################
#------------------------------------------------------------------------------#

## ANNOTATE GENES
Annot <- toTable(hgu133plus2SYMBOL)
rownames(Annot) <- Annot$probe_id

allscores <- jscores('hgu133plus2')

allscores <- allscores[rownames(allscores) %in% Annot$probe_id,]
allscores_only_scores <- as.data.frame(allscores[,7])
rownames(allscores_only_scores) <- rownames(allscores)
colnames(allscores_only_scores) <- "scores"


Annot_score <- merge(Annot, allscores_only_scores, by = 'row.names', all = TRUE)
rownames(Annot_score) <- Annot_score$probe_id
Annot_score$Row.names <- NULL

data_eset_anot <- merge(Annot_score, data_eset, by.x=0, by.y=0, all.y=T)
row.names(data_eset_anot) <- data_eset_anot$Row.names
data_eset_anot$Row.names <- NULL
data_eset_anot$probe_id <- rownames(data_eset_anot)


#------------------------------------------------------------------------------#
##################              PCA (Pre-Combat)             ###################
#------------------------------------------------------------------------------#

only.exprs <- data_eset_anot[,-c(1:3)]
dataset <- data.frame("dataset" = pheno.data$first_author)
rownames(dataset) <- row.names(pheno.data)

pca = prcomp(t(only.exprs), scale=TRUE, center=TRUE)

PCs <- as.data.frame(pca$x)
PCs$Sample <- rownames(PCs)
PCs <- merge(PCs, dataset, by=0)

Variance <- round(summary(pca)$importance[2,]*100, digits=1)

plot2 <- ggplot(PCs, aes(PC1, 
                         PC2, 
                         color=dataset))

plot2 <- plot2 +
  geom_point(size=1,alpha=0.8) +
  xlab(paste0("PC1", ": ", Variance[1], "% variance")) +
  ylab(paste0("PC2", ": ", Variance[2], "% variance")) +
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title=element_text(size=10,face="bold"),
        title=element_text(size=12,face="bold"),
        axis.text.x = element_text(size=10, hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1))+
  ggtitle("PCA All datasets (fRMA(Dataset))") 

xlim_vector2 = c(1.2*min(PCs$PC1), 1.2*max(PCs$PC1))

ylim_vector2 = c(1.2*min(PCs$PC2), 1.2*max(PCs$PC2))

plot2 <- plot2 +
  xlim(xlim_vector2) +
  ylim(ylim_vector2)

pdf("/FINAL_PLOTS/PCA_All_datasets_fRMA_Dataset.pdf", 
    height = 5, width = 6)
plot(plot2)
dev.off()


#------------------------------------------------------------------------------#
##################   COMBAT (Remove database batch effect)   ###################
#------------------------------------------------------------------------------#

combat.data <- ComBat(dat = data_eset,
                      batch = pheno.data$dataset)

exprs.data <- merge(Annot_score, combat.data, by.x=0, by.y=0, all.y=T)
row.names(exprs.data) <- exprs.data$Row.names
exprs.data$Row.names <- NULL
exprs.data$probe_id <- rownames(exprs.data)


#------------------------------------------------------------------------------#
##################             PCA (Post-Combat)             ###################
#------------------------------------------------------------------------------#

only.exprs <- exprs.data[,-c(1:3)]
dataset <- data.frame("dataset" = pheno.data$first_author)
rownames(dataset) <- row.names(pheno.data)

pca = prcomp(t(only.exprs), scale=TRUE, center=TRUE)

PCs <- as.data.frame(pca$x)
PCs$Sample <- rownames(PCs)
PCs <- merge(PCs, dataset, by=0)

Variance <- round(summary(pca)$importance[2,]*100, digits=1)

plot2 <- ggplot(PCs, aes(PC1, 
                         PC2, 
                         color=dataset))

plot2 <- plot2 +
  geom_point(size=1,alpha=0.8) +
  xlab(paste0("PC1", ": ", Variance[1], "% variance")) +
  ylab(paste0("PC2", ": ", Variance[2], "% variance")) +
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title=element_text(size=10,face="bold"),
        title=element_text(size=12,face="bold"),
        axis.text.x = element_text(size=10, hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1))+
  ggtitle("PCA All datasets (fRMA(Dataset))") 

xlim_vector2 = c(1.2*min(PCs$PC1), 1.2*max(PCs$PC1))

ylim_vector2 = c(1.2*min(PCs$PC2), 1.2*max(PCs$PC2))

plot2 <- plot2 +
  xlim(xlim_vector2) +
  ylim(ylim_vector2)

pdf("/FINAL_PLOTS/PCA_All_datasets_fRMA_Dataset_ComBat.pdf", 
    height = 5, width = 6)
plot(plot2)
dev.off()


#------------------------------------------------------------------------------#
###############################  CMS classifier ################################
#------------------------------------------------------------------------------#


exprs.data.mean <- cbind(exprs.data[c(1:3)], "mean" = rowMeans(exprs.data[,-c(1:3)]))
exprs.data.mean <- exprs.data.mean %>%
  group_by(symbol) %>%
  dplyr::slice(which.max(mean)) 
exprs.data.selected <- exprs.data[exprs.data$probe_id %in% exprs.data.mean$probe_id,]
exprs.data.selected <- exprs.data.selected[which(!is.na(exprs.data.selected$symbol)),]

# FROM SYMBOL TO ENTREZ
hs <- org.Hs.eg.db
my.symbols <- exprs.data.selected$symbol
symbol_entrez <- select(hs, 
                        keys = my.symbols,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")
colnames(symbol_entrez) <- c("symbol", "entrezid")

exprs.data.selected <- merge(symbol_entrez, exprs.data.selected, by ="symbol")
rownames(exprs.data.selected) <- exprs.data.selected$entrezid
exprs.data.selected <- exprs.data.selected[,-c(1:4)]

#CMS
set.seed(123)
Rfcms <- CMSclassifier::classifyCMS(exprs.data.selected,method="RF")[[3]]
SScms <- CMSclassifier::classifyCMS(exprs.data.selected,method="SSP")[[3]]

finalModel
# Call:
#  randomForest(x = x[idxs, ], y = y[idxs], ntree = 1L, importance = TRUE) 
#                Type of random forest: classification
#                      Number of trees: 100
# No. of variables tried at each split: 16

pheno.data$cms <- Rfcms$RF.predictedCMS
table(pheno.data$cms)

#------------------------------------------------------------------------------#
###############################   Save .RData   ################################
#------------------------------------------------------------------------------#
 
save(exprs.data, pheno.data, file = "CRC_metacohort_August_2025.RData")

