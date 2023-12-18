#Libraries ----
library(readxl)
library(vioplot)
library(RankProd)
library(reshape2)
library(readr)
library(penalized)
library(dplyr)
library("GSVA")
library("msigdbr")
library("RankProd")
library(FactoMineR)
library(factoextra)
library(ComplexHeatmap)
library(circlize)

#loading data ----
#miRNA
agilentmiRNA <- read_excel("AgilentData/RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.xls", sheet = "Results", skip = 10)
  #cleaning up miRNa data
  rows <- agilentmiRNA$`Identifier c`
  agilentmiRNA <- agilentmiRNA[,-c(1:9)]#removing unneeded columns
  rownames(agilentmiRNA) <- rows
#mRNA
agilentmRNA <- read_excel("AgilentData/RNA__Agilent_mRNA_log2.xls", skip = 10)
  #cleaning up mRNA data
  agilentmRNA <- agilentmRNA[-which(agilentmRNA$`Entrez gene id e` == '0'),] #removing genes with no gene id
  agilentmRNAdupRemoved <- agilentmRNA[!duplicated(agilentmRNA$`Gene name d`), ]
    rows <- agilentmRNAdupRemoved$`Gene name d`
    agilentmRNAdupRemoved <- agilentmRNAdupRemoved[,-c(1:7)]
    rownames(agilentmRNAdupRemoved) <- rows
  rows <- agilentmRNA$`Identifier c`
  agilentmRNA <- agilentmRNA[,-c(1, 3:7)]#removing unneeded columns
  rownames(agilentmRNA) <- rows
  egfrData <- agilentmRNA[which(agilentmRNA$`Gene name d` == 'EGFR'),2:ncol(agilentmRNA)]
#Meta Data
  metaData <- read_excel("AgilentData/NCI60_CELL_LINE_METADATA.xls", skip = 7)
  metaData <- metaData[1:60,] #removing the bottom sheet excels
  rownames(metaData) <- colnames(agilentmiRNA)
#copy number
  cnvData <- read_excel("AgilentData/DNA__aCGH_Agilent_44K_Copy_number_estimate.xls", skip = 10)
  cnvData <- cnvData[cnvData$`Gene name d` == 'EGFR',-c(1:6)] #selecting for egfr data
#Visualization ----
  #PCA mrna data plot
  mrnaData <- PCA(t(agilentmRNA[,2:ncol(agilentmRNA)]),scale.unit = TRUE, graph = FALSE)
  fviz_pca_ind(mrnaData, habillage = as.factor(metaData$`tissue of origin a`), palette = "npg") + labs(title ="NC1-60 mrna data PCA plot", x = "PC1", y = "PC2")
  #PCA mirna data plot
  mirnaData <- PCA(t(agilentmiRNA),scale.unit = TRUE, graph = FALSE)
  fviz_pca_ind(mirnaData, habillage = as.factor(metaData$`tissue of origin a`), palette = "npg") + labs(title ="NC1-60 miRna data PCA plot", x = "PC1", y = "PC2")
  #mRNA data
  boxplot(agilentmRNA[,2:ncol(agilentmRNA)], xlab = "Samples", ylab= "mrna reads") 
  boxplot(agilentmiRNA[], xlab = "Samples", ylab= "miRNA reads") 
  vioplot(unlist(egfrData), xlab = "EGFR", ylab = "agilent reads")
#Data Analysis ----
  #making matrix to track which samples are what egfr status
  trackingMatrix <- data.frame(matrix(ncol = 3, nrow = length(colnames(agilentmiRNA))), row.names = (colnames(agilentmiRNA)))
  colnames(trackingMatrix) <- c("EGFR", "EGFR Pathway", "Amplification")
  #EGFR high vs low expression at median
  thr <- median(unlist(egfrData))
  pos <- as.factor(ifelse(egfrData > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  trackingMatrix$EGFR <- design_matrix[,2]
  fitmiRNA_1 <- eBayes(lmFit(agilentmiRNA, design = design_matrix))
  agilentEGFR <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  #EGFR pathway high vs low expression at median using ssgsea pathway
  allGeneSets <- msigdbr(species = "Homo sapiens")
  geneSet1 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGF_PATHWAY",]
  geneSet1 <- unique(geneSet1$gene_symbol)
  ssAnalysis <- gsva(as.matrix(agilentmRNAdupRemoved), list(geneSet1), method = "ssgsea", mx.diff = FALSE )
  thr <- median(as.vector(ssAnalysis))
  pos <- as.factor(ifelse(ssAnalysis > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  trackingMatrix$`EGFR Pathway` <- design_matrix[,2]
  fitmiRNA_1 <- eBayes(lmFit(agilentmiRNA, design = design_matrix))
  agilentPathway <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  #Amplification vs no Amplification
  pos <- as.factor(ifelse(cnvData > 2.25, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  trackingMatrix$Amplification<- design_matrix[,2]
  fitmiRNA_1 <- eBayes(lmFit(agilentmiRNA, design = design_matrix))
  agilentAmplification <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  #tracking matrix
  trackingMatrix <- cbind(trackingMatrix, rowSums(trackingMatrix))
  trackingMatrix <- cbind(trackingMatrix, c(metaData[,2:12]))
#RankProd ----
  common_ids <-Reduce('intersect', list(rownames(agilentEGFR), rownames(agilentPathway), rownames(agilentAmplification)))
  rankProd <- RP(cbind(agilentEGFR$logFC[match(common_ids, table=rownames(agilentEGFR))], agilentPathway$logFC[match(common_ids, table=rownames(agilentPathway))], agilentAmplification$logFC[match(common_ids, table=rownames(agilentAmplification))]), rep(1,3))
  common_ids <- t(as.data.frame(as.list(common_ids)))
  rankNCI60 <- topGene(rankProd,method="pfp",logged = TRUE, gene.names = common_ids, cutoff = 0.25)
