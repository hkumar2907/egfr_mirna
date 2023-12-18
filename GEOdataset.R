#Load appropriate libraries ----
library(GEOquery)
library(GSVA)
library(msigdbr)
library(SummarizedExperiment)
library(pheatmap)
library(tidyverse)
library(maftools)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(biomaRt)
library(stringr)
library(hugene21sttranscriptcluster.db)
library(miRBaseVersions.db)
library(miRNAtap) 
library(affy)
library(pheatmap)
library(limma)
library(hgu133plus2.db)
library(clusterProfiler)
library(corto)
library(gsva)
library(pd.mirna.4.0) #DNA sequences, not useful here? 
library(RankProd)

#Loading geo dataset ----
#Step 1: Query, download, and prep GEO data sets 
  GEOquery <- getGEO(GEO = "GSE90604", GSEMatrix = TRUE)
  gene_exprs <- exprs(GEOquery[[1]])
  gene_exprs <- gene_exprs[,1:16]
  miRNA_exprs <- exprs(GEOquery[[2]])
  colnames(miRNA_exprs) = pData(GEOquery[[2]])$title #Rename sample (column) names 
  colnames(miRNA_exprs) <- gsub("Glioblastoma Multiforme", "GBM", colnames(miRNA_exprs))
  miRNA_exprs <- miRNA_exprs[,1:16] #isolate GBM samples, remove any healthy samples 
#miRNA labels
  miRNA_exprs_t = as.data.frame(miRNA_exprs)
  miRNA_exprs_t$ID = rownames(miRNA_exprs)
  miRNA_annotation_table <- read.table("GEO GBM miRNA annotation.txt", header = TRUE, sep = "\t", fill = TRUE)
  miRNA_annotation_table = miRNA_annotation_table[,c(1,4)] #only want the ID and names from annotation file
  miRNA_exprs_t <- merge(miRNA_annotation_table, miRNA_exprs_t, by = "ID")  #change rownames to miRNA names based on annotation file 
  colnames(miRNA_annotation_table)[2] = 'miRNA'
  rownames(miRNA_exprs_t) <- miRNA_exprs_t$Transcript.ID.Array.Design.
  miRNA_exprs_t$ID <- NULL
  miRNA_exprs_t$Transcript.ID.Array.Design. <- NULL
  #converting miRNA base version to version 22
    list <- miRNAVersionConvert(rownames(miRNA_exprs_t),targetVersion="v22",exact=TRUE,verbose=TRUE) #converting to new version
    naIndex <- which(is.na(list$TargetName))
    miRNA_exprs_t <- miRNA_exprs_t[-naIndex,]
    rownames(miRNA_exprs_t) <- list$TargetName[-naIndex]
#mRNA labels
  table <- read.table("GBM GEO annotation 2.txt", header = TRUE, sep = "\t", fill = TRUE)
  gene_exprs<- cbind(rownames(gene_exprs), gene_exprs)
  colnames(gene_exprs)[1] ="ProbeName"
  gene_exprs <- merge(gene_exprs, table, by = "ProbeName")
  gene_exprs <- gene_exprs[!duplicated(gene_exprs$GeneSymbols), ]
  rownames(gene_exprs) <- gene_exprs$GeneSymbols
  gene_exprs <- gene_exprs[, 1:(ncol(gene_exprs) - 5)]
  gene_exprs <- gene_exprs[, -1]
#Creating the lists ----
#list1
  thr <- median(as.numeric(gene_exprs[which(rownames(gene_exprs) == 'EGFR'),]))
  pos <- as.factor(ifelse(gene_exprs[which(rownames(gene_exprs) == 'EGFR'),] > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  fitmiRNA_1b <- lmFit(miRNA_exprs_t, design = design_matrix)
  fitmiRNA_1b <- eBayes(fitmiRNA_1b)
  tablefitmiRNA_1b <- topTable(fitmiRNA_1b, coef = 2, , number=Inf)
#list2
  lowQuartile <- quantile(as.numeric(gene_exprs[which(rownames(gene_exprs) == 'EGFR'),]), probs = 0.25)
  highQuartile <- quantile(as.numeric(gene_exprs[which(rownames(gene_exprs) == 'EGFR'),]), probs = 0.75)
  pos <- as.factor(ifelse(gene_exprs[rownames(gene_exprs) == 'EGFR', ] > highQuartile, 'pos', ifelse(gene_exprs[rownames(gene_exprs) == 'EGFR', ] < lowQuartile, 'neg', NA)))
  design_matrix <- model.matrix(~(pos[!is.na(pos)]))
  fitmiRNA_2b <- lmFit(miRNA_exprs_t[,!is.na(pos)], design = design_matrix)
  fitmiRNA_2b <- eBayes(fitmiRNA_2b)
  tablefitmiRNA_2b <- topTable(fitmiRNA_2b, coef = 2, , number=Inf)
#list3
  #allGeneSets <- msigdbr(species = "Homo sapiens")
  geneSet1 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGF_PATHWAY",]
  geneSet1 <- unique(geneSet1$gene_symbol)
  for (i in 1:16){
    gene_exprs[,i] <- as.numeric(gene_exprs[,i])
  }
  ssAnalysis <- gsva(as.matrix(gene_exprs), list(geneSet1), method = "ssgsea", mx.diff = FALSE )
  thr <- median(as.vector(ssAnalysis))
  pos <- as.factor(ifelse(ssAnalysis > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  fitmiRNA_3b <- lmFit(miRNA_exprs_t, design = design_matrix)
  fitmiRNA_3b <- eBayes(fitmiRNA_3b)
  tablefitmiRNA_3b <- topTable(fitmiRNA_3b, coef = 2, , number=Inf)
#list4
  lowQuartile <- quantile(ssAnalysis, probs = 0.25)
  highQuartile <- quantile(ssAnalysis, probs = 0.75)
  pos <- as.factor(ifelse(ssAnalysis > highQuartile, 'pos', ifelse(ssAnalysis < lowQuartile, 'neg', NA)))
  design_matrix <- model.matrix(~(pos[!is.na(pos)]))
  fitmiRNA_4b <- lmFit(miRNA_exprs_t[,!is.na(pos)], design = design_matrix)
  fitmiRNA_4b <- eBayes(fitmiRNA_4b)
  tablefitmiRNA_4b <- topTable(fitmiRNA_4b, coef = 2, , number=Inf)
#RankProd ----
  library("RankProd")
  common_ids <-Reduce('intersect', list(rownames(tablefitmiRNA_1b), rownames(tablefitmiRNA_2b), rownames(tablefitmiRNA_3b), rownames(tablefitmiRNA_4b)))
  rankProd <- RP(cbind(tablefitmiRNA_1b$logFC[match(common_ids, table=rownames(tablefitmiRNA_1b))], tablefitmiRNA_2b$logFC[match(common_ids, table=rownames(tablefitmiRNA_2b))], tablefitmiRNA_3b$logFC[match(common_ids, table=rownames(tablefitmiRNA_3b))], tablefitmiRNA_4b$logFC[match(common_ids, table=rownames(tablefitmiRNA_4b))]), rep(1,4))
  common_ids <- t(as.data.frame(as.list(common_ids)))
  rankProdGeo <- topGene(rankProd,method="pval",logged = TRUE, gene.names = common_ids, cutoff = 0.05)
  View(rankProdGeo[["Table1"]])
  View(rankProdGeo[["Table2"]])
