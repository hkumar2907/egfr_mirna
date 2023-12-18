#installing libraries ----
  library(dplyr)
  library(limma)
  library("readxl")
  library("msigdbr")
  library("GSVA")
  library("RankProd")
#loading data ----
  miRNAData <- read_excel('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/Normalized counts_miRNA.xlsx')
  mRNAData <-  read_excel('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/Normalized counts_mRNA.xlsx')
  supplementData <- read_excel('noac220_suppl_supplementary_table_s3.xlsx', sheet = 2)
  supplementData <- subset(supplementData, Hugo_Symbol == 'EGFR')
  initial_patients_miRNA <- select(miRNAData, matches("_Initial_"))
  initial_patients_miRNA <- as.data.frame(log2(initial_patients_miRNA))
  rownames(initial_patients_miRNA) <- miRNAData$...1
  initial_patients_mRNA <- as.data.frame(select(mRNAData, matches("_Initial_")))
  rownames(initial_patients_mRNA) <- mRNAData$...1
  colnames(initial_patients_miRNA) <- sub("_Initial_.*", "", colnames(initial_patients_miRNA))
  colnames(initial_patients_mRNA) <- sub("_Initial_.*", "", colnames(initial_patients_mRNA))
  supplementData$Patient_ID <- sub("GBM0", "", supplementData$Patient_ID)
  #converting miRNA base version to version 22
    list <- miRNAVersionConvert(rownames(initial_patients_miRNA),targetVersion="v22",exact=TRUE,verbose=TRUE) #converting to new version
    naIndex <- which(is.na(list$TargetName))
    initial_patients_miRNA <- initial_patients_miRNA[-naIndex,]
    rownames(initial_patients_miRNA) <- list$TargetName[-naIndex]
#visualization ----
  boxplot(miRNAData[,2:ncol(miRNAData)])
  boxplot(mRNAData[,2:ncol(mRNAData)])
#creating the four lists ----
  #list1
  thr <- median(as.numeric(initial_patients_mRNA[which(rownames(initial_patients_mRNA) == 'EGFR'),]))
  pos <- as.factor(ifelse(initial_patients_mRNA[which(rownames(initial_patients_mRNA) == 'EGFR'),] > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  fitmiRNA_1c <- lmFit(initial_patients_miRNA, design = design_matrix)
  fitmiRNA_1c <- eBayes(fitmiRNA_1c)
  tablefitmiRNA_1c <- topTable(fitmiRNA_1c, coef = 2, , number=Inf)
  #list2
  lowQuartile <- quantile(as.numeric(initial_patients_mRNA[which(rownames(initial_patients_mRNA) == 'EGFR'),]), probs = 0.25)
  highQuartile <- quantile(as.numeric(initial_patients_mRNA[which(rownames(initial_patients_mRNA) == 'EGFR'),]), probs = 0.75)
  pos <- as.factor(ifelse(initial_patients_mRNA[rownames(initial_patients_mRNA) == 'EGFR', ] > highQuartile, 'pos', ifelse(initial_patients_mRNA[rownames(initial_patients_mRNA) == 'EGFR', ] < lowQuartile, 'neg', NA)))
  design_matrix <- model.matrix(~(pos[!is.na(pos)]))
  fitmiRNA_2c <- lmFit(initial_patients_miRNA[,!is.na(pos)], design = design_matrix)
  fitmiRNA_2c <- eBayes(fitmiRNA_2c)
  tablefitmiRNA_2c <- topTable(fitmiRNA_2c, coef = 2, , number=Inf)
  #list3
  #allGeneSets <- msigdbr(species = "Homo sapiens")
  geneSet1 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGF_PATHWAY",]
  geneSet1 <- unique(geneSet1$gene_symbol)
  ssAnalysis <- gsva(as.matrix(initial_patients_mRNA), list(geneSet1), method = "ssgsea", mx.diff = FALSE )
  thr <- median(as.vector(ssAnalysis))
  pos <- as.factor(ifelse(ssAnalysis > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  fitmiRNA_3c <- lmFit(initial_patients_miRNA, design = design_matrix)
  fitmiRNA_3c <- eBayes(fitmiRNA_3c)
  tablefitmiRNA_3c <- topTable(fitmiRNA_3c, coef = 2, , number=Inf)
  #list4
  lowQuartile <- quantile(ssAnalysis, probs = 0.25)
  highQuartile <- quantile(ssAnalysis, probs = 0.75)
  pos <- as.factor(ifelse(ssAnalysis > highQuartile, 'pos', ifelse(ssAnalysis < lowQuartile, 'neg', NA)))
  design_matrix <- model.matrix(~(pos[!is.na(pos)]))
  fitmiRNA_4c <- lmFit(initial_patients_miRNA[,!is.na(pos)], design = design_matrix)
  fitmiRNA_4c <- eBayes(fitmiRNA_4c)
  tablefitmiRNA_4c <- topTable(fitmiRNA_4c, coef = 2, , number=Inf)
#list comparing mutations with no mutations ----
  pos <- as.factor(ifelse(colnames(initial_patients_miRNA) %in% supplementData$Patient_ID, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  fitmiRNA_5b <- lmFit(initial_patients_miRNA, design = design_matrix)
  fitmiRNA_5b <- eBayes(fitmiRNA_5b)
  tablefitmiRNA_5b <- topTable(fitmiRNA_5b, coef = 2, , number=Inf)
#RankProd ----
  common_ids <-Reduce('intersect', list(rownames(tablefitmiRNA_1c), rownames(tablefitmiRNA_2c), rownames(tablefitmiRNA_3c), rownames(tablefitmiRNA_4c)))
  rankProd <- RP(cbind(tablefitmiRNA_1c$logFC[match(common_ids, table=rownames(tablefitmiRNA_1c))], tablefitmiRNA_2c$logFC[match(common_ids, table=rownames(tablefitmiRNA_2c))], tablefitmiRNA_3c$logFC[match(common_ids, table=rownames(tablefitmiRNA_3c))], tablefitmiRNA_4c$logFC[match(common_ids, table=rownames(tablefitmiRNA_4c))]), rep(1,4))
  common_ids <- t(as.data.frame(as.list(common_ids)))
  rankProdstudy <- topGene(rankProd,method="pval",logged = TRUE, gene.names = common_ids, cutoff = 0.05)
  View(rankProdstudy[["Table1"]])
  View(rankProdstudy[["Table2"]])
  