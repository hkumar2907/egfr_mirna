#Installing Libraries ----
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
library(biclust)

#Mutation information: ----
# mutation_type <- (maf_matrix[,9]=='Missense_Mutation')*1 +
#   (maf_matrix[,9]=='Nonsense_Mutation')*2 +
#   (maf_matrix[,9]=='Frame_Shift_Del')*3 +
#   (maf_matrix[,9]=='Frame_Shift_Ins')*4 +
#   (maf_matrix[,9]=='Splice_Site')*5 +
#   (maf_matrix[,9]=='Translation_Start_Site')*6 +
#   (maf_matrix[,9]=='Nonstop_Mutation')*7 +
#   (maf_matrix[,9]=='3\'UTR')*8 +
#   (maf_matrix[,9]=='5\'UTR')*9 +
#   (maf_matrix[,9]=='3\'Flank')*10 +
#   (maf_matrix[,9]=='5\'Flank')*11 +
#   (maf_matrix[,9]=='RNA')*12 +
#   (maf_matrix[,9]=='Targeted_Region')*13 +
#   (maf_matrix[,9]=='In_Frame_Del')*14 +
#   (maf_matrix[,9]=='In_Frame_Ins')*15 +
#   (maf_matrix[,9]=='IGR')*16 +
#   (maf_matrix[,9]=='Intron')*17 +
#   (maf_matrix[,9]=='Silent')*18
# this is how i did it
# so missense = 1  nonsense = 2, etc


#LOADING all samples ----
#Creating cancer list
  cancer_types_list <- list();
  cancer_types_list[[1]] <- c('BRCA','UCEC','HNSC')
  cancer_types_list[[2]] <- c('KIRC','LUAD','THCA')
  cancer_types_list[[3]] <- c('PRAD','LUSC','OV')
  cancer_types_list[[4]] <- c('STAD','BLCA','COAD')
  cancer_types_list[[5]] <- c('LIHC','CESC','KIRP')
  all_cancer_types <- melt(cancer_types_list)$value
#Loading miRNA datasets
  all_miRNA_datasets <- list();
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    fname_mrna <- paste0('miRNA\ Datasets/', cancer_type,'.miRseq_mature_RPM_log2.txt')
    all_miRNA_datasets[[cancer_type]] <- read.table(fname_mrna, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
    all_miRNA_datasets[[cancer_type]] <- all_miRNA_datasets[[cancer_type]][!apply(is.na(all_miRNA_datasets[[cancer_type]][-1, ]), 1, all), ]
    rownames(all_miRNA_datasets[[cancer_type]]) <- all_miRNA_datasets[[cancer_type]][,1] #getting the miRNA names
    colnames(all_miRNA_datasets[[cancer_type]]) <- substr(colnames(all_miRNA_datasets[[cancer_type]]), 1, 12) #making same format as mRNA
    all_miRNA_datasets[[cancer_type]] <- all_miRNA_datasets[[cancer_type]][, -1] #removing the miRNA column
    #removing rows with all NA
  }
#Loading mRNA datasets
  all_mRNA_datasets <- list();
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    if(cancer_type!='BRCA'){
      fname_mrna <- paste0('Reprocessed\ mRNA\ Datasets/',cancer_type,'_cleaned_mRNA.txt')
    }else{
      fname_mrna <- paste0('Reprocessed\ mRNA\ Datasets/',cancer_type,'_cleaned_mRNA_ductal.txt')
    }
    all_mRNA_datasets[[cancer_type]] <- read.table(fname_mrna, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
    #colnames(all_mRNA_datasets[[cancer_type]]) <- gsub('[.]','-',colnames(all_mRNA_datasets[[cancer_type]]))
    # want log2 data
    all_mRNA_datasets[[cancer_type]] <- log2(all_mRNA_datasets[[cancer_type]]+1)
    all_mRNA_datasets[[cancer_type]][!is.finite(as.matrix(all_mRNA_datasets[[cancer_type]]))] <- NA
  }
#Matching miR and mRNA samples
  matched_miRNA_datasets <- list();
  matched_mRNA_datasets <- list();
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    commons <- intersect(colnames(all_mRNA_datasets[[cancer_type]]), colnames(all_miRNA_datasets[[cancer_type]]))
    matched_miRNA_datasets[[cancer_type]] <- all_miRNA_datasets[[cancer_type]][commons]
    matched_mRNA_datasets[[cancer_type]] <- all_mRNA_datasets[[cancer_type]][,c(which(colnames(all_mRNA_datasets[[cancer_type]]) %in% colnames(matched_miRNA_datasets[[cancer_type]])))]
  }
# #Importing non-mutated samples
#   mutation_datasets <- list();
#   for (cancer_type in all_cancer_types){
#     print(cancer_type)
#     mut_data <- paste0('TCGAMutationData/',cancer_type,'.mutations.txt')
#     mut_data <- read.table(mut_data, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
#     #only selecting egfr
#     mutation_datasets[[cancer_type]] <- mut_data[rownames(mut_data) == 'EGFR',]
#     #removing 0 and 18 (no mutation and silent mutation)
#     mutation_datasets[[cancer_type]] <- mutation_datasets[[cancer_type]][which((mutation_datasets[[cancer_type]] != 0) & (mutation_datasets[[cancer_type]] != 18))]
#     mutation_datasets[[cancer_type]] <- mutation_datasets[[cancer_type]]intersect(names(mutation_datasets[[cancer_type]]), names(matched_miRNA_datasets[[cancer_type]]))
#     #mutation_datasets[[cancer_type]] <- mutation_datasets[[cancer_type]][which((mutation_datasets[[cancer_type]] != 0))]
#   }
#Importing CNV data
  amplificationData <- list()
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    cnv_data <- paste0('CNVdata/',cancer_type,'.cleaned_cnv.txt')
    cnv_data <- read.table(cnv_data, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
    #only selecting egfr
    amplificationData[[cancer_type]] <- cnv_data[rownames(cnv_data) == '1956',]
    #matching samples to mirna datasets and mrna datasets
    amplificationData[[cancer_type]] <- amplificationData[[cancer_type]][intersect(names(amplificationData[[cancer_type]]), names(matched_miRNA_datasets[[cancer_type]]))]
    #taking samples that are above 2.25 (so amplified)
    amplificationData[[cancer_type]] <- amplificationData[[cancer_type]][which(amplificationData[[cancer_type]] > 2.25)]
  } 
  # nonmutated_samples <- read_csv("nonmutated_samples.csv")
  # nonmutated_samples <- nonmutated_samples[,-1]
  # nonmutated_mRNA_datasets <- list();
  # for (cancer_type in all_cancer_types){
  #   if ((cancer_type != "THCA") && (cancer_type != "COAD") && (cancer_type != "KIRP")) {
  #     list <- unlist(nonmutated_samples[cancer_type])
  #     #nonmutated_mRNA_datasets[[cancer_type]] <- select(all_miRNA_datasets[[cancer_type]], list)
  #     nonmutated_mRNA_datasets[[cancer_type]] <- intersect(names(all_miRNA_datasets[[cancer_type]]), list)
  #   }
  # }
  
#DOING THE ANALYSIS ----
#EGFR high vs low expression at median
  #getting the EGFR data (gene id: 1956) from the mRNA data and then sorting
  tablefit_egfrexpr <- list();
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    thr <- median(as.numeric(matched_mRNA_datasets[[cancer_type]][which(rownames(matched_mRNA_datasets[[cancer_type]]) == '1956'),]))
    pos <- as.factor(ifelse(matched_mRNA_datasets[[cancer_type]][which(rownames(matched_mRNA_datasets[[cancer_type]]) == '1956'),] > thr, 'pos', 'neg'))
    design_matrix <- model.matrix(~pos)
    fitmiRNA_1 <- lmFit(matched_miRNA_datasets[[cancer_type]][,match(colnames(matched_mRNA_datasets[[cancer_type]]), table=colnames(matched_miRNA_datasets[[cancer_type]]))], design = design_matrix)
    fitmiRNA_1 <- eBayes(fitmiRNA_1)
    tablefit_egfrexpr[[cancer_type]] <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  }
#EGFR pathway high vs low expression at median using ssgsea pathway
  allGeneSets <- msigdbr(species = "Homo sapiens")
  #geneSet1 <- allGeneSets[allGeneSets$gs_name == "REACTOME_SIGNALING_BY_EGFR",]
  geneSet1 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGF_PATHWAY",]
  geneSet1 <- unique(geneSet1$entrez_gene)
  tablefit_ssgseaexpr <- list();
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    ssAnalysis <- gsva(as.matrix( matched_mRNA_datasets[[cancer_type]]), list(geneSet1), method = "ssgsea", mx.diff = FALSE )
    thr <- median(as.vector(ssAnalysis))
    pos <- as.factor(ifelse(ssAnalysis > thr, 'pos', 'neg'))
    design_matrix <- model.matrix(~pos)
    fitmiRNA_1 <- lmFit(matched_miRNA_datasets[[cancer_type]][,match(colnames(matched_mRNA_datasets[[cancer_type]]), table=colnames(matched_miRNA_datasets[[cancer_type]]))], design = design_matrix)
    fitmiRNA_1 <- eBayes(fitmiRNA_1)
    tablefit_ssgseaexpr[[cancer_type]] <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  }
# #EGFR mutation vs non-mutation
#   tablefit_mutationexpr <- list();
#   for (cancer_type in all_cancer_types) {
#     if ((cancer_type != "THCA") && (cancer_type != "COAD") && (cancer_type != "KIRP")) {
#       print(cancer_type)
#       mutatedList <- intersect(colnames(all_miRNA_datasets[[cancer_type]]), pull(nonmutated_samples[cancer_type]))
#       pos <- as.factor(ifelse(colnames(matched_mRNA_datasets[[cancer_type]]) %in% mutatedList, 'pos', 'neg'))
#       design_matrix <- model.matrix(~pos)
#       fitmiRNA_1 <- lmFit(matched_miRNA_datasets[[cancer_type]][,match(colnames(matched_mRNA_datasets[[cancer_type]]), table=colnames(matched_miRNA_datasets[[cancer_type]]))], design = design_matrix)
#       fitmiRNA_1 <- eBayes(fitmiRNA_1)
#       tablefit_mutationexpr[[cancer_type]] <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
#     }
#   }
#Amplification vs no Amplification
  tablefit_amplification <- list();
  for (cancer_type in all_cancer_types) {
    print(cancer_type)
    pos <- as.factor(ifelse(colnames(matched_miRNA_datasets[[cancer_type]]) %in% names(amplificationData[[cancer_type]]), 'pos', 'neg'))
    design_matrix <- model.matrix(~pos)
    fitmiRNA_1 <- lmFit(matched_miRNA_datasets[[cancer_type]][,match(colnames(matched_mRNA_datasets[[cancer_type]]), table=colnames(matched_miRNA_datasets[[cancer_type]]))], design = design_matrix)
    fitmiRNA_1 <- eBayes(fitmiRNA_1)
    tablefit_amplification[[cancer_type]] <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  }
#Amplification vs Deletion
#RANKPROD ----
  #previous
  # rankProdssgsea <- rankprod(tablefit_ssgseaexpr)
  # rankProdexpr <- rankprod(tablefit_egfrexpr)
  # rankProdMutation <- rankprod(tablefit_mutationexpr)
  # #function
  # rankprod <- function(tablefit) {
  #   common_ids <- rownames(tablefit[[i]])
  #   c <- tablefit[[1]]$t[match(common_ids, table=rownames(tablefit[[1]]))]
  #   for (i in 2:length(tablefit)) {
  #     common_ids <- intersect(rownames(tablefit[[i]]), common_ids) #PROBLEM
  #     c <- cbind(c, tablefit[[i]]$t[match(common_ids, table=rownames(tablefit[[i]]))])
  #   }
  #   rownames(c) <- common_ids
  #   rankProd <- RP(c, rep(1,length(tablefit)))
  #   common_ids <- t(as.data.frame(as.list(common_ids)))
  #   rankGenes <- topGene(rankProd,method="pval",logged = FALSE, gene.names = common_ids, cutoff = 0.05)
  #   return(rankGenes)
  # }
  # View(rankProdexpr[["Table1"]])
  # View(rankProdexpr[["Table2"]])
  # View(rankProdssgsea[["Table1"]])
  # View(rankProdssgsea[["Table2"]])
  # View(rankProdMutation[["Table1"]])
  # View(rankProdMutation[["Table2"]])
  #current
  rankprod <- function(ssgsea, expr, amplification) {
    common_ids <-Reduce('intersect', list(rownames(ssgsea), rownames(expr), rownames(amplification)))
    rankProd <- RP(cbind(ssgsea$t[match(common_ids, table=rownames(ssgsea))], expr$t[match(common_ids, table=rownames(expr))], amplification$t[match(common_ids, table=rownames(amplification))]), rep(1,3))
    common_ids <- t(as.data.frame(as.list(common_ids)))
    rankGenes <- topGene(rankProd,method="pfp",logged = TRUE, gene.names = common_ids, cutoff = 0.1)
    return(rankGenes)
  }
  rankProdAllCancers <- list();
  for (cancer_type in all_cancer_types) {
    print(cancer_type)
    rankProdAllCancers[[cancer_type]] <- rankprod(tablefit_ssgseaexpr[[cancer_type]], tablefit_egfrexpr[[cancer_type]], tablefit_amplification[[cancer_type]])
  }
#Finding commonalities between cancers ----
  #intersectedTable1 <- rownames(rankProdAllCancers[["BRCA"]][["Table1"]])
  #intersectedTable2 <- rownames(rankProdAllCancers[["BRCA"]][["Table2"]])
  combinedListTable1 <- list()
  combinedListTable2 <- list()
  num <- 50 #variable for selecting how many top mirs to choose
  for (cancer_type in all_cancer_types) {
    # intersectedTable1 <- intersect(intersectedTable1, rownames(rankProdAllCancers[[cancer_type]][["Table1"]]))
    # intersectedTable2 <- intersect(intersectedTable2, rownames(rankProdAllCancers[[cancer_type]][["Table2"]]))
    combinedListTable1 <- c(combinedListTable1, rownames(rankProdAllCancers[[cancer_type]][["Table1"]])[1:num])
    combinedListTable2 <- c(combinedListTable2, rownames(rankProdAllCancers[[cancer_type]][["Table2"]])[1:num])
  }
  minmirs <- 3
  counts <- table(t(as.data.frame(combinedListTable1)))
  counts2 <- table(t(as.data.frame(combinedListTable2)))
  counts2 <- counts2[counts2>=minmirs]
  counts <- counts[counts>=minmirs]
  #significantgenes <- counts[counts>5]
  #significantgenestable1 <- data.frame(matrix(ncol = 15, nrow = length(significantgenes)), row.names = names(significantgenes))
  #colnames(significantgenestable1) <- all_cancer_types
  allGenesTable1 <- data.frame(matrix(ncol = 15, nrow = length(counts)), row.names = names(counts))
  colnames(allGenesTable1) <- all_cancer_types
  allGenesTable2 <- data.frame(matrix(ncol = 15, nrow = length(counts2)), row.names = names(counts2))
  colnames(allGenesTable2) <- all_cancer_types
  FCGenesTable1 <- data.frame(matrix(ncol = 15, nrow = length(counts)), row.names = names(counts))
  colnames(FCGenesTable1) <- all_cancer_types
  FCGenesTable2 <- data.frame(matrix(ncol = 15, nrow = length(counts2)), row.names = names(counts2))
  colnames(FCGenesTable2) <- all_cancer_types
  # significantgenes2 <- counts2[counts2>5]
  # significantgenestable2 <- data.frame(matrix(ncol = 15, nrow = length(significantgenes2)), row.names = names(significantgenes2))
  # colnames(significantgenestable2) <- all_cancer_types
  for (cancer_type in all_cancer_types) {
    # #for significant mir table
    # matching <-  match(names(significantgenes), rownames(rankProdAllCancers[[cancer_type]][["Table1"]]))
    # significantgenestable1[cancer_type] <- rankProdAllCancers[[cancer_type]][["Table1"]][matching,3]
    # matching <-  match(names(significantgenes2), rownames(rankProdAllCancers[[cancer_type]][["Table2"]]))
    # significantgenestable2[cancer_type] <- rankProdAllCancers[[cancer_type]][["Table2"]][matching,3]
    #for all mirtable
    print(cancer_type)
    matching <-  match(names(counts), rownames(rankProdAllCancers[[cancer_type]][["Table1"]])[1:num])
    matching[is.na(matching)] <- 0 #replacing all na values with zero
    matching[matching > 0] <- 1 #replacing all other values with one
    allGenesTable1[cancer_type] <- matching #setting that as the value
    matching2 <-  match(names(counts2), rownames(rankProdAllCancers[[cancer_type]][["Table2"]])[1:num])
    matching2[is.na(matching2)] <- 0 #replacing all na values with zero
    matching2[matching2 > 0] <- 1 #replacing all other values with one
    allGenesTable2[cancer_type] <- matching2
    #TRYING with rankprod values
    matching <-  match(names(counts), rownames(rankProdAllCancers[[cancer_type]][["Table1"]])[1:num])
    matching[is.na(matching)] <- 0 #replacing all na values with zero
    matching[matching>0] <- 1/(rankProdAllCancers[[cancer_type]][["Table1"]][ matching[matching>0],2]) #taking 1/rankprod vals so lower values are highest
    FCGenesTable1[cancer_type] <- matching
    matching <-  match(names(counts2), rownames(rankProdAllCancers[[cancer_type]][["Table2"]])[1:num])
    matching[is.na(matching)] <- 0 #replacing all na values with zero
    matching[matching>0] <- 1/(rankProdAllCancers[[cancer_type]][["Table2"]][ matching[matching>0],2]) #taking 1/rankprod vals so lower values are highest
    FCGenesTable2[cancer_type] <- matching
  }
  # table1PCA <- PCA(t(significantgenestable1), graph = FALSE)
  # fviz_pca_ind(table1PCA) + labs(title ="Table 1 PCA")
  # table2PCA <- PCA(t(significantgenestable2), graph = FALSE)
  # fviz_pca_ind(table2PCA) + labs(title ="Table 2 PCA")
#ConsensusClusterPlus ----
  title=tempdir()
  results = ConsensusClusterPlus(as.matrix(FCGenesTable1),maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,clusterAlg="pam",distance="pearson",seed=1262118388.71279,plot="png")
#Creating cluster heatmap ----
  #creating combined tables
  commonRows <- intersect(rownames(FCGenesTable1), rownames(FCGenesTable2))
  commonVal <- FCGenesTable1[commonRows,] - FCGenesTable2[commonRows,] #just adding the vals?
  combinedFCGeneTable <- rbind(FCGenesTable1[!(row.names(FCGenesTable1) %in% commonRows),], -1*FCGenesTable2[!(row.names(FCGenesTable2) %in% commonRows),], commonVal)
  commonVal <- allGenesTable1[commonRows,] - allGenesTable2[commonRows,] #just adding the vals?
  combinedGeneTable <- rbind(allGenesTable1[!(row.names(allGenesTable1) %in% commonRows),], -1*allGenesTable2[!(row.names(allGenesTable2) %in% commonRows),], commonVal)
  #for binary heatmap
  col_fun = colorRamp2(c(0, 0.5, 1), c("white", "grey", "black"))
  Heatmap(allGenesTable1,col= col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 6), column_title = "Table 1 (mir in top 50 or not)")
  Heatmap(allGenesTable2,col= col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 6), column_title = "Table 2 (mir in top 50 or not)")
  col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "blue")) #for combined
  Heatmap(combinedGeneTable,col= col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 6), column_title = "Combined Tables (mir in top 50 or not)")
  #for 1/RP val non-binary heatmap
  col_fun = colorRamp2(c(0, 0.1, 0.2), c("white", "gray35", "black"))
  Heatmap(FCGenesTable1,col= col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 6), column_title = "Table 1 w/ RP(mir in top 50 or not)")
  Heatmap(FCGenesTable2,col= col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 6), column_title = "Table 2 w/ RP(mir in top 50 or not)")
  col_fun = colorRamp2(c(-0.15, 0, 0.15), c("red4", "white", "royalblue4")) #for combined
  Heatmap(combinedFCGeneTable,col= col_fun, border = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 6), column_title = "Combined Tables w/ RP(mir in top 50 or not)")
  #for PCA plot
  table1PCA <- PCA(t(allGenesTable1), graph = FALSE)
  fviz_pca_ind(table1PCA) + labs(title ="Table 1 PCA")
  table2PCA <- PCA(t(allGenesTable2), graph = FALSE)
  fviz_pca_ind(table2PCA) + labs(title ="Table 2 PCA")
  tableComPCA <- PCA(t(combinedGeneTable), graph = FALSE)
  fviz_pca_ind(tableComPCA) + labs(title ="Combined PCA")
#Biclustering ----
  biclustFC <- biclust(as.matrix(combinedFCGeneTable), method = BCCC())

  
  
  biclustFC <- biclust(data.matrix(combinedFCGeneTable), method=BCQuest(), ns=10, nd=10, sd=5, alpha=0.05, number=100)
  