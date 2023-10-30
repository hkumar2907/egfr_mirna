library(RankProd)
library(reshape2)
library(readr)
library(penalized)
library(dplyr)
library("GSVA")
library("msigdbr")
library("RankProd")

#mutation information: ----
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
#creating cancer list
cancer_types_list <- list();
cancer_types_list[[1]] <- c('BRCA','UCEC','HNSC')
cancer_types_list[[2]] <- c('KIRC','LUAD','THCA')
cancer_types_list[[3]] <- c('PRAD','LUSC','OV')
cancer_types_list[[4]] <- c('STAD','BLCA','COAD')
cancer_types_list[[5]] <- c('LIHC','CESC','KIRP')
all_cancer_types <- melt(cancer_types_list)$value
#loading miRNA datasets
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
#loading mRNA datasets
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
#matching miR and mRNA samples
matched_miRNA_datasets <- list();
matched_mRNA_datasets <- list();
for (cancer_type in all_cancer_types){
  print(cancer_type)
  commons <- intersect(colnames(all_mRNA_datasets[[cancer_type]]), colnames(all_miRNA_datasets[[cancer_type]]))
  matched_miRNA_datasets[[cancer_type]] <- all_miRNA_datasets[[cancer_type]][commons]
  matched_mRNA_datasets[[cancer_type]] <- all_mRNA_datasets[[cancer_type]][,c(which(colnames(all_mRNA_datasets[[cancer_type]]) %in% colnames(matched_miRNA_datasets[[cancer_type]])))]
}
#importing non-mutated samples
nonmutated_samples <- read_csv("nonmutated_samples.csv")
nonmutated_samples <- nonmutated_samples[,-1]
nonmutated_mRNA_datasets <- list();
for (cancer_type in all_cancer_types){
  if ((cancer_type != "THCA") && (cancer_type != "COAD") && (cancer_type != "KIRP")) {
    list <- unlist(nonmutated_samples[cancer_type])
    #nonmutated_mRNA_datasets[[cancer_type]] <- select(all_miRNA_datasets[[cancer_type]], list)
    nonmutated_mRNA_datasets[[cancer_type]] <- intersect(names(all_miRNA_datasets[[cancer_type]]), list)
  }
}
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
  geneSet1 <- allGeneSets[allGeneSets$gs_name == "HALLMARK_PI3K_AKT_MTOR_SIGNALING",]
  geneSet1 <- unique(geneSet1$human_entrez_gene)
  tablefit_ssgseaexpr <- list();
  for (cancer_type in all_cancer_types){
    print(cancer_type)
    thr <- median(as.numeric(matched_mRNA_datasets[[cancer_type]][which(rownames(matched_mRNA_datasets[[cancer_type]]) == '1956'),]))
    pos <- as.factor(ifelse(matched_mRNA_datasets[[cancer_type]][which(rownames(matched_mRNA_datasets[[cancer_type]]) == '1956'),] > thr, 'pos', 'neg'))
    design_matrix <- model.matrix(~pos)
    fitmiRNA_1 <- lmFit(matched_miRNA_datasets[[cancer_type]][,match(colnames(matched_mRNA_datasets[[cancer_type]]), table=colnames(matched_miRNA_datasets[[cancer_type]]))], design = design_matrix)
    fitmiRNA_1 <- eBayes(fitmiRNA_1)
    tablefit_ssgseaexpr[[cancer_type]] <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
  }
#EGFR mutation vs non-mutation
  tablefit_mutationexpr <- list();
  for (cancer_type in all_cancer_types) {
    if ((cancer_type != "THCA") && (cancer_type != "COAD") && (cancer_type != "KIRP")) {
      print(cancer_type)
      mutatedList <- intersect(colnames(all_miRNA_datasets[[cancer_type]]), pull(nonmutated_samples[cancer_type]))
      pos <- as.factor(ifelse(colnames(matched_mRNA_datasets[[cancer_type]]) %in% mutatedList, 'pos', 'neg'))
      design_matrix <- model.matrix(~pos)
      fitmiRNA_1 <- lmFit(matched_miRNA_datasets[[cancer_type]][,match(colnames(matched_mRNA_datasets[[cancer_type]]), table=colnames(matched_miRNA_datasets[[cancer_type]]))], design = design_matrix)
      fitmiRNA_1 <- eBayes(fitmiRNA_1)
      tablefit_mutationexpr[[cancer_type]] <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
    }
  }


#RANKPROD ----
  rankProdssgsea <- rankprod(tablefit_ssgseaexpr)
  rankProdexpr <- rankprod(tablefit_egfrexpr)
  rankProdMutation <- rankprod(tablefit_mutationexpr)
  #function
  rankprod <- function(tablefit) {
    common_ids <- rownames(tablefit[[i]])
    c <- tablefit[[1]]$t[match(common_ids, table=rownames(tablefit[[1]]))]
    for (i in 2:length(tablefit)) {
      common_ids <- intersect(rownames(tablefit[[i]]), common_ids) #PROBLEM
      c <- cbind(c, tablefit[[i]]$t[match(common_ids, table=rownames(tablefit[[i]]))])
    }
    rownames(c) <- common_ids
    rankProd <- RP(c, rep(1,length(tablefit)))
    common_ids <- t(as.data.frame(as.list(common_ids)))
    rankGenes <- topGene(rankProd,method="pval",logged = FALSE, gene.names = common_ids, cutoff = 0.05)
    return(rankGenes)
  }
  View(rankProdexpr[["Table1"]])
  View(rankProdexpr[["Table2"]])
  View(rankProdssgsea[["Table1"]])
  View(rankProdssgsea[["Table2"]])
  View(rankProdMutation[["Table1"]])
  View(rankProdMutation[["Table2"]])
  