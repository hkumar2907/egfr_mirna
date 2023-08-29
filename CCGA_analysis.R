# Installing Packages ----
  install.packages("dplyr")                          # Install dplyr package
  library("dplyr") # Load dplyr package
  install.packages("corrr")
  library('corrr')
  install.packages("ggcorrplot")
  library(ggcorrplot)
  install.packages("FactoMineR")
  library("FactoMineR")
  install.packages("factoextra")
  library("factoextra")
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2")
  library( "DESeq2" )
  BiocManager::install("limma")
  BiocManager::install("GSVA", force = TRUE)
  library("GSVA")
  install.packages("jsonlite")
  library(jsonlite)
  BiocManager::install("fgsea")
  library("fgsea")
  install.packages("survival")  # Install the package
  library(survival)
  BiocManager::install("RankProd")
  install.packages("msigdbr")
  library("msigdbr")
  install.packages("crosstalkr")
  BiocManager::install("STRINGdb")
# Loading CGGA Data -----
  #Dataset 2 mRNAseq with 693 samples NGS method
  clinical693 <- read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/CGGA.mRNAseq_693_clinical.20200506.txt', header = TRUE, fill = TRUE)
  clinical693 <- clinical693[, -c(4, 15:18)]
  names(clinical693) <- c("CGGA_ID", "PRS_type", "Histology", "Grade", "Gender", "Age", "OS", "Censor", "Radio_Status", "Chemo_status", "IDH_mutation_status", "X1p19q_codeletion_status", "MGMTp_methylation_status")
  mRNA693 <- read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/CGGA.mRNAseq_693.RSEM-genes.20200506.txt', header = TRUE, fill = TRUE)
  #Dataset 6 for microRNA arrray with 198 samples
  clinicalmiRNA <- read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/CGGA.microRNA_array_198_clinical.20200506.txt', header = TRUE, fill = TRUE)
  clinicalmiRNA <- clinicalmiRNA[, -c(4, 15:18)]
  names(clinicalmiRNA) <- c("CGGA_ID", "PRS_type", "Histology", "Grade", "Gender", "Age", "OS", "Censor", "Radio_Status", "Chemo_status", "IDH_mutation_status", "X1p19q_codeletion_status", "MGMTp_methylation_status")
  arraymiRNA <- read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/CGGA.microRNA_array_198_gene_level.20200506.txt', header = TRUE, fill = TRUE)
  #Dataset 4 for mRNA array data with 301 samples
  clinicalmRNAarray <- read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/CGGA.mRNA_array_301_clinical.20200506.txt', header = TRUE, fill = TRUE)
  clinicalmRNAarray <- clinicalmRNAarray[, -c(4:5, 15:18)]
  names(clinicalmRNAarray) <- c("CGGA_ID", "PRS_type", "Histology", "Grade", "Gender", "Age", "OS", "Censor", "Radio_Status", "Chemo_status", "IDH_mutation_status", "X1p19q_codeletion_status", "MGMTp_methylation_status")
  arraymRNA <- read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/CGGA.mRNA_array_301_gene_level.20200506.txt', header = TRUE, fill = TRUE)
  #finding out which samples are common for dataset 4&6
  library("dplyr") # Load dplyr package
  commonClinical = semi_join(clinicalmiRNA, clinicalmRNAarray, by = "CGGA_ID")
  commonPatients = intersect(colnames(arraymiRNA), colnames(arraymRNA))
  miRNACommon = arraymiRNA[, c(1, which(colnames(arraymiRNA) %in% commonPatients))]
  miRNACommon = na.omit(miRNACommon)
  #grabbing only primary GBM patients
  primaryGBMID <- commonClinical$PRS_type == "Primary" & commonClinical$Histology == "GBM"
  commonClinical <- commonClinical[primaryGBMID,] 
  miRNACommon <- miRNACommon[,c(TRUE,primaryGBMID)]
  mRNACommon <- arraymRNA[, c(1, which(colnames(arraymRNA) %in% commonClinical[,1]))]
# Visualization ----
  #creating a box plot of all the means for each patient (patient on x axis)
  boxplot(mRNACommon[,2:ncol(mRNACommon)]) #all means centered at zero
  library("FactoMineR")
  arrayPCA <- PCA(mRNACommon[,2:ncol(mRNACommon)], graph = FALSE)
  #performing principal component analysis by transposing the array
  patienPCA <- PCA(t(mRNACommon[,2:ncol(mRNACommon)]),scale.unit = TRUE, graph = FALSE)
  #graphing the individual pca
  library("factoextra")
  fviz_pca_ind(patienPCA)
  #plots the variances (%) for each PC 
  fviz_eig(patienPCA, addlabels = TRUE)
  #manually calculating the standard for the mRNA to see if some points should be eliminated
  patientmRNAsd <- apply(mRNACommon[,2:ncol(mRNACommon)],2,sd)
  patientmRNAsd <- patientmRNAsd[order(patientmRNAsd, decreasing = TRUE)]
  #MIRNA DATA
  #creating a box plot of all the means for each patient (patient on x axis)
  boxplot(miRNACommon[,2:ncol(miRNACommon)]) #all means centered at zero
  arrayPCAmi <- PCA(miRNACommon[,2:ncol(miRNACommon)], graph = FALSE)
  #performing principal component analysis by transposing the array
  patienPCAmi <- PCA(t(miRNACommon[,2:ncol(miRNACommon)]),scale.unit = TRUE, graph = FALSE)
  #graphing the individual pca
  fviz_pca_ind(patienPCAmi)
  #plots the variances (%) for each PC 
  fviz_eig(patienPCAmi, addlabels = TRUE)
# Further preprocessing and MicroRNA preprocessing + Visualization ----
  #need to first remove genes that have very low variance across all samples
  genemRNAsd <- apply(t(mRNACommon[,2:ncol(mRNACommon)]),2,sd)
  genemRNAmean <- apply(t(mRNACommon[,2:ncol(mRNACommon)]),2,mean)
  #creating scatter plot of sd vs mean
  plot(genemRNAsd, genemRNAmean)
  #histogram of variance
  hist(genemRNAsd, breaks = 40)
  library("limma")
  normalizedmiRNACommon = normalizeBetweenArrays(miRNACommon[, 2:ncol(miRNACommon)], method = "quantile")
  rownames(normalizedmiRNACommon) <- miRNACommon$microRNA_ID
  boxplot(normalizedmiRNACommon)
  miRNAPCA <- PCA(t(normalizedmiRNACommon),scale.unit = TRUE, graph = FALSE)
  fviz_pca_ind(miRNAPCA)
  fviz_eig(miRNAPCA, addlabels = TRUE)
  row_ids <- mRNACommon$Gene_Name
  mRNACommon <- data.matrix(mRNACommon[,-1])
  rownames(mRNACommon) <- row_ids
#List 1 miRNA Differential Analysis ----
  #getting the EGFR data from the mRNA data and then sorting
  thr <- median(as.numeric(mRNACommon[which(rownames(mRNACommon) == 'EGFR'),]))
  pos <- as.factor(ifelse(mRNACommon[which(rownames(mRNACommon) == 'EGFR'),] > thr, 'pos', 'neg'))
  design_matrix <- model.matrix(~pos)
  fitmiRNA_1 <- lmFit(normalizedmiRNACommon[,match(colnames(mRNACommon), table=colnames(normalizedmiRNACommon))], design = design_matrix)
  fitmiRNA_1 <- eBayes(fitmiRNA_1)
  tablefitmiRNA_1 <- topTable(fitmiRNA_1, coef = 2, , number=Inf)
#List 2 miRNA Differential Analysis ----
  #getting the EGFR data from the mRNA data and then sorting
  lowQuartile <- quantile(as.vector(mRNACommon[which(rownames(mRNACommon) == 'EGFR'),]), probs = 0.25)
  highQuartile <- quantile(as.vector(mRNACommon[which(rownames(mRNACommon) == 'EGFR'),]), probs = 0.75)
  pos <- as.factor(ifelse(mRNACommon[rownames(mRNACommon) == 'EGFR', ] > highQuartile, 'pos', ifelse(mRNACommon[rownames(mRNACommon) == 'EGFR', ] < lowQuartile, 'neg', NA)))
  pos2 <- pos[!is.na(pos)]
  #creating the design matrix from the 25% and 75% indices
  design_matrix2 <- model.matrix(~pos2)
  #conducting the t-test using the limma package
  fitmiRNA_2 <- lmFit(normalizedmiRNACommon[,match(colnames(mRNACommon[,!is.na(pos)]), table=colnames(normalizedmiRNACommon))], design = design_matrix2)
  fitmiRNA_2 <- eBayes(fitmiRNA_2)
  tablefitmiRNA_2 <- topTable(fitmiRNA_2, coef = 2, , number=Inf)
#List 3 miRNA Differential Analysis ----
  #our data is the mRNACommon data
  #inputting the gene set matrix that contains the EGFR pathway genes
  #this dataset is downloaded from https://maayanlab.cloud/Harmonizome/dataset/Wikipathways+Pathways
  #egfrPathwaySet = read.table('/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/EGFR_gene_set.txt', header = TRUE, fill = TRUE)
  #egfrPathwaySet <- read_json(file = "/Users/harshitakumar/Documents/Research/Dhawan Lab/EGFR Project/BIOCARTA_EGF_PATHWAY.v2023.1.Hs.json")
  #egfrPathwaySet <- as.matrix(egfrPathwaySet[["BIOCARTA_EGF_PATHWAY"]][["geneSymbols"]])
  allGeneSets <- msigdbr(species = "Homo sapiens")
  #Getting Gene Sets
    geneSet1 <- allGeneSets[allGeneSets$gs_name == "HALLMARK_PI3K_AKT_MTOR_SIGNALING",]
    geneSet1 <- unique(geneSet1$gene_symbol)
    geneSet2 <- allGeneSets[allGeneSets$gs_name == "REACTOME_SIGNALING_BY_EGFR",]
    geneSet2 <- unique(geneSet2$gene_symbol)
    geneSet3 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGFR_SMRTE_PATHWAY",]
    geneSet3 <- unique(geneSet3$gene_symbol)
    geneSet4 <- allGeneSets[allGeneSets$gs_name == "REACTOME_EGFR_DOWNREGULATION",]
    geneSet4 <- unique(geneSet4$gene_symbol)
    geneSet5 <- allGeneSets[allGeneSets$gs_name == "WP_EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE",]
    geneSet5 <- unique(geneSet5$gene_symbol)
  library("GSVA")
  ssAnalysis <- gsva(mRNACommon, list(geneSet1), method = "ssgsea", mx.diff = FALSE )
  thr <- median(as.vector(ssAnalysis))
  pos <- as.factor(ifelse(ssAnalysis > thr, 'pos', 'neg'))
  design_matrix_3 <- model.matrix(~pos)
  #conducting the t-test using the limma package
  fitmiRNA_3 <- lmFit(normalizedmiRNACommon[,match(colnames(mRNACommon), table=colnames(normalizedmiRNACommon))],design = design_matrix_3)
  fitmiRNA_3 <- eBayes(fitmiRNA_3)
  tablefitmiRNA_3 <- topTable(fitmiRNA_3, coef = 2, , number=Inf)
#List 4 miRNA Differential Analysis ----
  lowQuartile <- quantile(ssAnalysis, probs = 0.25)
  highQuartile <- quantile(ssAnalysis, probs = 0.75)
  pos <- as.factor(ifelse(ssAnalysis > highQuartile, 'pos', ifelse(ssAnalysis < lowQuartile, 'neg', NA)))
  pos2 <- pos[!is.na(pos)]
  #creating the design matrix from the 25% and 75% indices
  design_matrix_4 <- model.matrix(~pos2)
  #conducting the t-test using the limma package
  fitmiRNA_4 <- lmFit(normalizedmiRNACommon[,match(colnames(mRNACommon[,!is.na(pos)]), table=colnames(normalizedmiRNACommon))], design = design_matrix_4)
  fitmiRNA_4 <- eBayes(fitmiRNA_4)
  tablefitmiRNA_4 <- topTable(fitmiRNA_4, coef = 2, , number=Inf)
#miRNA survival correlation Analysis  ----
  miRNAList <- 'hsa-miR-219-2-3p'
  miRNAData <- normalizedmiRNACommon[which(rownames(normalizedmiRNACommon) == miRNAList),]
  medianmiRNA <- median(as.vector(miRNAData))
  pos <- as.factor(ifelse(miRNAData > medianmiRNA, 'pos', 'neg'))
  miRNADesignMatrix <- model.matrix(~pos)
  clinicalmiRNAData <- lmFit(commonClinical[match(commonClinical$CGGA_ID, names(miRNAData)),7], design = miRNADesignMatrix)
  clinicalmiRNAData <- eBayes(clinicalmiRNAData)
  tableClinicalmiRNAData <-  topTable(clinicalmiRNAData, coef = 2, , number=Inf)
  View(tableClinicalmiRNAData)
  #creating Kaplan Meier Curve
    surv_object <- Surv(time = commonClinical$OS, event = commonClinical$Censor)
    km_fit <- survfit(surv_object ~ miRNADesignMatrix[,2], data = commonClinical)
    plot(km_fit,  lty = c("solid", "dashed"), col = c("black", "grey"), main =  c(miRNAList,"Kaplan-Meier Survival Curve"), xlab = "Time", ylab = "Survival Probability")
    legend("topright", c("High", "Low"), lty = c("solid", "dashed"), col = c("black", "grey"))
    
#Rank Product Analysis ----
  # #Dividing lists by downregulated and upregulated
  #   list1Up <- tablefitmiRNA_1[tablefitmiRNA_1$logFC>0, ]
  #   list1Down <- tablefitmiRNA_1[tablefitmiRNA_1$logFC<0, ]
  #   list2Up <- tablefitmiRNA_1[tablefitmiRNA_2$logFC>0, ]
  #   list2Down <- tablefitmiRNA_1[tablefitmiRNA_2$logFC<0, ]
  #   list3Up <- tablefitmiRNA_1[tablefitmiRNA_3$logFC>0, ]
  #   list3Down <- tablefitmiRNA_1[tablefitmiRNA_3$logFC<0, ]
  #   list4Up <- tablefitmiRNA_1[tablefitmiRNA_4$logFC>0, ]
  #   list4Down <- tablefitmiRNA_1[tablefitmiRNA_4$logFC<0, ]
  #   colnames(list2Up) <- c("logFC","AveExpr","t2","P.Value","adj.P.Val","B")
  #   colnames(list3Up) <- c("logFC","AveExpr","t3","P.Value","adj.P.Val","B")
  #   colnames(list4Up) <- c("logFC","AveExpr","t4","P.Value","adj.P.Val","B")
  #   allListUp <- merge(merge(merge(list1Up, list2Up,by = 0), list3Up, by = 0), list4Up, by = 0)
  #   allListDown <- c(list1Down[1:20,],list2Down[1:20,],list3Down[1:20,],list4Down[1:20,])
  #union of these lists
    # miRNAList <- miRNACommon$microRNA_ID
    # positionLists <- data.frame(match(miRNAList, rownames(list1Up)), match(miRNAList, rownames(list2Up)), match(miRNAList, rownames(list3Up)), match(miRNAList, rownames(list4Up)))
    # rankProdUp <- RP.advance(positionLists, rep(1, length(miRNAList)),rep(1, length(miRNAList)), logged = FALSE)
    library("RankProd")
    common_ids <-Reduce('intersect', list(rownames(tablefitmiRNA_1), rownames(tablefitmiRNA_2), rownames(tablefitmiRNA_3), rownames(tablefitmiRNA_4)))
    rankProdFC <- RP(cbind(tablefitmiRNA_1$logFC[match(common_ids, table=rownames(tablefitmiRNA_1))], tablefitmiRNA_2$logFC[match(common_ids, table=rownames(tablefitmiRNA_2))], tablefitmiRNA_3$logFC[match(common_ids, table=rownames(tablefitmiRNA_3))], tablefitmiRNA_4$logFC[match(common_ids, table=rownames(tablefitmiRNA_4))]), rep(1,4))
    # rankProdTVal <- RP(cbind(tablefitmiRNA_1$t[match(common_ids, table=rownames(tablefitmiRNA_1))], tablefitmiRNA_2$t[match(common_ids, table=rownames(tablefitmiRNA_2))], tablefitmiRNA_3$t[match(common_ids, table=rownames(tablefitmiRNA_3))], tablefitmiRNA_4$t[match(common_ids, table=rownames(tablefitmiRNA_4))]), rep(1,4))
    common_ids <- t(as.data.frame(as.list(common_ids)))
    rankProdGenes <- topGene(rankProdFC,method="pval",logged = FALSE, gene.names = common_ids, cutoff = 0.05)
    View(rankProdGenes[["Table1"]])
    View(rankProdGenes[["Table2"]])
    # rankProdGenesTVal <- topGene(rankProdTVal,method="pval",logged = FALSE, gene.names = common_ids, cutoff = 0.05)
    # View(rankProdGenesTVal[["Table1"]])
    # View(rankProdGenesTVal[["Table2"]])
 # Comparison ----
    armandatatable1 <- read_excel("miRNAs Rank prod  (1).xlsx", sheet = 2)
    armandatatable2 <- read_excel("miRNAs Rank prod  (1).xlsx", sheet = 1)
