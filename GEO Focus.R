#GEO Datasets 
#Load appropriate libraries
library(TCGAbiolinks)
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


#Part A: Data collection and annotation 

#Step 1: Query, download, and prep GEO data sets 
GEOquery <- getGEO(
  GEO = "GSE90604",
  GSEMatrix = TRUE)

#Looking at metadata
metadata <- pData(phenoData(GEOquery[[1]]))
pData(GEOquery[[1]])
fData(GEOquery[[1]])
exprs(GEOquery[[1]])[1,]

gene_exprs <- exprs(GEOquery[[1]])
test_exprs <- log2(gene_exprs)
summary(gene_exprs) #checking if normalized already (they are normalize)
boxplot(gene_exprs, outline = F, main = "mRNA", ylab= "expression")
plot(density(gene_exprs[,1]), main="Density Plot of mRNA Data") #Visualize with Density plots 
for (i in 2:ncol(gene_exprs)) {
  lines(density(gene_exprs[,i]))}

miRNA_exprs <- exprs(GEOquery[[2]])
summary(miRNA_exprs) #checking if normalized already (they are normalized)
boxplot(miRNA_exprs, outline = F, main = "miRNA")
plot(density(miRNA_exprs[,1]), main="Density Plot of microRNA Data") #Visualize with Density plots 
for (i in 2:ncol(miRNA_exprs)) {
  lines(density(miRNA_exprs[,i]))}

colnames(gene_exprs) = pData(GEOquery[[1]])$title #Rename sample (column) names 
colnames(miRNA_exprs) = pData(GEOquery[[2]])$title #Rename sample (column) names 

#Visualize with heatmaps
genecorMatrix <- cor(gene_exprs, use="c")
annotation_col = dplyr::select(pData(GEOquery[[1]]), "gender:ch1", "location:ch1") #Gender and GBM location info
rownames(annotation_col) = pData(GEOquery[[1]])$title
pheatmap(genecorMatrix, annotation_col = annotation_col)

#Visualize with PCA plots 
pca <- prcomp(t(gene_exprs)) #transpose so gene list is in columns, samples are rownames 
sample_type = dplyr::select(pData(GEOquery[[1]]), "sample type:ch1") #sample type info
# Update the column names
sample_type[,1] <- gsub("fresh frozen", "", sample_type[,1])
rownames(sample_type) = pData(GEOquery[[1]])$title
pca = cbind(sample_type, pca$x)
rownames(pca) <- gsub("Glioblastoma Multiforme", "GBM", rownames(pca))
ggplot(pca, aes(x = PC1, y = PC2, col = pca[,1], label = rownames(pca))) + geom_point() + geom_text_repel()


#Differential Expression Analysis
design <- model.matrix(~ 0 + sample_type[,1])
cutoff <- median(gene_exprs) #median expression level
is_expressed <- gene_exprs > cutoff #table, label TRUE if gene expression above cutoff 
keep <- rowSums(is_expressed) > 11 #genes that are expressed above the median in more than 11 samples 
table(keep)
gene_exprs2 <- gene_exprs[keep,]
aw <- arrayWeights(gene_exprs2, design)
aw 

# Set the column names of the design matrix to the levels of `sample_type`
colnames(design) <- c("GBM_tumor", 
                       "healthy_brain", 
                       "fetal_human_astrocyte", 
                       "normal_human_astrocyte")

zero_variance_genes <- apply(gene_exprs2, 1, var) == 0
gene_exprs2 <- gene_exprs2[!zero_variance_genes, ]# Remove the genes with zero variance
aw <- arrayWeights(gene_exprs2, design) # Recalculate array weights


fit <- lmFit(gene_exprs2, design, weights = aw)
head(fit$coefficients)

# Define the contrast
contrasts <- makeContrasts(
  GBM_tumor - healthy_brain, 
  #GBM_tumor - fetal_human_astrocyte,
  #GBM_tumor - normal_human_astrocyte, 
  levels = design
)

fit2<- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

topTable(fit2)
topTable1 <- topTable(fit2, coef=1)

summary(decideTests(fit2))
table(decideTests(fit2))


#Annotate
anno <- fData(GEOquery[[1]])
head(anno)
anno <- dplyr::select(anno,ID,GB_ACC)
fit2$genes <- anno
topTable(fit2)


#Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()

#change according to  needs
p_cutoff <- 0.01
fc_cutoff <- 1
full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()


#filter(full_results1, GB_ACC == "NM_201283")
#full_results2 <- full_results1[ order(as.numeric(row.names(full_results1))), ]

table <- read.table("/Users/armanzadeh/Downloads/GBM GEO annotation 2.txt", header = TRUE, sep = "\t", fill = TRUE)
rownames(gene_exprs)[rownames(gene_exprs) %in% gene_symbols$affy_hugene_2_1_st] <- gene_symbols$external_gene_name
filtered_gene_exprs <- gene_exprs[rownames(gene_exprs) %in% egfr_genes, ]
which(full_results1$ID %in% table$ProbeName)
colnames(table) = c("ID", "GeneSymbols", "GeneNames", "GOTerms", "GemmaIDs", "NCBIids")
full_results1$ID = rownames(full_results1)
table <- table[,1:2]
merged_data <- merge(table, full_results1, by = "ID")



#Examine differentialy expressed miRNAs between high and low EGFR expressing samples
#Step 1: Determine which samples from mRNA dataset show high EGFR Expression
colnames(gene_exprs) <- gsub("Glioblastoma Multiforme", "GBM", colnames(gene_exprs) )
gene_exprs_t <- gene_exprs #labeled data with genes. non-t is numeric data only 
gene_exprs_t <- as.data.frame(gene_exprs_t)
gene_exprs_t$ID <- rownames(gene_exprs_t) 
gene_exprs_t <- merge(table, gene_exprs_t, by = "ID")  #change rownames to gene names based on annotation file 
rows_with_EGFR <- which(gene_exprs_t$GeneSymbols == "EGFR")
print(rows_with_EGFR)
cutoff <- median(gene_exprs[rows_with_EGFR,1:16])  #filter for EGFR median from GBM samples only (no healthy samples)
high_EGFR_exprs_samples = gene_exprs[,which(gene_exprs[rows_with_EGFR,1:16] > cutoff)]
low_EGFR_exprs_samples = gene_exprs[,which(gene_exprs[rows_with_EGFR,1:16] < cutoff)]
expression_type = c()
for(i in 1:length(gene_exprs[1,1:16])){
  if(gene_exprs[rows_with_EGFR,i] > cutoff){
    expression_type = c(expression_type, "High EGFR Expression")
  }
  else{
    expression_type = c(expression_type, "Low EGFR Expression")
  }
}

#Step X: DEA on miRNA

colnames(miRNA_exprs) <- gsub("Glioblastoma Multiforme", "GBM", colnames(miRNA_exprs))
miRNA_exprs_t = as.data.frame(miRNA_exprs)
miRNA_exprs_t$ID = rownames(miRNA_exprs)
miRNA_annotation_table <- read.table("/Users/armanzadeh/Downloads/GEO GBM miRNA annotation.txt", header = TRUE, sep = "\t", fill = TRUE)
miRNA_annotation_table = miRNA_annotation_table[,c(1,4)] #only want the ID and names from annotation file
miRNA_exprs_t <- merge(miRNA_annotation_table, miRNA_exprs_t, by = "ID")  #change rownames to miRNA names based on annotation file 
colnames(miRNA_annotation_table)[2] = 'miRNA'

#Visualize with heatmaps
miRNAcorMatrix <- cor(miRNA_exprs, use="c")
#annotation_col = dplyr::select(pData(GEOquery[[2]]))
#rownames(annotation_col) = pData(GEOquery[[1]])$title
pheatmap(miRNAcorMatrix) #, annotation_col = annotation_col)



# Set the column names of the design2 matrix to the levels of expression_type determined from mRNA
miRNA_exprs <- exprs(GEOquery[[2]]) #Reset miRNA data table
colnames(miRNA_exprs) = pData(GEOquery[[2]])$title #Rename sample (column) names 
colnames(miRNA_exprs) <- gsub("Glioblastoma Multiforme", "GBM", colnames(miRNA_exprs))
miRNA_exprs <- miRNA_exprs[,1:16] #isolate GBM samples, remove any healthy samples 

design2 <- model.matrix(~ 0 + expression_type)
colnames(design2) <- c("High_EGFR_expression", "Low_EGFR_expression")
cutoff <- median(miRNA_exprs) #median expression level
is_expressed <- miRNA_exprs > cutoff #table, label TRUE if miRNA expression above cutoff 
keep <- rowSums(is_expressed) > 7 #miRNAs that are expressed above the median in more than 7 (~50%) samples 
table(keep)
miRNA_exprs2 <- miRNA_exprs[keep,]
aw <- arrayWeights(miRNA_exprs2, design2)
aw 
zero_variance_miRNA <- apply(miRNA_exprs2, 1, var) == 0
miRNA_exprs2 <- miRNA_exprs2[!zero_variance_miRNA, ]# Remove the miRNAs with zero variance
aw <- arrayWeights(miRNA_exprs2, design2) # Recalculate array weights
aw
fit <- lmFit(miRNA_exprs2, design2, weights = aw) #make sure to update design variable
head(fit$coefficients)

# Define the contrast
contrasts <- makeContrasts(
  High_EGFR_expression - Low_EGFR_expression,
  levels = design2 #make sure to update design variable
)

fit2<- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

topTable(fit2)
topTable1 <- topTable(fit2, coef=1)

#Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()

#change according to  needs
p_cutoff <- 0.05
fc_cutoff <- 1
full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

#test = full_results1 %>% 
 # mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) 

#Make DEA table for miRNA 
full_results1 <- cbind(ID = rownames(full_results1), full_results1)
merged_miRNA_data <- merge(miRNA_annotation_table, full_results1, by = "ID") #only 474/3250 in annotation table? need to find a new annotation table
sig_merged_miRNA_data <- merged_miRNA_data[which(merged_miRNA_data$P.Value < 0.05),] #basically just miRNAs that have a label in our dataset, and are significant (very few)

expression_type_vis = as.data.frame(expression_type, row.names = colnames(miRNA_exprs))






#check mir10 expression "hsa-miR-10" 
mir10_IDs <- c("20500438", "20500440", "20500441")
mir10_rows <- which(full_results1$ID %in% mir10_IDs)
full_results1[mir10_rows,]





#ssGSEA 
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

ssgsea_input <- gene_exprs[,1:16]
ssgsea_input <- as.data.frame(ssgsea_input)
ssgsea_input$ID <- rownames(ssgsea_input)
ssgsea_input <- merge(table, ssgsea_input, by = "ID")
for (i in 1:length(ssgsea_input[,1])){
  if(ssgsea_input[i,2] == ""){
    ssgsea_input[i,2] = paste0("NA", i)
  }
}
table(duplicated(ssgsea_input$GeneSymbols))
duplicates <- which(duplicated(ssgsea_input$GeneSymbols))
ssgsea_input <- ssgsea_input[-duplicates,]
table(duplicated(ssgsea_input$GeneSymbols))
rownames(ssgsea_input) = ssgsea_input$GeneSymbols
ssgsea_input <- ssgsea_input[,-c(1,2)]


groupGeneSet <- list(geneSet1, geneSet2, geneSet3, geneSet4, geneSet5)
names(groupGeneSet) <- c("HALLMARK_PI3K_AKT_MTOR_SIGNALING", "REACTOME_SIGNALING_BY_EGFR", 
                        "BIOCARTA_EGFR_SMRTE_PATHWAY", "REACTOME_EGFR_DOWNREGULATION",
                        "WP_EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE")

ssgsea <- ssgsea(inmat = ssgsea_input, groups = groupGeneSet, scale = TRUE, minsize = 10)
colnames(ssgsea) = colnames(ssgsea_input)
#test = z2p(ssgsea)

cutoff_ssgsea <- median(ssgsea)
high_EGFR_pathway_exprs_samples = ssgsea[,which(ssgsea > cutoff)]
low_EGFR_pathway_exprs_samples = ssgsea[,which(ssgsea < cutoff)]
pathway_type = c()

#for(i in 1:length(gene_exprs[1,1:16])){
#  if(gene_exprs[rows_with_EGFR,i] > cutoff){
#    expression_type = c(expression_type, "High EGFR Expression")
#  }
#  else{
#    expression_type = c(expression_type, "Low EGFR Expression")
#  }
#}



#Rank Prod 

  




































































###Scratch###



#MsigDB EGFR gene pathway list 
gene_set <- msigdbr()
egfr_genes <- gene_set[gene_set$gs_name == "HALLMARK_PI3K_AKT_MTOR_SIGNALING",] #filter for the specific gene set
egfr_genes <- unique(egfr_genes$gene_symbol) #isolate the gene names

# Convert probe IDs to gene symbols
table <- read.table("/Users/armanzadeh/Downloads/GBM GEO annotation 2.txt", header = TRUE, sep = "\t", fill = TRUE)
rownames(gene_exprs)[rownames(gene_exprs) %in% gene_symbols$affy_hugene_2_1_st] <- gene_symbols$external_gene_name
filtered_gene_exprs <- gene_exprs[rownames(gene_exprs) %in% egfr_genes, ]






#xx <- as.list(hugene21sttranscriptclusterALIAS2PROBE)
#gene_symbols <- mapIds(hugene21sttranscriptcluster.db, keys=rownames(gene_exprs), column="SYMBOL", keytype="PROBEID")

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Convert Affymetrix IDs to Ensembl gene IDs
gene_symbols <- getBM(attributes = c('affy_hugene_2_1_st_v1', 'ensembl_gene_id', 'external_gene_name'),
                      filters = "affy_hugene_2_1_st_v1",
                      values = rownames(gene_exprs),
                      mart = mart)
#rownames(gene_exprs) <- gene_symbols$external_gene_name
rownames(gene_exprs)[rownames(gene_exprs) %in% gene_symbols$affy_hugene_2_1_st] <- gene_symbols$external_gene_name
filtered_gene_exprs <- gene_exprs[rownames(gene_exprs) %in% egfr_genes, ]



#getGene( id = rownames(gene_exprs), type = "affy_hugene_2_1_st_v1", mart)


# Convert probe set names to miRNA IDs (probably not needed for now?)
miRNA_IDs <- mapIds(miRBaseVersions.db, keys=rownames(miRNA_exprs), column="miRNA_ID", keytype="PROBEID")






#GEOdat = GEOquery[["GSE90604-GPL17692_series_matrix.txt.gz"]]@assayData[["exprs"]]
#GEOdat_tumor = GEOdat[,1:16] #tumor samples from GEO dataset 
#GEOdat_normal = GEOdat[,17:25] #not to be used for now 























