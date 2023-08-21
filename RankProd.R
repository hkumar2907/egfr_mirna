#Rank Prod 
  # List 1 = significant miRNA’s between high and low EGFR patients (using the median)
  # List 2 = top 25% patients and bottom 25% patients instead
  # List 3 = median like the first list, but instead we use the ssgsea score to split patients
  # List 4 = ssgsea score, and like list 2, splits by 25% and 25% patients

#Step 0: Add Gene symbols and filter for EGFR expression across samples 
gene_exprs_t <- gene_exprs #labeled data with genes. non-t is numeric data only 
gene_exprs_t <- as.data.frame(gene_exprs_t)
gene_exprs_t$ID <- rownames(gene_exprs_t) 
gene_exprs_t <- merge(table, gene_exprs_t, by = "ID")  #change rownames to gene names based on annotation file 
rows_with_EGFR <- which(gene_exprs_t$GeneSymbols == "EGFR")
print(rows_with_EGFR)

#LIST 1
#Step 1.1: Create EGFR threshold from GBM samples only 
cutoff_List1 <- median(gene_exprs[rows_with_EGFR,1:16])  #filter for EGFR median from GBM samples only (no healthy samples)
high_EGFR_exprs_samples = gene_exprs[,which(gene_exprs[rows_with_EGFR,1:16] > cutoff_List1)]
low_EGFR_exprs_samples = gene_exprs[,which(gene_exprs[rows_with_EGFR,1:16] < cutoff_List1)]

#Step 1.2: Separate 16 GBM samples into low or high based on thresholds 
expression_type = c()
for(i in 1:length(gene_exprs[1,1:16])){
  if(gene_exprs[rows_with_EGFR,i] > cutoff_List1){
    expression_type = c(expression_type, "High EGFR Expression")
  }
  else{
    expression_type = c(expression_type, "Low EGFR Expression")
  }
}

#Step 1.3: DEA - miRNAs differentialy expressed based on EGFR expression high vs low 
miRNA_exprs_GBM <- miRNA_exprs[,1:16]
#is_expressed <- miRNA_exprs_GBM > median(miRNA_exprs_GBM) #table, label TRUE if miRNA expression above cutoff 
#keep <- rowSums(is_expressed) > 7 #miRNAs that are expressed above the median in more than 7 (~50%) samples 
#table(keep)
#miRNA_exprs_GBM <- miRNA_exprs_GBM[keep,] #remove any miRNAs that are not diff. expressed in greater than 7 samples 
zero_variance_miRNA <- apply(miRNA_exprs_GBM, 1, var) == 0
miRNA_exprs_GBM <- miRNA_exprs_GBM[!zero_variance_miRNA, ]# Remove the miRNAs with zero variance

design_list1 <- model.matrix(~ 0 + expression_type)
colnames(design_list1) <- c("High_EGFR_expression", "Low_EGFR_expression")

aw <- arrayWeights(miRNA_exprs_GBM, design_list1) # Recalculate array weights
aw
fit <- lmFit(miRNA_exprs_GBM, design_list1, weights = aw) #make sure to update design variable
head(fit$coefficients)


contrasts <- makeContrasts(
  High_EGFR_expression - Low_EGFR_expression,
  levels = design_list1) #make sure to update design variable

fit2<- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable1 <- topTable(fit2, coef=1)

#Step 1.4: Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()

p_cutoff <- 0.05
fc_cutoff <- 1.0
full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

#Step 1.5: Merge miRNA tables and view significant ones 
full_results1 <- cbind(ID = rownames(full_results1), full_results1)
merged_miRNA_data <- merge(miRNA_annotation_table, full_results1, by = "ID") #only 474/3250 in annotation table? need to find a new annotation table
sig_merged_miRNA_data <- merged_miRNA_data[which(merged_miRNA_data$P.Value < 0.05),] #basically just miRNAs that have a label in our dataset, and are significant (very few)

#Step 1.6: collect FC values 
list1_fullFCs <- merged_miRNA_data[,2:3]
list1_sigFCs <- sig_merged_miRNA_data[,2:3]



#LIST 2
#Step 2.1: Create cutoffs for top25% and bottom25%, filter out samples
cutoff_List2_q3 <- as.numeric(summary(gene_exprs[rows_with_EGFR,1:16])[5])  #filter for EGFR q1 from GBM samples only (no healthy samples)
cutoff_List2_q1 <- as.numeric(summary(gene_exprs[rows_with_EGFR,1:16])[2])  #filter for EGFR q3 from GBM samples only (no healthy samples)
topQuartile_EGFR_exprs_samples = gene_exprs[,which(gene_exprs[rows_with_EGFR,1:16] > cutoff_List2_q3)]
bottomQuartile_EGFR_exprs_samples = gene_exprs[,which(gene_exprs[rows_with_EGFR,1:16] < cutoff_List2_q1)]

#Step 2.2: Reset expression_type. Create list for design_list
expression_type = c("High EGFR Expression", "High EGFR Expression", "High EGFR Expression", "High EGFR Expression",
                    "Low EGFR Expression", "Low EGFR Expression", "Low EGFR Expression", "Low EGFR Expression")
#for(i in 1:length(gene_exprs[1,1:16])){
#  if(gene_exprs[rows_with_EGFR,i] > cutoff_List2_q3){
#    expression_type = c(expression_type, "High EGFR Expression")
#  }
#  if(gene_exprs[rows_with_EGFR,i] < cutoff_List2_q1){
#    expression_type = c(expression_type, "Low EGFR Expression")
#  }
#}


#Step 2.3: DEA - miRNAs differentialy expressed based on EGFR expression high vs low 
miRNA_exprs_GBM <- miRNA_exprs[,c(which(gene_exprs[rows_with_EGFR,1:16] > cutoff_List2_q3), which(gene_exprs[rows_with_EGFR,1:16] < cutoff_List2_q1))]
#is_expressed <- miRNA_exprs_GBM > median(miRNA_exprs_GBM) #table, label TRUE if miRNA expression above cutoff 
#keep <- rowSums(is_expressed) > 7 #miRNAs that are expressed above the median in more than 7 (~50%) samples 
#table(keep)
#miRNA_exprs_GBM <- miRNA_exprs_GBM[keep,] #remove any miRNAs that are not sig. expressed in greater than 7 samples 
zero_variance_miRNA <- apply(miRNA_exprs_GBM, 1, var) == 0
miRNA_exprs_GBM <- miRNA_exprs_GBM[!zero_variance_miRNA, ]# Remove the miRNAs with zero variance

design_list1 <- model.matrix(~ 0 + expression_type)
colnames(design_list1) <- c("High_EGFR_expression", "Low_EGFR_expression")

aw <- arrayWeights(miRNA_exprs_GBM, design_list1) # Recalculate array weights
aw
fit <- lmFit(miRNA_exprs_GBM, design_list1, weights = aw) #make sure to update design variable
head(fit$coefficients)


contrasts <- makeContrasts(
  High_EGFR_expression - Low_EGFR_expression,
  levels = design_list1) #make sure to update design variable

fit2<- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable1 <- topTable(fit2, coef=1)

#Step 2.4: Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()

p_cutoff <- 0.05
fc_cutoff <- 1
full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

#Step2.5
full_results1 <- cbind(ID = rownames(full_results1), full_results1)
merged_miRNA_data <- merge(miRNA_annotation_table, full_results1, by = "ID") #only 474/3250 in annotation table? need to find a new annotation table
sig_merged_miRNA_data <- merged_miRNA_data[which(merged_miRNA_data$P.Value < 0.05),] #basically just miRNAs that have a label in our dataset, and are significant (very few)

#Step 2.6: collect FC values 
list2_fullFCs <- merged_miRNA_data[,2:3]
list2_sigFCs <- sig_merged_miRNA_data[,2:3]





#LIST 3
#Step 3.1: create threshold for ssgsea 
ssgsea_1 <- ssgsea[1,] #just HALLMARK_PI3K_AKT_MTOR_SIGNALING pathway 
cutoff_List3 <- median(ssgsea_1)
high_EGFR_pathway_exprs_samples = ssgsea_1[which(ssgsea_1 > cutoff_List3)]
low_EGFR_pathway_exprs_samples = ssgsea_1[which(ssgsea_1 < cutoff_List3)]

#Step 3.2: Separate 16 GBM samples into low or high based on thresholds 
expression_type = c()
for(i in 1:16){
  if(ssgsea_1[i] > cutoff_List3){
    expression_type = c(expression_type, "High EGFR Pathway Expression") #not technically EGFR, its the pathway
  }
  else{
    expression_type = c(expression_type, "Low EGFR Pathway Expression")
  }
}


#Step 3.3: DEA - miRNAs differentialy expressed based on EGFR expression high vs low 
miRNA_exprs_GBM <- miRNA_exprs[,1:16]
#is_expressed <- miRNA_exprs_GBM > median(miRNA_exprs_GBM) #table, label TRUE if miRNA expression above cutoff 
#keep <- rowSums(is_expressed) > 7 #miRNAs that are expressed above the median in more than 7 (~50%) samples 
#table(keep)
#miRNA_exprs_GBM <- miRNA_exprs_GBM[keep,] #remove any miRNAs that are not diff. expressed in greater than 7 samples 
zero_variance_miRNA <- apply(miRNA_exprs_GBM, 1, var) == 0
miRNA_exprs_GBM <- miRNA_exprs_GBM[!zero_variance_miRNA, ]# Remove the miRNAs with zero variance

design_list1 <- model.matrix(~ 0 + expression_type)
colnames(design_list1) <- c("High_EGFR_pathway_expression", "Low_EGFR_pathway_expression")

aw <- arrayWeights(miRNA_exprs_GBM, design_list1) # Recalculate array weights
aw
fit <- lmFit(miRNA_exprs_GBM, design_list1, weights = aw) #make sure to update design variable
head(fit$coefficients)


contrasts <- makeContrasts(
  High_EGFR_pathway_expression - Low_EGFR_pathway_expression,
  levels = design_list1) #make sure to update design variable

fit2<- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable1 <- topTable(fit2, coef=1)

#Step 3.4: Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()

p_cutoff <- 0.05
fc_cutoff <- 1
full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

#Step 3.5
full_results1 <- cbind(ID = rownames(full_results1), full_results1)
merged_miRNA_data <- merge(miRNA_annotation_table, full_results1, by = "ID") #only 474/3250 in annotation table? need to find a new annotation table
sig_merged_miRNA_data <- merged_miRNA_data[which(merged_miRNA_data$P.Value < 0.05),] #basically just miRNAs that have a label in our dataset, and are significant (very few)

#Step 3.6: collect FC values 
list3_fullFCs <- merged_miRNA_data[,2:3]
list3_sigFCs <- sig_merged_miRNA_data[,2:3]



#LIST 4
#Step 4.1: reset thresholds for list 4
cutoff_List4_q3 <- as.numeric(summary(ssgsea_1)[5])  
cutoff_List4_q1 <- as.numeric(summary(ssgsea_1)[2])  
topQuartile_EGFR_pathway_exprs_samples = ssgsea_1[which(ssgsea_1 > cutoff_List4_q3)]
bottomQuartile_EGFR_pathway_exprs_samples = ssgsea_1[which(ssgsea_1 < cutoff_List4_q1)]

#Step 4.2: make expression_type table for DEA
expression_type = c("High EGFR Expression", "High EGFR Expression", "High EGFR Expression", "High EGFR Expression",
                    "Low EGFR Expression", "Low EGFR Expression", "Low EGFR Expression", "Low EGFR Expression")


#Step 4.3: DEA - miRNAs differentialy expressed based on EGFR SSGSEA expression high vs low 
miRNA_exprs_GBM <- miRNA_exprs[,c(which(ssgsea_1 > cutoff_List4_q3), 
                                  which(ssgsea_1 < cutoff_List4_q1))]

zero_variance_miRNA <- apply(miRNA_exprs_GBM, 1, var) == 0
miRNA_exprs_GBM <- miRNA_exprs_GBM[!zero_variance_miRNA, ]# Remove the miRNAs with zero variance

design_list1 <- model.matrix(~ 0 + expression_type)
colnames(design_list1) <- c("High_EGFR_pathway_expression", "Low_EGFR_pathway_expression")

aw <- arrayWeights(miRNA_exprs_GBM, design_list1) # Recalculate array weights
aw
fit <- lmFit(miRNA_exprs_GBM, design_list1, weights = aw) #make sure to update design variable
head(fit$coefficients)


contrasts <- makeContrasts(
  High_EGFR_pathway_expression - Low_EGFR_pathway_expression,
  levels = design_list1) #make sure to update design variable

fit2<- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable1 <- topTable(fit2, coef=1)

#Step 4.4: Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()

p_cutoff <- 0.05
fc_cutoff <- 1
full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

#Step 4.5
full_results1 <- cbind(ID = rownames(full_results1), full_results1)
merged_miRNA_data <- merge(miRNA_annotation_table, full_results1, by = "ID") #only 474/3250 in annotation table? need to find a new annotation table
sig_merged_miRNA_data <- merged_miRNA_data[which(merged_miRNA_data$P.Value < 0.05),] #basically just miRNAs that have a label in our dataset, and are significant (very few)

#Step 4.6: collect FC values 
list4_fullFCs <- merged_miRNA_data[,2:3]
list4_sigFCs <- sig_merged_miRNA_data[,2:3]




#RankProd 
#Step 1: Merge all 4 lists 
mergedList = merge(list1_fullFCs, list2_fullFCs, by = "miRNA") 
mergedList = merge(mergedList, list3_fullFCs, by = "miRNA") 
mergedList = merge(mergedList, list4_fullFCs, by = "miRNA") 
colnames(mergedList) = c("miRNA", "List1", "List2", "List3", "List4")
rownames(mergedList) = mergedList$miRNA
mergedList = mergedList[,-1]

#Step 2: RP analysis 
RPresults = RP(mergedList, cl = rep(1,4))
plotRP(RPresults)

rankProdGenes <- topGene(RPresults, cutoff = 0.05, method= "pval", gene.names = rownames(mergedList))
View(rankProdGenes[["Table1"]])
View(rankProdGenes[["Table2"]])

rankProdtable2 <- rankProdGenes[["Table2"]]
write.csv(rank)


Harsh_Table2 <- read_excel("/Users/armanzadeh/Downloads/Correct rank prod.xlsx", sheet = 2)
Harshita_miRNAs <- Harsh_Table2$name
Armans_miRNAs <- rownames(rankProdtable2)
which(Harshita_miRNAs %in% Armans_miRNAs)
which(Armans_miRNAs %in% Harshita_miRNAs)

#223, 27a, 34a, 30a, 204





