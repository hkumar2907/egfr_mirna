#fold change stuff ----
#Geo dataset fold change check
# egfrPatients <- miRNA_exprs_t[,order(gene_exprs[which(rownames(gene_exprs) == 'EGFR'),])]
# miRofInterest <- egfrPatients[rownames(egfrPatients) == 'hsa-miR-204-5p', ]
# n <- ncol(miRofInterest)
# foldChange <- mean(as.numeric(miRofInterest[,(n/2+1):n]))/mean(as.numeric(miRofInterest[,1:(n/2)]))
library(ggplot2)

miRPlotEGFR <- function(miRName, mRNATable, miRTable, geneSetName) {
  egfrSpot <- order(mRNATable[which(rownames(mRNATable) == 'EGFR'),])
  egfrmRNA <- mRNATable[which(rownames(mRNATable) == 'EGFR'), egfrSpot]
  miRofInterest <- miRTable[rownames(miRTable) == miRName,egfrSpot]
  dataEGFR <- data.frame(egfr = list(t(egfrmRNA)), miR = list(t(miRofInterest)))
  #dataEGFR <- data.frame(egfr = list(t(egfrmRNA[1,])), miR = list(t(miRofInterest[1,])))
  colnames(dataEGFR) <- c("EGFR", "miR")
  ggplot(data = dataEGFR,aes(x=miR, y=EGFR)) + geom_point() + labs(title=paste(miRName, geneSetName)) + geom_text(aes(label = rownames(dataEGFR)), nudge_x = 0.1) + geom_smooth(method = "lm", se = FALSE, color = "red") 
}

miRPlotEGFR('hsa-miR-204', mRNACommon, normalizedmiRNACommon, "CGGA")

#venn diagram ----

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, category.names = c("CGGA" , "Geo" , "Study"),filename = NULL, ...)
  grid.draw(venn_object)
}

library(VennDiagram)
CGGA_t2 <- rownames(rankProdGenes[["Table2"]])[1:50]
Geo_t2 <- rownames(rankProdGeo[["Table2"]])[1:50]
WES_t2 <- rownames(rankProdstudy[["Table2"]])[1:50]
venn.diagram(x = list(G, Geo_t2, WES_t2),category.names = c("CGGA" , "Geo" , "Study"),filename = "table2venn1.png")
display_venn(list(CGGA_t2, Geo_t2, WES_t2))

CGGA_t1 <- rownames(rankProdGenes[["Table1"]])[1:50]
Geo_t1 <- rownames(rankProdGeo[["Table1"]])[1:50]
WES_t1 <- rownames(rankProdstudy[["Table1"]])[1:50]
venn.diagram(x = list(CGGA_t1, Geo_t1, WES_t1),category.names = c("CGGA" , "Geo" , "Study"),filename = "table2venn3.png",  col = "black", output = TRUE, overlap = TRUE)

library("clipr")
write_clip(CGGA_t1)
write_clip(CGGA_t2)
write_clip(WES_t1)
write_clip(WES_t2)
write_clip(Geo_t1)
write_clip(Geo_t2)
# miRNA Differential Analysis WRONG ---- 
# examining differentially expressed miRNAs between high and low EGFR expressing samples
# just getting the data
EGFRData <- mRNACommon[which(mRNACommon$Gene_Name == 'EGFR'),-1]
EGFR_indices <- order(EGFRData)
#creating the design matrix
design_matrix <- model.matrix(~ as.factor(rep(c(rep(1, 79), rep(0, 80)))))
#finding out which samples express more EGFR than others
EGFRData <- EGFRData[EGFR_indices]
design_matrix <- design_matrix[EGFR_indices,]
lowEGFR <- EGFRData[1, 1:(ncol(EGFRData)/2)]
highEGFR <- EGFRData[1, (ncol(EGFRData)/2 + 1): ncol(EGFRData)]
#splitting the miRNA samples based on high and low egfr
lowEGFRgene = intersect(colnames(miRNACommon), colnames(highEGFR))
#lowEGFRgene = miRNACommon[, c(1, which(colnames(miRNACommon) %in% lowEGFRgene))]
highEGFRgene = intersect(colnames(miRNACommon), colnames(lowEGFR))
#highEGFRgene = miRNACommon[, c(1, which(colnames(miRNACommon) %in% highEGFRgene))]
#conducting ttest
fitmiRNA <- lmFit(normalizedmiRNACommon, design = design_matrix)
fitmiRNA <- eBayes(fitmiRNA)
tablefitmiRNA <- topTable(fitmiRNA, coef = 2, , number=Inf)
#visualizing the p-Values
hist(tablefitmiRNA$P.Value, breaks = 20, col = "blue", main = "Histogram of p-values", xlab = "p-value")
#getting the genenames
common_rows <- intersect(rownames(tablefitmiRNA), rownames(miRNACommon))
miRNACommon <- miRNACommon[common_rows,]
#normalizing the microRNA EGFR data
normalizedmiRNACommon = normalizeBetweenArrays(miRNACommon[, 2:ncol(miRNACommon)], method = "cyclicloess")
boxplot(normalizedmiRNACommon)