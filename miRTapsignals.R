library(miRNAtap)
library(miRNAtap.db)
library(topGO)
library(org.Hs.eg.db)
library(annotate)
library(ComplexHeatmap)
library(circlize)
library(biclust)

#primary analysis ----
#these are mirs that are good targets from table 2
  #mirs <- list("30a", "24", "148a", "223", "424", "221", "34a", "30e", "195", "338")
  #mirs <- list("30a", "24", "148a", "223", "424", "221", "34a", "30e", "195")
  mirs <- list("30a", "24", "148a", "223", "424", "221", "34a", "30e", "195", "338", "152", "21", "222", "140", "708", "92b", "320a", "10b", "484", "328", "423")
  mirs <- lapply(mirs, function(x) paste("mir-", x, sep = ""))
  
  table2Targets <- names(table2GBM)[table2GBM > 1] #getting the targets from table 2 from GBM list (all atleast common in two datasets)
  mirs <- str_sub(table2Targets, 5, str_length(table2Targets))#setting the values to mirs and also removing the hsa- from the start
  #from top 50 targets
  predictedTargets <- list()
  for (item in mirs) {
    print(item)
    predictedTargets[[item]] <- unique(getSYMBOL(rownames(getPredictedTargets(item, species = 'hsa', method = 'geom', min_src = 2)), data='org.Hs.eg'))
  }
#loading analysis ----
  #getting genes that are in the top 75% expression of all genes (finding vals in total that are in bottom 25%, and then removing genes that have those for more than 50% of the samples)
  filteringGenes <- function(mrna) {
    G_filtered <- abs(mrna)
    G_quantile <- quantile(as.vector(G_filtered))[["25%"]]
    G_filtered[G_filtered < G_quantile] <- NA
    naSums <- sort(rowSums(is.na(G_filtered)), decreasing = TRUE)
    naSums <- naSums[naSums <= (0.5*ncol(G_filtered))]
    G_filtered <- mrna[names(naSums), ]
    return(rownames(G_filtered))
  }
  cggaFilter <- filteringGenes(mRNACommon) #CGGA
  gene_exprs <- as.matrix(gene_exprs)
  mode(gene_exprs) = "numeric"
  geoFilter <- filteringGenes(gene_exprs) #Geo Dataset
  studyFilter <- filteringGenes(as.matrix(initial_patients_mRNA)) #Whole exome study PROBLEM
  
  #EGFR pathway
  geneSet1 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGF_PATHWAY",]
  geneSet1 <- unique(geneSet1$gene_symbol)
  geneSet4 <- allGeneSets[allGeneSets$gs_name == "HALLMARK_PI3K_AKT_MTOR_SIGNALING",]
  geneSet4 <- unique(geneSet4$gene_symbol)
  geneSet2 <- allGeneSets[allGeneSets$gs_name == "REACTOME_SIGNALING_BY_EGFR",]
  geneSet2 <- unique(geneSet2$gene_symbol)
  geneSet3 <- allGeneSets[allGeneSets$gs_name == "BIOCARTA_EGFR_SMRTE_PATHWAY",]
  geneSet3 <- unique(geneSet3$gene_symbol)
  G_EGFR <- union(union(union(geneSet1, geneSet4), geneSet3), geneSet2)
  combined <- union(geoFilter, cggaFilter)
#running through the the different mirs ----
  #length(predictedTargets
  heatmapList <- list()
  for (miR in 1:length(predictedTargets)) {
    #running through the actual list of predicted targets
    mirName <- names(predictedTargets[miR])
    print(mirName)
    #heatmapList[[mirName]] <- NaN * seq(length(predictedTargets[[miR]]))
    for (gene in 1:length(predictedTargets[[miR]])) {
      geneName <- predictedTargets[[miR]][[gene]]
      #print(predictedTargets[[miR]][[gene]])
      #genes only EGFR targets
      if (predictedTargets[[miR]][[gene]] %in% G_EGFR) { #EGFR target so a val of 1 given
        heatmapList[[mirName]][[geneName]] <- 1
      } else if (predictedTargets[[miR]][[gene]] %in% combined) { #not egfr target so a val of -1 given
        heatmapList[[mirName]][[geneName]]<- -1
      } else { #not a target at all so a val of 0 given
        heatmapList[[mirName]][[geneName]] <- 0
      }
    }
    names(heatmapList[[mirName]]) <- predictedTargets[[miR]]
   # print(predictedTargets[[miR]])
  }
#finding genes between the different lists ----
  genesUnion <- names(heatmapList[[1]])
  genesIntersect <- names(heatmapList[[1]])
  for (mir in 2:length(heatmapList)) {
    genesUnion <- union(genesUnion, names(heatmapList[[mir]]))
    genesIntersect <- intersect(genesIntersect, names(heatmapList[[mir]]))
    #print(length(genesUnion))
  }  
#creating dataframe ----
  unionDataframe <- as.data.frame(matrix(data = NA, nrow = length(genesUnion), ncol = length(heatmapList)), row.names = genesUnion)
  colnames(unionDataframe) <- mirs
  for (mir in 1:length(heatmapList)) {
    num <- match(names(heatmapList[[mir]]), genesUnion, nomatch = 0)
    unionDataframe[num, mir] <- unname(unlist(heatmapList[[mir]]))
  }
  #sorting datasets to have genes most prevalent at the top
  unionDataframe <- unionDataframe[order(rowSums(abs(unionDataframe), na.rm = TRUE), decreasing = TRUE), ]
#blue over red calculation + Plotting ----
  calcFrame <- as.data.frame(matrix(0, nrow = 6, ncol = length(heatmapList)), row.names = c("EGFR", "non-egfr", "none", "Sum egfr + non-egfr", "All combined", "ratio"))
  colnames(calcFrame) <- mirs
  for (mir in 1:length(heatmapList)) {
    calcFrame[1, mir] <- length(which(unionDataframe[mir] == 1))
    calcFrame[2, mir] <- length(which(unionDataframe[mir] == -1))
    calcFrame[3, mir] <- length(which(unionDataframe[mir] == 0))
    calcFrame[4, mir] <- calcFrame[1, mir] + calcFrame[2, mir]
    calcFrame[5, mir] <- calcFrame[4, mir] + calcFrame[3, mir]
    calcFrame[6, mir] <- calcFrame[1, mir]/calcFrame[2, mir]
  }
  sorting <- order(calcFrame[5,], decreasing = TRUE)
  calcFrame <- calcFrame[,sorting]
  calcFrame <- format(calcFrame, digits = 1, scientific=F) #making it pretty
  unionDataframe <- unionDataframe[,sorting]
  unionDataframe <- unionDataframe[order(rowSums(unionDataframe, na.rm = TRUE), decreasing = TRUE), ]
  unionDataframe2 <- unionDataframe[rowSums(abs(unionDataframe), na.rm = TRUE) > 4, ]
  unionDataframe2 <- unionDataframe[order(rowSums(is.na(unionDataframe))), ] #sorting based on how many nas present in the data
  col_fun = colorRamp2(c(-1, 0, 1), c("red", "black", "blue")) #for combined
  Heatmap(unionDataframe[1:250,], col = col_fun, na_col = "white", border = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, heatmap_legend_param = list(title = "EGFR and MicroRNAs", at = c(1, 0, -1), labels = c("EGFR Genes", "none", "Expressed Genes")))
  Heatmap(unionDataframe2[1:500,], col = col_fun, na_col = "white", border = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, heatmap_legend_param = list(title = "EGFR and MicroRNAs", at = c(1, 0, -1), labels = c("EGFR Genes", "none", "Expressed Genes")))
  
#creating EGFR graph network?
  #all mirs with only egfr gene targets
  miREGFRDF <- list()
  for (miR in 1:length(predictedTargets)) {
    allGenes <- intersect(predictedTargets[[miR]], G_EGFR) #finding out which genes for that miR are also EGFR genes
    miREGFRDF <- rbind(miREGFRDF, cbind(rep(names(predictedTargets)[miR], length(allGenes)), allGenes))
  }
  