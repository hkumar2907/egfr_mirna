library(stringr)

#geoDataset
#Removing 3p and 5p
  geo2 <- removing3p5p(rownames(rankProdGeo[["Table2"]]))
  geo1 <- removing3p5p(rownames(rankProdGeo[["Table1"]]))
  cgga2 <- removing3p5p(rownames(rankProdGenes[["Table2"]]))
  cgga1 <- removing3p5p(rownames(rankProdGenes[["Table1"]]))
  wes2 <- removing3p5p(rownames(rankProdstudy[["Table2"]]))
  wes1 <- removing3p5p(rownames(rankProdstudy[["Table1"]]))

  num <- 50
  table2GBMrem <- c(c(geo2[1:num], cgga2[1:num]), wes2[1:num])
  table2GBMrem <- table(table2GBMrem)
  table1GBMrem <- c(c(geo1[1:num], cgga1[1:num]), wes1[1:num])
  table1GBMrem <- table(table1GBMrem)  
  
#Not removing 3p or 5p
  geo2 <- (rownames(rankProdGeo[["Table2"]]))
  geo1 <- (rownames(rankProdGeo[["Table1"]]))
  cgga2 <- (rownames(rankProdGenes[["Table2"]]))
  cgga1 <- (rownames(rankProdGenes[["Table1"]]))
  wes2 <- (rownames(rankProdstudy[["Table2"]]))
  wes1 <- (rownames(rankProdstudy[["Table1"]]))

  num <- 50
  table2GBM <- c(c(geo2[1:num], cgga2[1:num]), wes2[1:num])
  table2GBM <- table(table2GBM)
  table1GBM <- c(c(geo1[1:num], cgga1[1:num]), wes1[1:num])
  table1GBM <- table(table1GBM)

removing3p5p = function(rownames) {
  for (i in 1:length(rownames)) {
    if (grepl("-3p",rownames[i], fixed = TRUE) || grepl("-5p",rownames[i], fixed = TRUE)){
      rownames[i] <- str_sub(rownames[i], 1, str_length(rownames[i])-3)
    }
  }
  return(rownames)
}
