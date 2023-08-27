#Crosstalkr ----
library("crosstalkr")
library(dplyr)
library(STRINGdb)
string_db <-  STRINGdb$new( version="11.5", species=9606,score_threshold=200, network_type="full", input_directory="")
exampleMap <- string_db$map( rankProdGenes[["Table1"]], "gene.index", removeUnmappedRows = TRUE )
egfrdb <- string_db$get_neighbors('egfr')
exp <- tablefitmiRNA_1$logFC
names(exp) <- rownames(tablefitmiRNA_1)
g <- gfilter(method = "value",g=string_db, cache = NULL, val=exp, val_name = "expression", use_ppi = FALSE, desc = TRUE,n=100)
