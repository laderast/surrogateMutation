buildGeneIntTable <-
function(intome){
  
  #build an ID-gene table for interactome (HRPD) dataset
  #IDs in interactome are GENE(entrezid), need to map to gene
  geneIntTable <- NULL
  for(gene in names(interactome@edgeL)){
    geneName <- strsplit(gene,"\\(")[[1]][1]
    entrezID <- strsplit(gene,"\\(")[[1]][2]
    entrezID <- gsub(")", "", entrezID)
    geneIntTable <- rbind(geneIntTable,c(gene,geneName, entrezID))
  }
  rownames(geneIntTable) <- geneIntTable[,2]
  colnames(geneIntTable) <- c("NodeName", "Gene", "EntrezID")
  
  geneIntTable <- data.frame(geneIntTable)
  return(geneIntTable)
}
