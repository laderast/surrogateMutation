plotAllGraphs <-
function(surrogateTable, surrogateResults, intome) {
  source("~/Code/surrogate/plot-results.R")
  for(i in 1:nrow(surrogateTable)){
    NodeName <- as.character(surrogateTable[i,"NodeName"])
    Sample <- as.character(surrogateTable[i,"Sample"])
    GeneName <- as.character(surrogateTable[i,"Gene"])
    print(paste(NodeName,Sample))
    renderSubNetSimple(NodeName, Sample, GeneName, intome, gisticCopyCalls=NULL, 
                       resultObj=surrogateResults)
  }
}
