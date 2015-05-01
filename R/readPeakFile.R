readPeakFile <-
function(gisticPeakFile){
  gisticFile <- read.delim(gisticPeakFile, skip = 3, stringsAsFactors=FALSE) 
  gisticFile <- gisticFile[,-1]
  
  gisticList <- as.list(gisticFile)
  gisticList <- lapply(gisticList, function(x){ifelse(nchar(x)>0,x,NA)})  
  gisticList <- lapply(gisticList, function(x){na.omit(x)})
  gisticGenes <- do.call("c", gisticList)
  gisticGenes
}
