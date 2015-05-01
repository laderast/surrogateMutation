renderSubNetSimple <-
function(NodeName, sampleName, GeneName, intome=intome, 
                         gisticCopyCalls=NULL, resultObj,
                         fileOut=NULL){
  
  library(Rgraphviz)
  library(BioNet)
  
  #grab gistic status
  if(!is.null(gisticCopyCalls)){
  cellFrame <- gisticCopyCalls[,c("Gene.Symbol", "NodeName", sampleName)]
  }
  
  mutCopyFrame <- resultObj$mutCopyFrames[[sampleName]]
  cellResults <- resultObj$cellResults
  
  nodenet <- c(NodeName,intersect(as.character(graph::inEdges(NodeName, intome)[[1]]), 
                                  as.character(mutCopyFrame$NodeName)))  
  
  #grab those nodes with mutations and those nodes with copy alterations
  brcamut <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$mutations!=0,"NodeName"]))
  #brcacopy <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$copy!=0,"NodeName"]))
  
  #pull labels
  labeltab <- geneIntTable[geneIntTable$NodeName %in% nodenet,]
  labels <- as.character(labeltab$Gene)
  names(labels) <- labeltab$NodeName
    
  nodenetwork <- subNetwork(nodenet, intome)
  
  nodeRenderInfo(nodenetwork) <- list(shape="circle", 
                                      iheight=.5, iwidth= .5, 
                                      fixedsize=FALSE, label = as.list(labels))
  #inRPPA <- nodenet[nodenet %in% RPPANodes]
  
  #make Node protein a diamond
  RPPAshape <- list()
  RPPAshape[[NodeName]] <- "diamond"  
  #for(RP in inRPPA){RPPAshape[[RP]] <- "diamond"}
  nodeRenderInfo(nodenetwork) <- list(shape=RPPAshape)
  
  if(!is.null(gisticCopyCalls)){
  nodecopyframe <- cellFrame[as.character(cellFrame$NodeName) %in% nodenet,]
  rownames(nodecopyframe) <- nodecopyframe$NodeName
  colnames(nodecopyframe) <- c("Gene", "NodeName", "copyStatus")
  }
  else{
    nodecopyframe <- mutCopyFrame[,c("Gene", "NodeName", "copyStatus")]
  }
  
  nodecolors <- list()
  
  for(nd in rownames(nodecopyframe)){
    ndcol <- "white"
    if(nodecopyframe[nd, "copyStatus"] > 0){ndcol <- "pink"}
    if(nodecopyframe[nd, "copyStatus"] < 0){ndcol <- "lightgreen"}
    nodecolors[[nd]] <- ndcol  
  }
  
  for(nd in brcamut) {nodecolors[[nd]] <- "lightblue"}  
  
  cellRes <- cellResults[[sampleName]]
  
  NodDegree <- cellRes[cellRes$NodeName ==NodeName,"degree"]
  NeighborMuts <- cellRes[cellRes$NodeName == NodeName,"neighborVec"]
  pval <- cellRes[cellRes$NodeName == NodeName, "pvalue"]
  
  nodenetwork <- layoutGraph(nodenetwork)
  if(NeighborMuts > 0 & length(nodes(nodenetwork)) > 1){
    nodeRenderInfo(nodenetwork) <- list(fill = nodecolors)
  if(is.null(fileOut)){
    fileOut <- paste(sampleName,"-",GeneName,".svg", sep="")
  }
  svg(height=7, width=7, filename=fileOut)
  #plot.new()
  #layout(matrix(c(1,2),2,1), heights=c(2,1))
  renderGraph(nodenetwork)
  plotTitle <- paste(NodeName, " (",sampleName," ",", p=", pval, ", n=", NeighborMuts, ", d=", NodDegree, ")", sep="")
  title(plotTitle)
  #hist(distMat[,NodeName], main=plotTitle, xlab="Neighbor Mutations", ylab = "Frequency")
  dev.off()
  }  
}
