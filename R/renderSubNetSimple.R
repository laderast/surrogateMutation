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

summarizeNetworks <- function(IDvec, intome, mutCopyFrames, surrogateTable, geneIntTable, fileOut="www/netgraph.svg", filterNodes=NULL){
  #library(gplots)
  if(length(IDvec) < 1){return(NULL)}
  
  nodeTable <- surrogateTable[surrogateTable$ID %in% IDvec,]
  
  fullNet <- NULL
  surrNodes <- unique(as.character(nodeTable$NodeName))
  #for each row of the table, grab neighbors of the surrogate node
  for(i in 1:nrow(nodeTable)){
    Sample <- nodeTable[i,"Sample"]
    NodeName <- as.character(nodeTable[i,"NodeName"])
    nodenet <- c(NodeName,intersect(as.character(inEdges(NodeName, intome)[[1]]), 
                                    as.character(mutCopyFrames[[Sample]]$NodeName)))
    fullNet <- c(fullNet, nodenet)
  }
  
  #count the number of times a node is mutated
  fullNet <- table(fullNet)
  
  if(!is.null(filterNodes)){
    fullNet <- fullNet[fullNet > filterNodes]
  }
  
  surrNodeCount <- table(as.character(nodeTable$NodeName)) +1
  
  #subtract 1 from the count for the surrogate nodes, since we already count them
  fullNet[names(surrNodeCount)] <- fullNet[names(surrNodeCount)] - surrNodeCount 
  
  #don't graph if the network only has 1 node or less  
  if(sum(fullNet) < 2 | length(fullNet) < 2){return(NULL)}
  
  numColors <- max(fullNet) + 1
  #colorPanel <- topo.colors(numColors)
  colorPanel <- colorpanel(n=numColors, low="white", high="steelblue")
  
  #assign colors according to degres
  nodeColorList <- list()
  for(nd in names(fullNet)){
    nodeColorList[[nd]] <- colorPanel[fullNet[nd]]
  }
  
  #extract the relevant subgraph from the interactome
  fullGraph <- subGraph(names(fullNet), intome)
  
  labeltab <- geneIntTable[geneIntTable$NodeName %in% names(fullNet),]
  labels <- as.character(labeltab$Gene)
  names(labels) <- labeltab$NodeName
  
  #print(labels)
  
  #default node render properties
  nodeRenderInfo(fullGraph) <- list(shape="circle", 
                                    iheight=.5, iwidth= .5, 
                                    fixedsize=FALSE, label = as.list(labels))
  #inRPPA <- nodenet[nodenet %in% RPPANodes]
  
  #make surrogate node proteins a diamond
  shapeList <- list()
  for(surr in surrNodes){
    shapeList[[surr]] <- "diamond"  
  }
  #for(RP in inRPPA){RPPAshape[[RP]] <- "diamond"}
  nodeRenderInfo(fullGraph) <- list(shape=shapeList)
  
  #color each node by degree of mutation
  fullGraph <- layoutGraph(fullGraph, layoutType="fdp")
  
  nodeRenderInfo(fullGraph) <- list(fill = nodeColorList)
  if(!is.null(fileOut)){
    svg(height=7, width=7, filename = fileOut)
    renderGraph(fullGraph)
    #title("Summary Network")
    dev.off()
  }
  else{renderGraph(fullGraph)
  }
  return(fullNet)
}

getSubnetworkForSamples <- function(surrogateNode, intome, surrResult, sampList=NULL){
  require(graph)
  
  intomeNodes <- inEdges(surrogateNode, intome)[[1]]
  intomeNodes <- c(surrogateNode, intomeNodes)
  
  mutCopyFrames <- surrResult$mutCopyFrames
  
  if(!is.null(sampList)){
    mutCopyFrames <- mutCopyFrames[[sampList]]
  }
  
  allMutList <- lapply(mutCopyFrames, function(x){as.character(x$NodeName)})
  
  allMuts <- Reduce(union, allMutList)
  #print(allMuts)
  
  allMutCount <- rep(0, length(allMuts))
  names(allMutCount) <- allMuts
  #need function for counting mutations here
  countMuts <- lapply(allMutList, function(x){allMutCount[x] <- allMutCount[x] +1 
                                              return(allMutCount)})
  allMutCount <- Reduce("+", countMuts)
  
  print(allMutCount)
  
  inBoth <- intersect(intomeNodes,allMuts)
  if(!surrogateNode %in% inBoth){
  inBoth <- c(inBoth, surrogateNode)}
  allMutCount <- allMutCount[inBoth]
  
  subNet <- subGraph(inBoth, intome)
  nodeDataDefaults(subNet, "counts") <- 0 
    
  nodeData(subNet, attr="counts") <- allMutCount
  
  nodeNames <- unlist(lapply(nodes(subNet), function(x){strsplit(x, "\\(")[[1]][1]}))
  nodes(subNet) <- nodeNames
  
  return(subNet)
}