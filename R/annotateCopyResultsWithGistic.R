#' Given Copy Number Results, annotate with GISTIC files
#'
#' @param mutCopyFrames 
#' @param gisticNodeMapped 
#'
#' @return
#' @export
#'
#' @examples
annotateCopyResultsWithGistic <-
function(mutCopyFrames, gisticNodeMapped){
  for(sampleName in names(mutCopyFrames)){
    gisticSampleFrame <- gisticNodeMapped[,c("Gene.Symbol", "NodeName", sampleName)]
    colnames(gisticSampleFrame)[3] <- "Sample"
    rownames(gisticSampleFrame) <- gisticSampleFrame$NodeName
    mCF <- mutCopyFrames[[sampleName]]
    #alterNodes <- mCF[mCF$copy==1,"NodeName"]
    gSFCopy <- gisticSampleFrame[gisticSampleFrame$NodeName %in% mCF$NodeName,]
    mGVec <- rep(0, nrow(mCF))
    names(mGVec) <- mCF$NodeName
    for(nam in mCF$NodeName){
      if(nam %in% gisticSampleFrame$NodeName){
        if(gisticSampleFrame[nam,"Sample"] > 0) {mGVec[nam] <- 1}
        if(gisticSampleFrame[nam, "Sample"] < 0) {mGVec[nam] <- -1}
      }
    }
    mCF <- data.frame(mCF, copyStatus=mGVec)
    mutCopyFrames[[sampleName]] <- mCF
  }
  return(mutCopyFrames)
}
