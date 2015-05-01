buildSurrogateTable <-
function(surrogateResult){
  cellRes <- surrogateResult$cellResults
  names(cellRes) <- make.names(names(cellRes))
  cellRes2 <- lapply(names(cellRes), function(x){ tab <- cellRes[[x]]
                                                  ID <- paste(x,as.character(tab$Gene), sep="-")
                                                  Sample <- rep(x, nrow(tab))
                                                  out <- data.frame(ID, Gene=tab$Gene, Sample, NodeName = tab$NodeName,neighbor=tab$neighborVec,
                                                                    degree=tab$degree, pvalue=tab$pvalue, isMutated=tab$isMutated)
                                                  rownames(out) <- ID
                                                  return(out)
  })
  cellResTable <- do.call(rbind, cellRes2)
  return(cellResTable)
}
