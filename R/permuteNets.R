permuteNets <- function(sampleSize, numPermutes=10000, intAdjmat=intAdjmat, cellLine="", prefix="",
                        fileout=paste(prefix,cellLine, "permute-test-", sampleSize, ".txt", sep="")){
  
  #permutation script using adjacency matrix
  intPermute <- intome
  #RPPANodes <- unique(RPPANodes)
  
  #numPermutes <- 10000
  #sampleSize is number of mutations in network
  #sampleSize <- 368
  
  #intAdjmat <- as(intome, "graphAM")
  #intAdjmat <- intAdjmat@adjMat
  #rownames(intAdjmat) <- colnames(intAdjmat)
  #intAdjmat <- as.data.frame(intAdjmat)
  #intAdjmat <- intAdjmat[,RPPANodes]
  
  
  #distribution 1 probability is uniform distribution
  distribution1 <- rep(1,length(rownames(intAdjmat)))
  distribution1 <- distribution1/length(rownames(intAdjmat))
  names(distribution1) <- rownames(intAdjmat)
  
  #print(distribution1)
  
  #cntMuts is function for apply to look for intersection in adjacency matrix and mutations
  cntMuts <- function(x){
    ed <- which(x>0)
    idx <- ed[ed %in% names(freqvec)]
    #print(idx)
    mutcount <- sum(freqvec[as.character(idx)])
    mutcount
  }
  
  intNames <- rownames(intAdjmat)
  
  distribution1res <- matrix(ncol=ncol(intAdjmat), nrow=numPermutes, data=NA)
  #print("prename")
  colnames(distribution1res) <- colnames(intAdjmat)
  #print("stop here")
  
  for(i in 1:numPermutes){
    #print(i)
    #sample with replacement for mutations
    samp <- sample(names(distribution1), size=sampleSize, prob = distribution1, replace=TRUE)
    samp <- table(samp)
    sampframe <- data.frame(row.names(samp),samp)
    colnames(sampframe) <- c("NodeNam", "NodeName", "Freq")
    
    #index in IntNames of all nodes with mutations 
    inds <- which(intNames %in% names(samp))  
    #pull counts with mutation index
    indframe <- data.frame(inds, NodeName = intNames[inds])
    
    #merge frame with mutation index to build vector of counts with index numbers as names
    freqframe <- merge(sampframe, indframe, by.x = "NodeName", by.y="NodeName")
    freqvec <- freqframe$Freq
    names(freqvec) <- freqframe$inds
    
    cntMutEnv <- new.env()
    assign("freqvec", freqvec, cntMutEnv)
    environment(cntMuts) <- cntMutEnv
    
    #count number of mutations
    mutcounts <- apply(intAdjmat, 2, cntMuts)
    
    #print(mutcounts)
    #inds <- which(names(distribution1) %in% samp)
    
    distribution1res[i,] <- mutcounts
    
    #print(i)
    
    if(i %% 1000 == 0){
      write.table(distribution1res, fileout, quote=F, sep="\t")
      print(paste(cellLine, i))
    }
  }
  distribution1res
  
}
