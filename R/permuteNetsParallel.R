permuteNetsParallel <-
function(sampleSize, numPermutes=10000, intAdjmat=intAdjmat, cellLine="",
                                cores = 4, prefix="", fileout=paste(prefix,cellLine, "permute-test-", 
                                                         sampleSize, ".txt", sep="")){
  
  library(utils)
  library(foreach)
  #library(iterators)
  library(doMC)
  registerDoMC(cores=cores)
  
  perms <- iter(1:numPermutes)
  
  #permutation script parallelized and using adjacency matrix
  #intPermute <- intome
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
  


  
  
  intNames <- rownames(intAdjmat)
  
  #distribution1res <- matrix(ncol=length(RPPANodes), nrow=numPermutes, data=NA) 
  
  distribution1res <- foreach(i = 1:numPermutes, .combine='rbind') %dopar% {
  #for(i in 1:numPermutes){
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

    #cntMuts is function for apply to look for intersection in adjacency matrix and mutations
    cntMuts <- function(x){
      ed <- which(x>0)
      idx <- ed[ed %in% names(freqvec)]
      #print(idx)
      mutcount <- sum(freqvec[as.character(idx)])
      mutcount
    }
    
    cntMutEnv <- new.env()
    assign("freqvec", freqvec, cntMutEnv)
    environment(cntMuts) <- cntMutEnv
    
    #count number of mutations
    mutcounts <- apply(intAdjmat, 2, cntMuts)
    
    mutcounts
    
    #distribution1res[i,] <- mutcounts
    
  }

  
  #colnames(distribution1res) <- RPPANodes
  write.table(distribution1res, fileout, quote=F, sep="\t", row.names=F)
  distribution1res
  
}
