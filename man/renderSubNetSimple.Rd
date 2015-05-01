\name{renderSubNetSimple}
\alias{renderSubNetSimple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
renderSubNetSimple
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
renderSubNetSimple(NodeName, sampleName, GeneName, intome = intome, gisticCopyCalls = NULL, resultObj, fileOut = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NodeName}{
%%     ~~Describe \code{NodeName} here~~
}
  \item{sampleName}{
%%     ~~Describe \code{sampleName} here~~
}
  \item{GeneName}{
%%     ~~Describe \code{GeneName} here~~
}
  \item{intome}{
%%     ~~Describe \code{intome} here~~
}
  \item{gisticCopyCalls}{
%%     ~~Describe \code{gisticCopyCalls} here~~
}
  \item{resultObj}{
%%     ~~Describe \code{resultObj} here~~
}
  \item{fileOut}{
%%     ~~Describe \code{fileOut} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (NodeName, sampleName, GeneName, intome = intome, gisticCopyCalls = NULL, 
    resultObj, fileOut = NULL) 
{
    library(Rgraphviz)
    library(BioNet)
    if (!is.null(gisticCopyCalls)) {
        cellFrame <- gisticCopyCalls[, c("Gene.Symbol", "NodeName", 
            sampleName)]
    }
    mutCopyFrame <- resultObj$mutCopyFrames[[sampleName]]
    cellResults <- resultObj$cellResults
    nodenet <- c(NodeName, intersect(as.character(inEdges(NodeName, 
        intome)[[1]]), as.character(mutCopyFrame$NodeName)))
    brcamut <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$mutations != 
        0, "NodeName"]))
    labeltab <- geneIntTable[geneIntTable$NodeName \%in\% nodenet, 
        ]
    labels <- as.character(labeltab$Gene)
    names(labels) <- labeltab$NodeName
    nodenetwork <- subNetwork(nodenet, intome)
    nodeRenderInfo(nodenetwork) <- list(shape = "circle", iheight = 0.5, 
        iwidth = 0.5, fixedsize = FALSE, label = as.list(labels))
    RPPAshape <- list()
    RPPAshape[[NodeName]] <- "diamond"
    nodeRenderInfo(nodenetwork) <- list(shape = RPPAshape)
    if (!is.null(gisticCopyCalls)) {
        nodecopyframe <- cellFrame[as.character(cellFrame$NodeName) \%in\% 
            nodenet, ]
        rownames(nodecopyframe) <- nodecopyframe$NodeName
        colnames(nodecopyframe) <- c("GeneName", "NodeName", 
            "copyStatus")
    }
    else {
        nodecopyframe <- mutCopyFrame[, c("GeneName", "NodeName", 
            "copyStatus")]
    }
    nodecolors <- list()
    for (nd in rownames(nodecopyframe)) {
        ndcol <- "white"
        if (nodecopyframe[nd, "copyStatus"] > 0) {
            ndcol <- "pink"
        }
        if (nodecopyframe[nd, "copyStatus"] < 0) {
            ndcol <- "lightgreen"
        }
        nodecolors[[nd]] <- ndcol
    }
    for (nd in brcamut) {
        nodecolors[[nd]] <- "lightblue"
    }
    cellRes <- cellResults[[sampleName]]
    NodDegree <- cellRes[cellRes$NodeName == NodeName, "degree"]
    NeighborMuts <- cellRes[cellRes$NodeName == NodeName, "neighborVec"]
    pval <- cellRes[cellRes$NodeName == NodeName, "pvalue"]
    nodenetwork <- layoutGraph(nodenetwork)
    if (NeighborMuts > 0 & length(nodes(nodenetwork)) > 1) {
        nodeRenderInfo(nodenetwork) <- list(fill = nodecolors)
        if (is.null(fileOut)) {
            fileOut <- paste(sampleName, "-", GeneName, ".svg", 
                sep = "")
        }
        svg(height = 7, width = 7, filename = fileOut)
        renderGraph(nodenetwork)
        plotTitle <- paste(NodeName, " (", sampleName, " ", ", p=", 
            pval, ", n=", NeighborMuts, ", d=", NodDegree, ")", 
            sep = "")
        title(plotTitle)
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
