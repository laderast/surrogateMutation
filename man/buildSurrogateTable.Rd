\name{buildSurrogateTable}
\alias{buildSurrogateTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build melted surrogate table from surrogateResults.
}
\description{
Given a surrogateResults object from calcNeighborMutations(), generated a 'melted' table that summarizes all surrogateMutation/sample pairs.
}
\usage{
buildSurrogateTable(surrogateResult)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surrogateResult}{
  A result object from calcNeighborMutations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A melted data frame with the following columns:

\item{ID}{The ID for the row, which is the SampleName and Surrogate gene pasted together. For example, for the sample TCGA.A1.A0SD and the gene ABL1 pair, the id would "TCGA.A1.A0SD-ABL1"}
\item{Gene}{The Gene Symbol for the surrogate gene.}
\item{Sample}{The Sample name. Can be patient or cell line id.}
\item{NodeName}{Node Identifier of Surrogate Gene in HPRD network.}
\item{neighbor}{Number of neighboring mutations observed in the Sample for the surrogate gene}
\item{degree}{Total number of connections in the interaction network for the surrogate gene.}
\item{pvalue}{P-value from permutation analysis. Note this p-value is currently unadjusted for multiple comparisions.}
\item{isMutated}{Is the surrogate gene itself mutated?}

This surrogate table is used as one of the inputs to the Shiny interactive heatmap.  For an example of the interactive heatmap, please visit: https://tladeras.shinyapps.io/surrogateShiny/

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
##load surrogate result object
data(surrogateResult)
surrogateTable <- buildSurrogateTable(surrogateResult)

##filter the surrogate table based on p-value and return surrogate mutations greater than 5.
filteredSurrTable <- surrogateTable[surrogateTable$pvalue < 0.05 & surrogateTable$neighbor > 3,]
filteredSurrTable

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
