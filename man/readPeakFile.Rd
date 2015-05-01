\name{readPeakFile}
\alias{readPeakFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readPeakFile
}
\description{
Given a gistic peak file (amp_genes or del_genes) that contains columns of genes organized by gene peak, extract the gene lists from all peaks as a single vector.
}
\usage{
readPeakFile(gisticPeakFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gisticPeakFile}{
  A GISTIC peak file supplied from running GISTIC.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character vector with all genes extracted from each peak column.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ted Laderas
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
Please refer to the vignette for an example.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
