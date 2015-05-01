\name{buildGeneIntTable}
\alias{buildGeneIntTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build Mapping Table}
\description{Builds a mapping table between gene symbol and node names in the HPRD interactome network.
}
\usage{
buildGeneIntTable(intome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intome}{
  interactome (HPRD (Human Protein Reference Database) derived from the DLBCL package.)
}
}
\details{ This function takes the HPRD protein protein network that is available in the DLBCL package and produces a mapping table that is needed by functions in the surrogateMutation package.
}
\value{
A mapping table with the following columns:
\item{NodeName}{Node Names in the interactome network. These node names are in the format of GeneSymbol(EntrezID)}
\item{Gene}{Gene Symbols, extracted from the NodeName.}
\item{EntrezID}{Entrez Gene IDs, extracted from the NodeName.}
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
library(DLBCL)
data(interactome)
geneIntTable <- buildGeneIntTable(interactome)
geneIntTable[1:5,]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
