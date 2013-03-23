\name{merge}
\alias{merge}
\alias{merge-methods}
\alias{merge,eSet,AnnotatedDataFrame-method}
\alias{merge,AnnotatedDataFrame,eSet-method}

\title{Generic Function for Merging eSet and AnnotatedDataFrame}
\description{
  Generic function for the merging an \code{eSet} and an \code{AnnotatedDataFrame};
  that is, the informatiton given in the \code{AnnotatedDataFrame} is merged into
  the appropriate slot of the given \code{eSet}. This can be slot \code{phenoData}
  (default), \code{featureData}, or \code{protocolData}.
}
\usage{
merge(x, y, ...)
\S4method{merge}{eSet,AnnotatedDataFrame}(x, y, eSet.slot = "phenoData", 
                   by = intersect(names(pData(x)), names(pData(y))),
                   by.x = by, by.y = by, 
                   all = FALSE, all.x = all, all.y = all,
                   sort = FALSE, suffixes = c(".x",".y"), 
                   incomparables = NULL, ...)
\S4method{merge}{AnnotatedDataFrame,eSet}(x, y, eSet.slot = "phenoData", 
                   by = intersect(names(pData(x)), names(pData(y))),
                   by.x = by, by.y = by, 
                   all = FALSE, all.x = all, all.y = all,
                   sort = FALSE, suffixes = c(".x",".y"), 
                   incomparables = NULL, ...)
}
\arguments{
  \item{x}{ object of class \code{eSet} or \code{AnnotatedDataFrame} }
  \item{y}{ object of class \code{eSet} or \code{AnnotatedDataFrame} }
  \item{eSet.slot}{ name of the slot of the given \code{eSet} which is
  merged with the given \code{AnnotatedDataFrame}. This can be 
  \code{"phenoData"} (default), \code{"featureData"}, or \code{"protocolData"}}
  \item{by}{specifications of the columns used for merging.}
  \item{by.x}{specifications of the columns used for merging.}
  \item{by.y}{specifications of the columns used for merging.}
  \item{all}{logical; \code{all = L} is shorthand for \code{all.x = L} and \code{all.y = L}, where L is either \code{TRUE} or \code{FALSE}.}
  \item{all.x}{logical; if \code{TRUE}, then extra rows will be added to the output, one for each row in x that has no matching row in y. These rows will have NAs in those columns that are usually filled with values from y. The default is \code{FALSE}, so that only rows with data from both x and y are included in the output.}
  \item{all.y}{logical; analogous to \code{all.x}.}
  \item{sort}{logical. Should the result be sorted on the by columns?}
  \item{suffixes}{a character vector of length 2 specifying the suffixes to be used for making unique the names of columns in the result which not used for merging (appearing in by etc).}
  \item{incomparables}{values which cannot be matched. See \code{\link[base]{match}}.}
  \item{\dots}{ additional arguments to be passed to or from methods. }
}
\details{For details on the arguments see \code{\link[base]{merge}}. }
\value{ Object of class \code{"CyclesSet"}. }
\section{Methods}{
\describe{
  \item{x = "eSet", y = "AnnotatedDataFrame":}{ merge \code{y} into specified 
  slot of \code{x}. }

  \item{x = "AnnotatedDataFrame", y = "eSet":}{ merge \code{x} into specified
  slot of \code{y}. }
}}
%\references{ ~put references to the literature/web site here ~ }
\author{Nor Izayu Abdul Rahman, Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{ ~~further notes~~ }
\seealso{\code{\link[base]{merge}}, \code{\link[Biobase:class.eSet]{eSet}}}
\examples{
  path <- system.file("exData", package = "ReadqPCR")
  
  LC480.example <- file.path(path, "LC480_Example.txt")
  ## Read in the raw qPCR data from file "LC480_Example.txt"
  ## with maximum cycle to be read in the values is 45 (default).
  cycData <- read.LC480(file = LC480.example)

  LC480.SamInfo <- file.path(path, "LC480_Example_SampleInfo.txt")
  # Read in the sample information data from file "LC480_Example_SampleInfo.txt". 
  samInfo <- read.LC480SampleInfo(LC480.SamInfo)

  cycData1 <- merge(cycData, samInfo)
}
\concept{functional}
\concept{truncated moment}
\concept{moment}
\keyword{distribution}
\keyword{methods}
