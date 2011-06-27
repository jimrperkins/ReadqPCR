\name{read.qPCR}
\alias{read.qPCR}
\alias{read_qPCRBatch}
\alias{checkValidqPCRFilename}
\title{ Read user formatted qPCR data and produce a qPCRBatch }
\description{
  Reads RT-qPCR data in format specified in the ReadqPCR vignette and uses the data to populate an object of class
  \code{"qPCRBatch"}.
}
\usage{
read.qPCR(filename = character(0),
          phenoData = new("AnnotatedDataFrame"),
          notes = "",
          verbose = FALSE)
}
\arguments{
  \item{filename}{file name (must be formatted as shown in vignette).}
  \item{phenoData}{an \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
     object, a \code{character} of length one, or a \code{data.frame}.}
  \item{notes}{notes.}
  \item{verbose}{verbosity flag. If true more messages are given to the user on the processing steps}
}
\details{
  Permits the user to read in qPCR Ct value data in a predefined format 
  (more details on this format in the \code{ReadqPCR} package vignette), alongside 
  phenotypic data and further notes about the data. If \code{phenoData} is a \code{data.frame}, 
  it is converted to an \code{AnnotatedDataFrame}. If it is \code{NULL} then a default object
  of class \code{AnnotatedDataFrame} is created, whose \code{pData} is a
  \code{data.frame} with rownames being the names of the samples, and
  with one column \code{sample} with an integer index. 
  More details on how technical replicates are handled in the \code{ReadqPCR} package vignette
}
\value{Object of class \code{"qPCRBatch"}.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\note{
%}
\seealso{ \code{\link[Biobase]{ExpressionSet-class}} }
\examples{
  path <- system.file("exData", package = "ReadqPCR")
  qPCR.example <- file.path(path, "qPCR.example.txt")
  qPCRBatch.qPCR <- read.qPCR(qPCR.example)
}
\keyword{classes}
