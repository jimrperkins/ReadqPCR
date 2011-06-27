\name{read.taqman}
\alias{read.taqman}
\alias{read_TaqBatch}
\alias{checkValidTaqmanFilenames}
\title{ Read Taqman qPCR data and produce a qPCRBatch }
\description{
  Reads Taqman RT-qPCR data and uses the data to populate an object of class
  \code{"qPCRBatch"}.
}
\usage{
read.taqman(\dots,
           filenames = character(0),
           phenoData = new("AnnotatedDataFrame"),
           notes = "", 
           verbose = FALSE)
}
\arguments{
  \item{\dots}{file names separated by comma.}
  \item{filenames}{file names in a character vector.}
  \item{phenoData}{an
  \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
  object, a \code{character} of length one, or a \code{data.frame}.}
  \item{notes}{notes.}
  \item{verbose}{verbosity flag. If true more messages are given to the user on the processing steps}
}
\details{
  Permits the user to read in qPCR Ct value data from an sds output file, 
  alongside phenotypic data and further notes about the data. If \code{phenoData} 
  is a \code{data.frame}, it is converted to an \code{AnnotatedDataFrame}. If it
  is \code{NULL} then a default object  of class \code{AnnotatedDataFrame} is
  created, whose \code{pData} is a \code{data.frame} with rownames being the
  names of the , and with one column \code{sample} with an integer index.
  More details on how technical replicates are handled in the \code{ReadqPCR}
package vignette
}
\value{Object of class \code{"qPCRBatch"}.}
%\references{ ~put references to the literature/web site here ~ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk} }
%\note{
%}
\seealso{ \code{\link[Biobase]{ExpressionSet-class}} }
\examples{
  path <- system.file("exData", package = "ReadqPCR")
  # read 1 file
  taqman.example <- file.path(path, "example.txt")
  qPCRBatch.taq <- read.taqman(taqman.example)
  # read 2 files
  taqman.example.second.file <- file.path(path, "example2.txt")
  qPCRBatch.taq.two.files <- read.taqman(taqman.example, 
                             taqman.example.second.file)
}
\keyword{classes}
