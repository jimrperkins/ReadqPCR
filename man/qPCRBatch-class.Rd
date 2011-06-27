\name{qPCRBatch-class}
\docType{class}
\alias{qPCRBatch-class}
\alias{class:qPCRBatch}
\alias{qPCRBatch}
\alias{qSet}
\alias{exprs.well.order}
\alias{exprs.well.order<-}
\alias{exprs.well.order,qPCRBatch-method}
\alias{exprs.well.order<-,qPCRBatch-method}
\title{
  Class qPCRBatch
}

\description{
  Class to Contain and Describe raw and normalised qPCR Data, as Ct or delta-Ct values. Extends eSet
}
\section{Creating Objects}{
  \code{new("qPCRBatch")}
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{ Object of class \code{AssayData} containing the raw data,
      which will be at minimum a matrix of Ct values. This
      slot can also hold a matrix of well.info values if these are present in the input file read 
      in by \code{read.qPCR} or \code{read.taqman}
    }
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame}
      containing phenotypic data for the samples.
    }
    \item{\code{annotation}}{A character string identifying the
      annotation that may be used for the \code{ExpressionSet}
      instance.
    }
    \item{\code{protocolData}:}{Object of class \code{AnnotatedDataFrame}
      containing protocol data for the samples.
    }
    \item{\code{featureData}}{Object of class \code{AnnotatedDataFrame}
      containing feature-level (e.g., probeset-level) information.
    }
    \item{\code{experimentData}:}{Object of class "MIAME" containing
      experiment-level information.
    } 
    \item{\code{.__classVersion__}:}{Object of class \code{Versions}
      describing the R and Biobase version number used to create the
      instance. Intended for developer use.
    } 
  }
}
\section{Methods}{
  \describe{
    \item{exprs}{\code{signature(object = "qPCRBatch")}: extracts the Ct
      expression matrix.}
    \item{exprs<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the Ct expression matrix.}
    \item{exprs.well.order}{\code{signature(object = "qPCRBatch")}: extracts the Ct
      well order matrix (if it exists).}
    \item{exprs.well.order<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the Ct well order matrix.}
  }
}
\note{This class is better described in the vignette.}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet}}
}
\author{James Perkins}


\keyword{classes}

