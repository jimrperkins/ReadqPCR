\name{qPCRBatch-class}
\docType{class}
\alias{qPCRBatch-class}
\alias{class:qPCRBatch}
\alias{qPCRBatch}
\alias{qSet}
\alias{exprs,qPCRBatch-method}
\alias{exprs<-,qPCRBatch-method}
\alias{se.exprs,qPCRBatch-method}
\alias{se.exprs<-,qPCRBatch-method}
\alias{exprs.well.order}
\alias{exprs.well.order<-}
\alias{exprs.well.order,qPCRBatch-method}
\alias{exprs.well.order<-,qPCRBatch-method}
\alias{effs}
\alias{effs<-}
\alias{effs,qPCRBatch-method}
\alias{effs<-,qPCRBatch-method}
\alias{se.effs}
\alias{se.effs<-}
\alias{se.effs,qPCRBatch-method}
\alias{se.effs<-,qPCRBatch-method}
\title{
  Class qPCRBatch
}

\description{
  Class to Contain and Describe raw and normalised qPCR Data, as Cq or delta-Cq values. Extends eSet
}
\section{Creating Objects}{
  \code{new("qPCRBatch")}
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{ Object of class \code{AssayData} containing the raw data,
      which will be at minimum a matrix of Cq values. This
      slot can also hold a matrix of well.info values if these are present in the input file read 
      in by \code{read.qPCR} or \code{read.taqman}
    }
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame}
      containing phenotypic data for the samples.
    }
    \item{\code{annotation}}{A character string identifying the
      annotation that may be used for the \code{qPCRBatch}
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
    \item{exprs}{\code{signature(object = "qPCRBatch")}: extracts the Cq
      expression matrix.}
    \item{exprs<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the Cq expression matrix.}
    \item{se.exprs}{\code{signature(object = "qPCRBatch")}: extracts the
      expression matrix with SDs of Cq values.}
    \item{se.exprs<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the expression matrix with SDs of Cq values.}
    \item{exprs.well.order}{\code{signature(object = "qPCRBatch")}: extracts the Cq
      well order matrix (if it exists).}
    \item{exprs.well.order<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the Cq well order matrix.}
    \item{effs}{\code{signature(object = "qPCRBatch")}: extracts the efficiency
      matrix (if it exists).}
    \item{effs<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the efficiency matrix.}
    \item{se.effs}{\code{signature(object = "qPCRBatch")}: extracts the matrix 
      with the standard errors/deviations of the efficiencies (if it exists).}
    \item{se.effs<-}{\code{signature(object = "qPCRBatch", value = "matrix")}:
      replaces the matrix with the standard errors/deviations of the efficiencies.}
  }
}
\note{This class is better described in the vignette.}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet}}
}
\references{
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{James Perkins \email{jimrperkins@gmail.com}}


\keyword{classes}

