\name{CyclesSet-class}
\docType{class}
\alias{CyclesSet-class}
\alias{class:CyclesSet}
\alias{CyclesSet}
\alias{exprs,CyclesSet-method}
\alias{exprs<-,CyclesSet-method}
\title{
  Class CyclesSet
}

\description{
  Class to contain and describe raw fluorescence data. Extends eSet
}
\section{Creating Objects}{
  \code{new("CyclesSet")}
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{ Object of class \code{AssayData} containing the raw data,
      which will be a matrix of fluorescence values.
    }
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame}
      containing phenotypic data for the samples.
    }
    \item{\code{annotation}}{A character string identifying the
      annotation that may be used for the \code{CyclesSet}
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
    \item{exprs}{\code{signature(object = "CyclesSet")}: extracts the matrix with 
      the fluorescence values.}
    \item{exprs<-}{\code{signature(object = "CyclesSet", value = "matrix")}:
      replaces the matrix with the fluorescence values.}
  }
}
%\note{This class is better described in the vignette.}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet}}
}
\references{
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{Nor Izayu Abdul Rahman, Matthias Kohl \email{Matthias.Kohl@stamats.de}}


\keyword{classes}

