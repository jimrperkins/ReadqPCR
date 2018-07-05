\name{ReadqPCR-package}
\alias{ReadqPCR-package}
\alias{ReadqPCR}
\docType{package}
\title{Read qPCR data}
\description{
The package provides functions to read raw RT-qPCR data of different platforms.
}
\details{
\tabular{ll}{
Package: \tab ReadqPCR\cr
Type: \tab Package\cr
Version: \tab 1.5.3\cr
Date: \tab 2013-03-23\cr
Depends: \tab R(>= 2.14.0), Biobase, methods, affy\cr
Imports: \tab Biobase\cr
Suggests: \tab qpcR\cr
License: \tab LGPL-3\cr
LazyLoad: \tab yes\cr
}
library(ReadqPCR)
}
\author{
James Perkins <jimrperkins@gmail.com>\cr
Matthias Kohl <Matthias.Kohl@stamats.de>\cr
Nor Izayu Abdul Rahman

Maintainer: James Perkins <jimrperkins@gmail.com>
}
\references{
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## some examples are given in the vignette
\dontrun{
library(ReadqPCR)
vignette("ReadqPCR")
}
}
\keyword{ package }
