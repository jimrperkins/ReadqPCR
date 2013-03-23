\name{qpcR2CyclesSet}
\alias{qpcR2CyclesSet}
\title{ Transform a dataset from qpcR Package into CyclesSet}
\description{
  Transform a dataset from package \pkg{qpcR} into an object of class \code{"CyclesSet"}.
}
\usage{
qpcR2CyclesSet(x, cyc = 1, cycleThreshold)
}
\arguments{
  \item{x}{the name of the dataset from package \pkg{qpcR} to be transformed.}
  \item{cyc}{the column number which contains cycle data.}
  \item{cycleThreshold}{maximum number of cycles which the data will be counted and transformed.}
}
\details{
  Allows the user to transform a dataset from package \pkg{qpcR} into an object of class 
  \code{"CyclesSet"} class, alongside phenotypic data. 
  %More details on how Cq and amplification efficiency values are calculated in the \code{NormqPCR} package vignette.
}
\value{Object of class \code{"CyclesSet"}.}
\references{
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ Nor Izayu Abdul Rahman, Matthias Kohl \email{Matthias.Kohl@stamats.de} }
%\note{
%}
\seealso{ \code{\link[qpcR]{reps}}, \code{\link{CyclesSet-class}} }
\examples{
  library(qpcR)

  ## Transform the reps dataset from qpcR Package 
  ## with maximum cycle to be read in the value is 45.
  repsdata <- qpcR2CyclesSet(reps, cyc=1, cycleThreshold=45)
  
  ## Transform the batsch1 dataset from qpcR Package 
  ## with maximum cycle to be read in the value is 40.
  batsch1data <- qpcR2CyclesSet(batsch1, cyc=1, cycleThreshold=40)
}
\keyword{classes}