\name{read.LC480}
\encoding{UTF8}
\alias{read.LC480}
\title{ Read in raw qPCR data of Light Cycler 480}
\description{
  Reads in raw qPCR data of Light Cycler 480 and uses the data to populate 
  an object of class \code{"CyclesSet"}.
}
\usage{
read.LC480(file, colNames = c("Sample position", "Sample name",
                              "Program number", "Segment number", 
                              "Cycle number", "Acquisition time", 
                              "Acquisition temperature", 
                              "Fluorescence data"),
           cycleThreshold = 45, fileType = "txt", skip = 1,
           header = TRUE, sep = "\t", quote = "\"", dec = ".",
           fill = TRUE, comment.char = "")
}
\arguments{
  \item{file}{the name of the file to be read in.}
  \item{colNames}{a character vector of names to be assumed for the columns.}
  \item{cycleThreshold}{maximum number of cycles which will be read in.}
  \item{fileType}{the type of the file.}
  \item{skip}{integer: the number of lines of the data file to skip before beginning to read data.}
  \item{header}{a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.}
  \item{sep}{the field separator character. Values on each line of the file are separated by this character. If \code{sep = ""} (the default for \code{\link[utils]{read.table}}) the separator is ‘white space’, that is one or more spaces, tabs, newlines or carriage returns.}
  \item{quote}{the set of quoting characters. To disable quoting altogether, use \code{quote = ""}. See \code{\link[base]{scan}} for the behaviour on quotes embedded in quotes. Quoting is only considered for columns read as character, which is all of them unless colClasses is specified.}
  \item{dec}{the character used in the file for decimal points.}
  \item{fill}{logical. If TRUE then in case the rows have unequal length, blank fields are implicitly added. 
  See \code{\link[utils]{read.table}}.}
  \item{comment.char}{character: a character vector of length one containing a single character or an empty string. Use \code{""} to turn off the interpretation of comments altogether.}
}
\details{
  Allows the user to read in qPCR fluorescence data from Light Cycler 480 which has been exported to a txt-file, 
  alongside phenotypic data. 
%  More details on how to read in sample information data from Light Cycler in the \code{ReadqPCR} package vignette and also on how Cq and amplification efficiency values are calculated in the \code{NormqPCR} package vignette.
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
\seealso{ 
\code{\link[utils]{read.table}}

\code{\link{read.LC480SampleInfo}} for reading in sample information of qPCR data from Light Cycler and produce an object of Class \code{"CyclesSet"}.

\code{\link{merge}} for merging the phenodata from \code{read.LC480} and \code{read.LC480SampleInfo} and produce an object of Class \code{"CyclesSet"}.

%\code{\link[NormqPCR]{CqValues}} for computing Cq values, amplification efficiencies and other important values for quantitative PCR.

%\code{\link[NormqPCR]{Norm}} for technical replication and normalization of the data from this function.

\code{\link{read.qPCR}} and \code{\link{read.taqman}} for reading in the RT-qPCR data consisting of Cq values.

\code{\link[Biobase]{eSet}}, \code{\link{CyclesSet-class}}

}
\examples{
  path <- system.file("exData", package = "ReadqPCR")
  LC480.example <- file.path(path, "LC480_Example.txt")

  ## Read in the raw qPCR data from file "LC480_Example.txt"
  ## with maximum cycle to be read in the values is 45 (default).
  cycData <- read.LC480(file = LC480.example)
  
  ## Read in the data from file "LC480_Example.txt" 
  ## with maximum cycle to be read in the values is 50.
  rawdata <- read.LC480(file=LC480.example, cycleThreshold=50)
}
\keyword{classes}
