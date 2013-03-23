\name{read.LC480SampleInfo}
\encoding{UTF8}
\alias{read.LC480SampleInfo}
\title{ Read sample information data of a qPCR experiment from Light Cycler 480 }
\description{
  Reads sample information data of a qPCR experiment from Light Cycler 480 which is in txt-file and uses the data to populate an object of Class \code{"AnnotatedDataFrame"}.
}
\usage{
read.LC480SampleInfo(file, removeEmptyCols = TRUE,
                     header = TRUE, sep = "\t", quote = "\"", 
                     dec = ".",fill = TRUE, comment.char = "",
                     skip = 0)
}
\arguments{
  \item{file}{the name of the file to be read in.}
  \item{removeEmptyCols}{a logical value which indicates whether the empty column(s) should be removed or not.}
  \item{header}{a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.}
  \item{sep}{the field separator character. Values on each line of the file are separated by this character. If \code{sep = ""} (the default for \code{\link[utils]{read.table}}) the separator is ‘white space’, that is one or more spaces, tabs, newlines or carriage returns.}
  \item{quote}{the set of quoting characters. To disable quoting altogether, use \code{quote = ""}. See \code{\link[base]{scan}} for the behaviour on quotes embedded in quotes. Quoting is only considered for columns read as character, which is all of them unless colClasses is specified.}
  \item{dec}{the character used in the file for decimal points.}
  \item{fill}{logical. If TRUE then in case the rows have unequal length, blank fields are implicitly added. 
  See \code{\link[utils]{read.table}}.}
  \item{comment.char}{character: a character vector of length one containing a single character or an empty string. Use \code{""} to turn off the interpretation of comments altogether.}
  \item{skip}{integer: the number of lines of the data file to skip before beginning to read data.}
}
\details{
  Allows the user to read in sample information data of a qPCR experiment from Light Cycler 480 which is in 
  txt-file. 
  %More details on how to merge the phenodata from sample information and the existing phenodata in the \code{ReadqPCR} package vignette and also on how Cq and amplification efficiency values are calculated in the \code{NormqPCR} package vignette.
}
\value{Object of class \code{"AnnotatedDataFrame"}.}
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

\code{\link{read.LC480}} for reading in the qPCR data from Light Cycler and produce an object of Class \code{"CyclesSet"}.

\code{\link{merge}} for merging the phenodata from \code{read.LC480} and \code{read.LC480SampleInfo} and produce an object of Class \code{"CyclesSet"}.

\code{\link{read.qPCR}} and \code{\link{read.taqman}} for reading in the RT-qPCR data.

%\code{\link[NormqPCR]{CqValues}} for computing Cq values, amplification efficiencies and other important values for quantitative PCR.

%\code{\link[NormqPCR]{Norm}} for technical replication and normalization of the data from this function.

\code{\link[Biobase]{AnnotatedDataFrame-class}}
}
\examples{
  path <- system.file("exData", package = "ReadqPCR")
  LC480.SamInfo <- file.path(path, "LC480_Example_SampleInfo.txt")

  # Read in the sample information data from file "LC480_Example_SampleInfo.txt". 
  samInfo <- read.LC480SampleInfo(LC480.SamInfo)
}
\keyword{classes}