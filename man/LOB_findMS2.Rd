\name{LOB_findMS2}
\alias{LOB_findMS2}
\title{
Find MS2 spectra in a rawSpec object based on retention time and precursor m/z
}
\description{
This function allows the user to quickly search a rawSpec object for MS2 cooresponding to either an inputed m/z and retention time or search with multiple MS1 peaks at once by inputing a peaklist.
}
\usage{
LOB_findMS2(rawSpec,data=NULL,mz,rt,rtspan=175,ppm=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawSpec}{
  A \code{rawSpec} object to be searched contianing both MS1 and MS2 spectra.
}
  \item{data}{
  A LOBSTAHS peak list containing info for multiple peaks for which coresponding MS2 spectra will be searched for.
}
  \item{mz}{
  A single numeric equal to the m/z of a peak. 
}
  \item{rt}{
  A single numeric equal to the rt of a peak. 
}
  \item{rtspan}{
  A single numeric equal setting the retention time window around the peak 'rt' (defaults to ±175 seconds).
}
  \item{rt}{
  A single numeric equal setting the mass range around the peak 'm/z' (defaults to ±2.5 ppm). 
}
}
\details{
Using the 'data' argumnet overides the 'mz' and 'rt' arguments. Use 'data' to search with multiple peaks at once. Use 'mz' and 'rt' to search with a single peak.
}
\value{
Returns a \code[data.frame] or a list of \code[data.frame] if multipule peaks were searched with at once. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Henry Holm, hholm@mit.edu
}
\note{
Currently in development.
}
