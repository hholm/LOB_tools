\name{LOB_findMS2}
\alias{LOB_findMS2}
\title{
Find MS2 spectra in a rawSpec object based on retention time and precursor m/z
}
\description{
This function allows the user to quickly search a rawSpec object for MS2 cooresponding to either an inputed m/z and retention time or search with multiple MS1 peaks at once by inputing a peaklist.
}
\usage{
LOB_findMS2(XCMSnExp, peakdata = NULL, mz, rt, rtspan = 175, ppm_pre = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XCMSnExp}{
  A \code{XCMSnExp} object to be searched contianing both MS1 and MS2 spectra.
}
  \item{peakdata}{
  A \code{data.frame} generated from LOBSTAHS peaklist containing info for multiple peaks for which coresponding MS2 spectra will be searched. The input 'peakdata' must have columns 'peakgroup_rt', 'LOBdbase_mz', and 'compound_name' with rows for each lipid.
}
  \item{mz}{
  A single numeric vector equal to the m/z of a peak.
}
  \item{rt}{
  A single numeric vector equal to the rt of a peak.
}
  \item{rtspan}{
  A single numeric vector equal setting the retention time window around the peak 'rt' (defaults to ±175 seconds).
}
  \item{ppm}{
  A single numeric vector equal setting the mass range around the peak 'm/z' for which to search for scans (defaults to ±100 ppm).
}
}
\details{
Using the 'peakdata' argumnet overides the 'mz' and 'rt' arguments. Use 'peakdata' to search using multiple annotations at once. Use 'mz' and 'rt' to search with a single peak.
}
\value{
Returns a \code{data.frame} or a list of \code{data.frame} if multipule peaks were searched with at once.
}

\examples{
#Example - find MS2 of DNPPE
DNPPE_ms2 <- LOB_findMS2(XCMSnExp = rawSpec,
                         mz = 875.550487,
                         rt = 1000,
                         ppm = 2.5)
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Henry Holm, hholm@mit.edu \cr
Dan Lowenstien, dlowenstein@whoi.edu
Max Jahns, mjahns@mit.edu
Shavonna Bent, sbent@mit.edu
}
\note{
Currently in development.
}
