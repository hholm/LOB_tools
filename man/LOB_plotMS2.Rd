\name{LOB_plotMS2}
\alias{LOB_plotMS2}
\title{
Find and plot MS2 spectra in a rawSpec object using lipid annotations.
}
\description{
This function allows the user to quickly search and plot a MS2 spectra from a rawSpec object cooresponding to either an inputed m/z and retention time or search with multiple MS1 peaks at once by inputing a peaklist.
}
\usage{
LOB_plotMS2(XCMSnExp, peakdata = NULL, plot_file = "closest_scan", mz = NULL, rt = NULL, rtspan = 175, ppm_pre = 100, ppm = 2.5, window = 1, diagnostic = NULL, diagnostic_ppm = 15, NL = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XCMSnExp}{
  A \code{XCMSnExp} object to be searched contianing both MS1 and MS2 spectra.
}
  \item{peakdata}{
  A LOBSTAHS peaklist \code{data.frame} containing info for multiple peaks for which coresponding MS2 spectra will be searched. The input 'peakdata' must have columns 'peakgroup_rt', 'LOBdbase_mz', and 'compound_name' with rows for each lipid.
}
  \item{plot_file}{
  A charactor vector reading either 'most_scans','highest_int', or 'closest_scan' controlling which file is chosen for plotting (defaults to 'closest_scans'). To use a specific file, a character string exactly matching a sample name from the XCMSnExp object can be used (use \code{sampleNames(XCMSnExp)} to see files included in the object). Value 'most_scans' plots spectra from the file with the most MS2 scans found, 'highest_int' selects the file that contains the scan with the highest precursor intensity, 'closest_scan' plots spectra from the file with a scan closest to the specified rt.
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
  \item{ppm_pre}{
  A single numeric equal setting the mass range around the peak 'm/z' for which to search for scans (defaults to ±100 ppm).
}
  \item{ppm}{
  A single numeric equal setting the mass tolerence for the plotted chromatogram (defaults to ±2.5 ppm).
}
  \item{window}{
  A single numeric equal setting the zoom of the MS1 chromatogram window plotted with 1 equalling 100% of the rt region searched (defaults to '1', i.e. 100%). For example, setting this too '2' will double the size of the plotted rt-window. This does not effect the rt region searched which is defined by 'rt' and 'rtspan'.
}
  \item{diagnostic}{
  A named numeric list of fragment masses to be highlighted in the plotted spectra.
}
  \item{diagnostic_ppm}{
  A single numeric equal setting the mass tolerence for the highlighted fragments (defaults to ±15 ppm).
}
  \item{NL}{
  A named numeric list of NL masses to be highlighted in the plotted spectra.
}

}
\details{
Using the 'peakdata' argumnet overides the 'mz' and 'rt' arguments. Use 'peakdata' to search for multipule features spectra at once. Use 'mz' and 'rt' to search for a single features spectra.
}
\value{
Returns a \code{data.frame} or a list of \code{data.frame} if multipule peaks were searched with at once.
}

\examples{
#Example - find MS2 of DNPPE
DNPPE_ms2 <- LOB_findMS2(rawSpec = rawSpec,
                         mz = 875.550487,
                         rt = 1000,
                         ppm = 2.5)
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Henry Holm, hholm@mit.edu \cr
Dan Lowenstien, dlowenstein@whoi.edu \cr
Max Jahns, mjahns@mit.edu \cr
Shavonna Bent, sbent@mit.edu
}
\note{
Currently in development.
}
