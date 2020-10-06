\name{LOB_viewdata}
\alias{LOB_viewdata}
\title{
Vizualize lipids annotated by LOBSTAHS with live graphs in shiny
}
\description{
This function launches a \code{shiny} app to graph detected lipids in retention time verse m/z space. Graphs can be colored based on multipul screening criteria in order visualize data and identify correct or suspect annoations.

}
\usage{
LOB_viewdata(peakdata, rawSpec = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peakdata}{
  A object of class \code{LOBSet} or \code{data.frame} generated with \code{getLOBpeaklist()} that contains the peakdata of identifed lipids.
}
  \item{rawSpec}{
  If the user chooses, a \code{XCMSnExp} object to be searched contianing both MS1 and MS2 spectra that the peakdata was generated from. \code{LOB_findMS2()} can be run within the \code{shiny} app to detect MS2 spectra that may correspond with annotated peaks.
}

}
\details{
Runs a \code{shiny} app that allows for multipul visualizations of the LOBSTAHS generated lipid annotations. Annotations can be viewed by class, double bonds, carbon number, oxidation, retantion time, screening criteria, and more. In addition to this, the GUI allows for easy passing of mz and rt info to \code{LOB_findMS2()}. Peaks may be selected by clicking to display their metadata and multipul annotations in a table. MS2 spectra of annotated peaks may be identifed elsewhere using the the output of \code{LOB_findMS2()} which is displayed in app.
}

\value{
A \code{shiny} app.
}

\examples{
# View PtH2O2lipids annotation LOBSet
LOBSet <- PtH2O2lipids::ptH2O2lipids$LOBSet
LOB_viewdata(LOBSet)
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

