\name{LOB_RTFsort}
\alias{LOB_RTFsort}

\title{
Screen annoations against standards retention time factors and previous MS2 verifications.
}

\description{
This function allows the user screen lipid annotations against known retention time factors (RTFs), that is the ratio between the lipids retention time and a standards retention time. While retention times in LC-MS may have large shifts, RTFs are largely fixed as peaks usually drift in similar directions with instrument changes. Thus, RTFs provide a metric of retention time consitancy that may be used between datasets for identify compounds repeatably.
}

\usage{
LOB_RTFsort(peakdata, RT_Factor_Dbase, choose_class = NULL, \cr            plot_data = FALSE, save_plots = FALSE, data_title)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peakdata}{
  A object of class \code{LOBSet} or \code{data.frame} generated with \code{getLOBpeaklist()} that contains the peakdata of identifed lipids.
}
  \item{RT_Factor_Dbase}{
  A \code{data.frame} containing standard information and RTFs.
}
  \item{choose_class}{
  A \code{vector} of classes in the peakdata that should be screened. If NULL all classes are screened. 
}
  \item{plot_data}{
  A \code{logical} indicating whether the data should be plotted (TRUE) or should not be (FALSE).
}
  \item{save_plots}{
  A \code{logical} indicating whether the plots should be saved (TRUE) or should not be (FALSE) to the current working directory.
}
  \item{data_title}{
  NA
}
}
\details{

}
\value{
Returns a \code{data.frame} or a \code{LOBSet} object matching the \code{peakdata} argument but with columns containing RTF screening data appended. When run on a \code{LOBSet} object this screening information appended to the \code{data.frame} in the @peakdata slot of the object. The column names are as follows:

  \item{DNPPE_Factor}{
  RTF for the annotation: The retention time of the peakgroup in question divided by the retention time of DNP-PE.
}

  \item{DBase_DNPPE_RF}{
  Database RTF for this annotation: A previously observed RTF for this annotation.
}

  \item{Flag}{
  A flag denoting how close the actual RTF is too the database value. If the two values are within 5\% of eachother the tag "5\%_rtv" is used. Also, if an appropriate MS2 spectra for this compound has been seen at this RTF previously, then it receives the tag "ms2v". If the peakgroup is within 5-10\% of the database RTF than it is given the tag "10\%_rtv". If the RTF is greater than 10\% removed from the database value it is given the tag "Red". A tag of "Unknown" is given to every peakgroup for which there is no database value.
}

}

\examples{
# Screen a LOBSet
LOBSet <- PtH2O2lipids::ptH2O2lipids$LOBSet
LOB_RTFsort(LOBSet, RT_Factor_Dbase, choose_class = NULL, plot_data = TRUE, save_plots = FALSE, data_title)

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
