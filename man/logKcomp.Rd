\encoding{UTF-8}
\name{logKcomp}
\alias{logKcomp}
\alias{readlogK}
\title{Read and Compare Equilibrium Constants in GWB Data Files}
\description{
Read values of logK from a data file, and make plots comparing equilibrium constants from two files.
}

\usage{
  logKcomp(file1, file2, type = c("redox", "aqueous"), iTP = 2,
    lab1 = NULL, lab2 = NULL, xlim = NULL, ylim = NULL, plot.it = TRUE)
  readlogK(file, quiet = FALSE)
}

\arguments{
  \item{file1}{GWB data file}
  \item{file2}{GWB data file}
  \item{type}{Type of species (\samp{redox}, \samp{aqueous}, \samp{mineral}, or \samp{gas})}
  \item{iTP}{Which pair of \T, \P values to use (1..8)}
  \item{lab1}{Title for \code{file1} in plot labels}
  \item{lab2}{Title for \code{file2} in plot labels}
  \item{xlim}{Limits for \emph{x} axis}
  \item{ylim}{Limits for \emph{y} axis}
  \item{plot.it}{Logical, make a plot?}
  \item{file}{GWB data file}
  \item{quiet}{Logical, reduce messages and printed output?}
}

\details{

\code{logKcomp} compares \logK values from two GWB thermodynamic data files.
The values are taken from the same index of \T, \P pairs in the files.
The default for \code{iTP} corresponds to 25 \degC and 1 bar for commonly used data files, but it is the user's responsibility to ensure that the files have comparable \T, \P grids.

The values for each indicated \code{type} of species are read from \code{file1} and \code{file2} and combined.
By default, data for both redox couples and aqueous species are used, but \code{type} can be changed to one (or a combination of) \code{redox}, \code{aqueous}, \code{electron}, \code{mineral}, and \code{gas}.

The function includes only those species that have the same names in both files, or whose \emph{mapped} names (see \code{\link{logKcalc}}) are the same.
Species are also excluded from the plot (or the returned data frame) if their dissociation reactions are not identical in the two files.

The function returns a \code{\link[ggplot2]{ggplot}} object, which normally results in a plot being made on screen when run interactively.
Set \code{plot.it} to FALSE to return a data frame of \logK values from the two files instead of making a plot.
If \code{file2} is missing, the function doesn't make a plot and instead returns a numeric vector of the \logK values from \code{file1}.

A scatterplot is made by plotting values from \code{file1} on the \emph{x} axis and the difference between values in \code{file1} and \code{file2} on the \emph{y} axis.
Points with differences greater than 1 \logK unit are labeled.
The axis labels are generated from the \code{\link{basename}} of the files, but this can be overriden with \code{lab1} and \code{lab2}.

\code{readlogK} reads all of the \logK values from the specified \code{file} and returns a list with numeric values of \T and \P and lists of \logK values for each type of species.

}

\examples{
# Read logKs from a data file
file <- system.file("extdata/thermo_24elements.tdat", package = "logKcalc")
LOGK <- readlogK(file)
# Show the number of entries
sapply(LOGK, length)

# Compare logKs in two data files
file1 <- system.file("extdata/thermo_24noredox.tdat", package = "logKcalc")
file2 <- system.file("extdata/ThermoGWB_15_6_2020.tdat", package = "logKcalc")
logKcomp(file1, file2)
}
