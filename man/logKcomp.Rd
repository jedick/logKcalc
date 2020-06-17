\encoding{UTF-8}
\name{logKcomp}
\alias{logKcomp}
\alias{readlogK}
\title{Read and Compare Equilibrium Constants in GWB Data Files}
\description{
Read values of logK from a data file, and make plots comparing equilibrium constants from two files.
}

\usage{
  readlogK(file, quiet = FALSE)
  logKcomp(file1, file2, type = "aqueous", iTP = 2,
    lab1 = NULL, lab2 = NULL, plot.it = TRUE)
}

\arguments{
  \item{file}{GWB data file}
  \item{quiet}{Logical, reduce messages and printed output?}
  \item{file1}{GWB data file}
  \item{file2}{GWB data file}
  \item{type}{Type of species (\samp{redox}, \samp{aqueous}, \samp{mineral}, or \samp{gas})}
  \item{iTP}{Which pair of \T, \P values to use (1..8)}
  \item{lab1}{Title for \code{file1} in plot labels}
  \item{lab2}{Title for \code{file2} in plot labels}
  \item{plot.it}{Logical, make a plot?}
}

\details{

\code{readlogK} reads \logK values from the specified \code{file} and returns a list with elements named T, P, redox, aqueous, electron, mineral, gas.
In this list, T and P are numeric values read from the file, and \logK values are listed for the different species types.

\code{logKcomp} compares \logK values from two files.
The values are read from \code{file1} and \code{file2}, and selected for the indicated \code{type} of species and pair of \T and \P values.
The default for \code{iTP} corresponds to 25 \degC and 1 bar for commonly used data files.

The function prints the names of non-matching species between the two files.
Species are also excluded from the plot (or the returned data frame) if their dissociation reactions are not identical in the two files.

A scatterplot is made by plotting values from \code{file1} on the x axis and the difference between values in \code{file1} and \code{file2} on the y axis.
Points with differences greater than 1 \logK unit are labeled.
The axis labels are generated from the \code{\link{basename}} of the files, but this can be overriden with \code{lab1} and \code{lab2}.

The function returns a \code{\link[ggplot2]{ggplot}} object, which normally results in a plot being made on screen when run interactively.
Set \code{plot.it} to FALSE to return a data frame of \logK values from the two files instead of making a plot.

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
