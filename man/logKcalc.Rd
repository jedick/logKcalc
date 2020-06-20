\encoding{UTF-8}
\name{logKcalc}
\alias{logKcalc}
\title{Calculate Equilibrium Constants for a GWB Data File}
\description{
Read dissociation reactions from a thermodynamic data file for The Geochemist's Workbench®, calculate equilibrium constants at specified temperature and pressure using the OBIGT database in CHNOSZ, and write a modified GWB file.
}

\usage{
  logKcalc(infile = "thermo.tdat", outfile = "thermo_OBIGT.tdat",
    T = NULL, P = "Psat", ispecies = NULL, maxprint = Inf)
}

\arguments{
  \item{infile}{Input GWB data file}
  \item{outfile}{Output GWB data file}
  \item{T}{Temperature (\degC)}
  \item{P}{Pressure (bar)}
  \item{ispecies}{Species index of specie(s) to add to the data file}
  \item{maxprint}{Maximum number of missing species to show in print statements}
}

\details{

This function reads the input GWB file and assembles lists of species in the file.
Species names are then converted to those used in OBIGT, by automatically matching the species name or by looking in a mapping file (see first example).
Species that are not available in OBIGT are omitted from the output file.

Species are also omitted if they available in OBIGT but their dissociation reactions in the GWB file are written in terms of some other species that are not available.
Two particular examples are \samp{As(OH)4-} and \samp{Sn++++}, which are basis species in the default \samp{thermo.tdat} file as distributed with GWB.
This means that \emph{no} As- or Sn-bearing species are available with the default file.

By default, the function obtains values of \T from the input GWB file and calculates values of \P corresponding to \Psat (i.e. pressures corresponding to liquid-vapor saturation).
This is done because some values of \P in existing files are slightly below \Psat, preventing calculation at those conditions.
To use values of \P taken from the input file, set \code{P} to NULL.
To use any other values, provide 8 numeric values in \code{T} and/or \code{P}.

Logarithms (base 10) of equilibrium constants are calculated for each of the following species types: redox, aqueous, electron, mineral, gas.
The corresponding values for each species are inserted into the output GWB file.
Any NA values (e.g. for temperatures above the melting point of a mineral) are represented by \samp{500} in the output file to indicate unavailable data (Bethke and Yeakel, 2010).

The mass balance of dissociation reactions is checked by using the chemical formulas of species in OBIGT.
Any unbalanced reactions produce a warning from \code{\link{subcrt}}.
Species with unbalanced reactions are listed in the output file, but their \logK values are set to \samp{500} if the difference for any element exceeds 0.001.
This includes the minerals antigorite, pyrrhotite, and scorodite in the GWB \samp{thermo.tdat} file, for which the chemical formulas (and therefore reactions) differ from those in OBIGT.
An unbalanced reaction warning is also produced for realgar, not because of a different chemical formula but because of rounded reaction coefficients to the nearest 0.001, so its \logK values are included in the output.

}

\section{Names Mapping}{

The names mapping step uses some heuristics to find species in OBIGT.
All names from the GWB file are converted to lower-case, suffixes \dQuote{(aq)} and \dQuote{,aq} are removed, and underscores are converted to spaces.
Repeated \samp{-} or \samp{+} signs at the ends of names are changed to numbers (e.g. \samp{-2} or \samp{+3}).
Leading \samp{1-} are deleted (e.g. 1-pentanol becomes pentanol).

The names of species in the GWB file, which correspond to chemical formulas for many aqueous species, are used for mapping.
However, any chemical formulas that may be listed in addition to names of minerals or other species are not used for mapping.
Both names and formulas of species in OBIGT are searched, but formulas are searched as-is, so formulas that differ from those used in CHNOSZ (e.g. CH3COOH instead of C2H4O2) aren't matched unless they are in the names mapping file.

Species names that are not mapped after these steps are then searched in the names mapping file.
This file was constructed by identifying unmapped species in \samp{thermo.tdat} and GWB files produced by K2GWB and DBCreate.
Some of the mapped names are for SUPCRT92 minerals that are not in the default database, so they won't be mapped to the output file unless \code{modOBIGT} is used add these minerals.

}

\section{Adding Species}{

By default, only species in the input file that can be mapped to OBIGT are included in the output file.
Other species can be added to the output using the \code{ispecies} argument.
The values in this argument are the species index for one or more species, as returned by \code{CHNOSZ::\link[CHNOSZ]{info}}.
Species can only be added if all of their elements are in the available basis species in the GWB file.
Species with the states \samp{aq}, \samp{cr}, and \samp{gas} in OBIGT are added to the aqueous, mineral, and gas sections of the output file, respectively; any other states produce an error.

The molecular weight of species is calculated using \code{CHNOSZ::\link[CHNOSZ:util.formula]{mass}}, which likely does not correspond precisely to the element masses in the input GWB file.
Ion sizes for all new aqueous species are set to 4.0 A.

}

\seealso{
\code{\link{modOBIGT}} for preparing the thermodynamic database for the \logK calculation.
}

\examples{
# Show the first few entries in the names mapping file
mapfile <- system.file("extdata/map_names.csv", package = "logKcalc")
head(read.csv(mapfile))

# Process a GWB file and save the output in a temporary file
infile <- system.file("extdata/ThermoGWB_15_6_2020.tdat", package = "logKcalc")
outfile <- tempfile("thermo", fileext = ".tdat")
modOBIGT(c("addSUPCRT", "steam", "realgar*4"))
logKcalc(infile, outfile)
}

\references{
Bethke, C. M. and Yeakel, S. (2010) The Geochemist's Workbench® Reference Manual (Release 8.0), Hydrogeology Program, University of Illinois.
}
