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

The mass balance of dissociation reactions is checked (in CHNOSZ::\code{\link{subcrt}}) using the chemical formulas of species in OBIGT.
To allow for rounding errors in fractional coefficients, reactions are considered unbalanced if the difference for any element exceeds 0.001.
Messages are shown for unbalanced reactions, and these species are omitted from the output file.
This includes the minerals antigorite and pyrrhotite in the GWB \samp{thermo.tdat} file, for which the chemical formulas (and therefore reactions) differ from those in OBIGT.
These minerals can be included by using \code{\link{modOBIGT}} and the \code{ispecies} argument, as shown in the example and \code{thermo.Rmd} vignette.

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
Ion sizes for new aqueous species are set to 4.5 A; this setting follows the behavior in UNITHERM for ions without a specified value (Shvarov and Bastrakov, 1999).
TODO: specify ion sizes for neutral aqueous species (Cleverley and Bastrakov, 2005)

}

\section{Activity and Fugacity Models}{
Coefficients in the extended Debye-Hückel equation (\samp{adh} and \samp{bdh} in the GWB file) are calculated as a function of \T and \P using \code{CHNOSZ::\link[CHNOSZ]{water}}.
The extended term (\samp{bdot} parameter) by default is calculated using a spline fit to values at 25-300 \degC from Helgeson (1969); values at \T > 300 \degC are set to zero.

TODO: The \code{DH.method} argument can be changed to \samp{bgamma} to use the b-gamma approximation, in which the ion size parameter is a constant for all ions, and the extended term parameter is a function of temperature and pressure that has been extrapolated to 1000 \degC and 30 kbar (Manning et al., 2013) and provisionally as high as 60 kbar (see Figure S5 of Dick, 2019).
The high-pressure extrapolation is intended to be used in conjunction with the Deep Earth Water (DEW) model for water.
}

\seealso{
\code{\link{modOBIGT}} for preparing the thermodynamic database for the \logK calculation.
}

\examples{
# Show the first few entries in the names mapping file
mapfile <- system.file("extdata/map_names.csv", package = "logKcalc")
head(read.csv(mapfile))

# Process a GWB file
# 1) define input and output files
infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
outfile <- tempfile("thermo", fileext = ".tdat")
# 2) make some modifications to OBIGT
# addSUPCRT - minerals from SUPCRT92 database
# steam - deprecated H2O(g)
modOBIGT(c("addSUPCRT", "steam"))
# 3) identify species to add to the output
ispecies <- info("pyrrhotite")
# 4) process the file
logKcalc(infile, outfile, ispecies = ispecies)
}

\references{
Bethke, C. M. and Yeakel, S. (2010) The Geochemist's Workbench® Reference Manual (Release 8.0), Hydrogeology Program, University of Illinois.

Cleverley, J. S. and Bastrakov, E. N. (2005) K2GWB: Utility for generating thermodynamic data files for The Geochemist's Workbench® at 0-1000 °C and 1-5000 bar from UT2K and the UNITHERM database. \emph{Computers & Geosciences} \bold{31}, 756--767. \doi{10.1016/j.cageo.2005.01.007}

Dick, J. M. (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \doi{10.2138/rmg.2013.75.5}

Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{http://pid.geoscience.gov.au/dataset/ga/25473}
}
