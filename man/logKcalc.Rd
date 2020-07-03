\encoding{UTF-8}
\name{logKcalc}
\alias{logKcalc}
\title{Calculate Equilibrium Constants for a GWB Data File}
\description{
Read dissociation reactions from a thermodynamic data file for The Geochemist's Workbench®, calculate equilibrium constants at specified temperature and pressure using the OBIGT database in CHNOSZ, and write a modified GWB file.
}

\usage{
  logKcalc(infile = "thermo.tdat", outfile = "thermo_OBIGT.tdat",
    T = NULL, P = "Psat", ispecies = NULL, a0_ion = NULL, a0_neutral = 0,
    update.formulas = TRUE, DH.method = "bdot", maxprint = Inf)
}

\arguments{
  \item{infile}{Input GWB data file}
  \item{outfile}{Output GWB data file}
  \item{T}{Temperature (\degC)}
  \item{P}{Pressure (bar)}
  \item{ispecies}{Species index of specie(s) to add to the data file}
  \item{a0_ion}{Ion size parameter for added aqueous ions}
  \item{a0_neutral}{Ion size parameter for added aqueous neutral species}
  \item{update.formulas}{Should formulas for aqueous species and minerals be updated?}
  \item{DH.method}{Method for generating Debye-Hückel parameters, \samp{bdot} or \samp{bgamma}}
  \item{maxprint}{Maximum number of missing species to show in print statements}
}

\details{

This function reads the input GWB file and assembles lists of species in the file.
Species names are then converted to those used in \link[CHNOSZ:thermo]{OBIGT}, by automatically matching the species name or by looking in a mapping file (see Names Mapping).
Species that are not available in OBIGT are omitted from the output file.

Species are also omitted if they available in OBIGT but their dissociation reactions in the GWB file are written in terms of basis species or redox couples that are not available.
Two examples are As(OH)\s4\S{-} and Sn\S{++++}, which are basis species in the \samp{thermo.tdat} file distributed with GWB but are not present in the default OBIGT database.
\code{\link{addOBIGT}} can be used to add these species to OBIGT.

Normally, the function obtains values of \T from the input GWB file and calculates values of \P corresponding to \Psat (i.e. pressures corresponding to liquid-vapor saturation).
This is done because some values of \P in existing files are slightly below \Psat, preventing calculation at those conditions.
To use values of \P taken from the input file, set \code{P} to NULL.
To use any other values, provide 8 numeric values in \code{T} and/or \code{P}.

Logarithms (base 10) of equilibrium constants are calculated for each of the following species types: redox, aqueous, electron, mineral, gas.
The corresponding values for each species are inserted into the output GWB file.
Any NA values (e.g. for temperatures above the melting point of a mineral) are represented by \samp{500} in the output file to indicate unavailable data (Bethke and Farrell, 2020).

One or two citations are added for each species based on the references in OBIGT.
If the CHNOSZ version is > 1.3.6, a reference list is added at the end of the output file.

}

\section{Chemical formulas and reactions}{

By default, the stoichiometric formulas in OBIGT are added to the \samp{formula} field in the output for aqueous species and minerals.
These formulas may be formatted differently from formulas in existing GWB data files, especially for charged species.
Set \code{update.formulas} to FALSE to not add or update formulas, but retain any formulas that may be present in the input file.
This might be required for dataset formats \samp{oct13} or earlier, which do not have formulas for the aqueous species (Bethke and Farrell, 2020).
Note that the \emph{names} of aqueous species, which may also be formulas, are never changed from the input file.

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
Both names and formulas of species in OBIGT are searched, but formulas are searched as-is, so formulas that differ from those used in CHNOSZ (e.g. CH{\s3}COOH instead of C{\s2}H{\s4}O{\s2}) aren't matched unless they are in the names mapping file.

Species names that are not mapped after these steps are then searched in the names mapping file.
This file was constructed by identifying unmapped species in \samp{thermo.tdat} and GWB files produced by K2GWB and DBCreate.
Some names (for feldspars in particular) have multiple entries in this file, to allow mapping to names used in the Berman or SUPCRT92 datasets, which can be selected using \code{\link{modOBIGT}}.

}

\section{Adding Species}{

By default, only species in the input file that can be mapped to OBIGT are included in the output file.
Other species can be added to the output using the \code{ispecies} argument.
The values in this argument are the species index for one or more species, as returned by \code{CHNOSZ::\link[CHNOSZ]{info}}.
Species can only be added if all of their elements are in the available basis species in the GWB file.
Species with the states \samp{aq}, \samp{cr}, and \samp{gas} in OBIGT are added to the aqueous, mineral, and gas sections of the output file, respectively; any other states produce an error.

The molecular weight of species is calculated using \code{CHNOSZ::\link[CHNOSZ:util.formula]{mass}}, which likely does not correspond precisely to the element masses in the input GWB file.
The default ion size parameter for added aqueous ions (\code{a0_ion}) is 4.5 \AA; this setting follows the behavior in UNITHERM for ions without a specified value (Shvarov and Bastrakov, 1999).

The ion size parameter for neutral aqueous species (\code{a0_neutral}) has a special meaning in the GWB thermodynamic data file (Bethke and Farrell, 2020).
The default value of 0 corresponds to unit activity coefficient (the default in the UT2K converter used by K2GWB; Cleverley and Bastrakov, 2005).
With a value of -0.5, the activity coefficient is equal to that of CO\s2 in an NaCl solution of the same ionic strength (this is used for B(OH)\s3, O\s2, SiO\s2, CH\s4, and H\s2 in \samp{thermo.tdat}).
For values \le -1, the activity coefficient is calculated as the product of the \samp{bdot} parameter and the true ionic strength.

The values of \code{a0_ion} and \code{a0_neutral} are recycled to the total number of species, including all aqueous ions and neutral species as well as non-aqueous species.
Therefore, to specify particular values of these parameters for different ions or neutral aqueous species, values must be specified for \emph{all} added species, but will only be used for the respective aqueous species.

}

\section{Activity and Osmotic Coefficients}{

Coefficients in the extended Debye-Hückel equation (\samp{adh} and \samp{bdh} in the GWB file) are calculated as a function of \T and \P using \code{CHNOSZ::\link[CHNOSZ]{water}}.
The extended term (\samp{bdot} parameter) by default is calculated using a spline fit to values at 25-300 \degC from Helgeson (1969); values at \T > 300 \degC are set to zero.

The \code{DH.method} argument can be changed to \samp{bgamma} to use the b-gamma approximation, in which the extended term parameter is a function of temperature and pressure that has been extrapolated to 1000 \degC and 30 kbar (Manning et al., 2013) and provisionally as high as 60 kbar (see Figure S5 of Dick, 2019).
The high-pressure extrapolation is intended to be used in conjunction with the Deep Earth Water (DEW) model for water.
The default for \samp{bgamma} is to set the ion size parameter for \emph{all} ions in the output file to 3.72 \AA, which is a characteristic value for an NaCl-dominated solution (Helgeson et al., 1981).
This differs from HCh, where ion size parameters for individual ions are calculated from their effective elecrostatic radii (Shvarov and Bastrakov, 1999).
Calling the function with the explicit argument \code{a0_ion = NULL} overrides this setting, so that no changes are made to ion size parameters from the input file.

Parameters for calculating the activity coefficient of CO\s2 in a pure NaCl solution and the osmotic coefficient of \H2O (i.e., \samp{co2} and \samp{h2o} blocks of the GWB file) are calculated as described for K2GWB (Cleverley and Bastrakov, 2005).
Because of a bug with calculating exponents in some versions of K2GWB, the CO\s2 parameters (especially in the first block) may not match those in GWB files produced by K2GWB for the same values of \T.

}

\section{Dataset Formats}{

The conversion is compatible with all existing dataset formats of GWB files.
For dataset formats \samp{feb94} and \samp{oct94}, \code{logKcalc} updates the \logK blocks for the Eh reaction and O\s2 gas solubility.
(Blocks for H\s2 and N\s2 gas solubility, which are not used by GWB, are not updated.)

These blocks are not present in dataset formats \samp{oct13} and later.
Instead, the free electron is present as a species; \code{logKcalc} updates its \logK values.

Dataset formats \samp{jul17} and \samp{jan19} specify a fugacity model (tsonopoulos, peng-robinson, or spycher-reed).
The parameters are located in the data blocks for individual gas species and are not modified by \code{logKcalc}.

}

\seealso{
\code{\link{modOBIGT}} for preparing the thermodynamic database for the \logK calculation.
}

\examples{
# Show the first few entries in the names mapping file
mapfile <- system.file("extdata/map_names.csv", package = "logKcalc")
head(read.csv(mapfile))

# Process a GWB file
# -- define input and output files
infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
outfile <- tempfile("thermo", fileext = ".tdat")
# -- make some modifications to OBIGT
# (add minerals from SUPCRT92 database and deprecated H2O(g))
modOBIGT(c("addSUPCRT", "steam"))
# -- identify species to add to the output
ispecies <- info("pyrrhotite")
# -- process the file
logKcalc(infile, outfile, ispecies = ispecies)
}

\references{
Bethke, C. M. and Farrell, B. (2020) The Geochemist's Workbench® Release 14 GWB Reference Manual, Aqeuous Solutions, LLC, Champaign, Illinois.

Cleverley, J. S. and Bastrakov, E. N. (2005) K2GWB: Utility for generating thermodynamic data files for The Geochemist's Workbench® at 0-1000 °C and 1-5000 bar from UT2K and the UNITHERM database. \emph{Computers & Geosciences} \bold{31}, 756--767. \doi{10.1016/j.cageo.2005.01.007}

Dick, J. M. (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Helgeson, H. C., Kirkham, D. H. and Flowers, G. C. (1981) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. IV. Calculation of activity coefficients, osmotic coefficients, and apparent molal and standard and relative partial molal properties to 600\degC and 5 Kb. \emph{Am. J. Sci.} \bold{281}, 1249--1516. \doi{10.2475/ajs.281.10.1249}

Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \doi{10.2138/rmg.2013.75.5}

Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{http://pid.geoscience.gov.au/dataset/ga/25473}
}
