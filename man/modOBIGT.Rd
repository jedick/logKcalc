\encoding{UTF-8}
\name{modOBIGT}
\alias{modOBIGT}
\alias{addOBIGT}
\title{Modify the Thermodynamic Database to Prepare for \logK Calculations}
\description{
These functions modify the thermodynamic database by adding data from optional data files in CHNOSZ or thermodynamic parameters fitted to equilibrium constants from a GWB file.
}

\usage{
  modOBIGT(mod)
  addOBIGT(species, formula = NULL, file = system.file("extdata/thermo_24swapped.tdat",
           package = "logKcalc"), tolerance = 0.1)
}

\arguments{
  \item{mod}{Name(s) of modifications to make to OBIGT}
  \item{species}{Name of a species to add from a GWB file}
  \item{formula}{Chemical formula of the species for OBIGT}
  \item{file}{GWB file that has the species as a redox couple, aqueous or mineral species}
  \item{tolerance}{Tolerance for checking equivalence of input and calculated \logK values}
}

\details{

Before running \code{\link{logKcalc}}, \code{modOBIGT} can be used to make preset modifications to the default OBIGT database, for example to use a different dataset for minerals.
If you require a redox, aqueous, or mineral species from a GWB file that is not present in OBIGT, \code{addOBIGT} can be used to add it.

\code{modOBIGT} takes one or more names of modfications in \code{mod} and runs them in the order given.
These are the available mods:

\describe{
  \item{\samp{addSUPCRT}}{
    The main source of mineral data for OBIGT is Berman (1988), but selected minerals from SUPCRT92, \dQuote{including native elements, sulfides, halides, sulfates, and carbonates and oxides that do not duplicate those in the Berman dataset}, are also included (Dick, 2019).
    This mod add the remaining minerals from SUPCRT92 that are not in OBIGT.
  }
  \item{\samp{allSUPCRT}}{
    This mod adds all minerals from SUPCRT92, possibly replacing the Berman ones.
  }
  \item{\samp{noBerman}}{
    This can be used before either of the above modifications to remove the Berman minerals from the database and generate a set of minerals that is mostly based on SUPCRT92.
  }
  \item{\samp{steam}}{
    This adds data for steam (\H2O(g)) taken from the \samp{sprons92.dat} file of SUPCRT92, which is referenced there as \dQuote{Johnson, J.W., Personal calculation, Earth Sciences Dept., Lawrence Livermore National Lab., Livermore, CA. Parameters given provide smooth metastable extrapolation of one-bar steam properties predicted by the Haar et al. (1984) equation of state to temperatures < the saturation temperature (99.632 C).}
    This entry was removed in CHNOSZ 1.1.3 based on the issues raised by Zavarin et al. (2016, p. 6).
    However, many existing GWB files have an entry for steam, and the data are needed to include it in the output.
  }
  \item{\samp{realgar*4}}{
    This mod scales the chemical formula and thermodynamic parameters of realgar by a factor of 4.
    This makes realgar in OBIGT compatible with the UNITHERM database in the HCh package (Shvarov and Bastrkov, 1999).
  }
  \item{\samp{antigorite/2}}{
    This mod scales the chemical formula and thermodynamic parameters of antigorite by a factor of 1/2.
    This makes antigorite in OBIGT compatible with the thermo.tdat file from GWB.
  }
}

\code{addOBIGT} reads the equilibrium constants for the dissociation reaction of one \code{species} from the indicated \code{file}.
The standard Gibbs energy of the species is computed at each temperature from the \logK values and the Gibbs energies of the other species in the reaction, which must be available in OBIGT.
Then, the standard Gibbs energy of the species is fit using the thermodynamic parameters \emph{G}, \emph{S}, and \emph{Cp}.

The species is then added to OBIGT; the upper temperature limit (corresponding to the \T range for the fitted \logK values) is set in the \samp{T} parameter for minerals or the \samp{abbrv} parameter for aqueous species; the latter is specially used by other functions in this package to assign NA values at appropriate places (i.e., \samp{500} in GWB data files).
An attempt is made to convert GWB-style formulas (e.g. Sn\S{++++}) to the format used in OBIGT (Sn\S{+4}) (see Names Mapping section in \code{\link{logKcalc}}), but if there is any doubt, the \code{formula} should be given explicitly.
For adding minerals, a valid formula must be provided in the \code{formula} argument.

Note that the equation assumes a constant heat capacity and no volume change, so in general does not perfectly fit the \logK values in the input GWB file and may be unsuitable for extrapolation beyond the original \T and \P range.
\code{\link{all.equal}} is used to test the near equivalence of the input and calculated equilibrium constants; if the differences exceed \code{tolerance}, the code generates an error.

}

\seealso{
See \code{\link[CHNOSZ]{add.OBIGT}} in \pkg{CHNOSZ} for more information about modifying the database.
}


\examples{
## Useful settings for processing the GWB thermo.tdat file
modOBIGT(c("addSUPCRT", "steam", "antigorite/2"))

# As(OH)4- and Sn++++ are basis species in the thermo.tdat GWB file,
# but aren't available in OBIGT, so we add them like this:
addOBIGT("As(OH)4-")
addOBIGT("Sn++++")
# Now As- and Sn-bearing species can be included in the logK calculation:
#logKcalc("thermo.tdat")
}

\references{
Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \doi{10.1093/petrology/29.2.445}

Dick, J. M. (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{http://pid.geoscience.gov.au/dataset/ga/25473}

Zavarin, M., Wolery, T. J. and Atkins-Duffin, C. (2016) M4FT-16LL080302052-Update to Thermodynamic Database Development and Sorption Database Integration. Lawrence Livermore National Laboratory, LLNL-TR-701407. \doi{10.2172/1325873}
}
