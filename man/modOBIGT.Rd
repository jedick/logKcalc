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
           package = "logKcalc"), tolerance = 0.05)
}

\arguments{
  \item{mod}{Name(s) of modifications to make to OBIGT}
  \item{species}{Name of a species to add from a GWB file}
  \item{formula}{Chemical formula of the species for OBIGT}
  \item{file}{GWB file that has the species as a redox or aqueous species}
  \item{tolerance}{Tolerance for checking equivalence of input and calculated \logK values}
}

\details{

Before running \code{\link{logKcalc}}, you may wish to use \code{modOBIGT} to include some minerals not in OBIGT by default.
The main source of mineral data for OBIGT is Berman (1988), but selected minerals from SUPCRT92, \dQuote{including native elements, sulfides, halides, sulfates, and carbonates and oxides that do not duplicate those in the Berman dataset}, are also included (Dick, 2019).
Use \code{modOBIGT("addSUPCRT")} to add the remaining minerals from SUPCRT92 that are not in OBIGT.
\code{modOBIGT("allSUPCRT")} adds all minerals from SUPCRT92, possibly replacing the Berman ones.
Either of these commands can be preceeded by \code{modOBIGT("noBerman")} to remove the Berman minerals from the database and generate a set of minerals that is mostly compatible with SUPCRT92.

\code{modOBIGT("steam")} adds data for steam (H2O(g)) taken from the \samp{sprons92.dat} file of SUPCRT92, which is referenced there as \dQuote{Johnson, J.W., Personal calculation, Earth Sciences Dept., Lawrence Livermore National Lab., Livermore, CA. Parameters given provide smooth metastable extrapolation of one-bar steam properties predicted by the Haar et al. (1984) equation of state to temperatures < the saturation temperature (99.632 C).}
This entry was removed in CHNOSZ 1.1.3 based on the issues raised by Zavarin et al. (2016, p. 6).
However, many existing GWB files have an entry for steam, and the data are needed to include it in the output.

\code{modOBIGT("realgar*4")} and \code{modOBIGT("antigorite/2")} scale the properties of realgar and antigorite by the indicated factors.
These changes are appropriate for the Unitherm database in the HCh package (Shvarov and Bastrkov, 1999) and for the thermo.tdat file from GWB, respectively.

If you require a redox or aqueous species from a GWB file that is not present in OBIGT, \code{addOBIGT} can be used to add it.
This function reads the equilibrium constants for the dissociation reaction of the species from the GWB file.
Then, the standard Gibbs energy of the species is computed at each temperature from the \logK values and the Gibbs energies of the other species in the reaction, which must be available in OBIGT.
The standard Gibbs energy of the species is fit using the thermodynamic parameters G, S, and Cp.

The species is then added to OBIGT.
An attempt is made to convert GWB-style formulas (e.g. \samp{Sn++++}) to an acceptable format for OBIGT (\samp{Sn+4}) (see name mapping heuristics in \code{\link{logKcalc}}), but if there is any doubt, the \code{formula} should be given explicitly.

Note that the equation assumes a constant Cp and no volume change, so is not appropriate for some types of calculations.
However, it is able to reproduce the equilibrium constants in the \samp{thermo.tdat} GWB file within a reasonable error.
\code{\link{all.equal}} is used to test the equivalence of the input and calculated equilibrium constants; if the differences exceed \code{tol}, the code generates an error.

}

\seealso{
See \code{\link[CHNOSZ]{add.obigt}} in \pkg{CHNOSZ} for more information about modifying the database.
}


\examples{
## Useful settings for processing the GWB thermo.tdat file
modOBIGT(c("addSUPCRT", "steam", "antigorite/2"))
# As(OH)4- and Sn++++ are basis species in the thermo.tdat GWB file,
# but aren't available in OBIGT, so we add them here
addOBIGT("As(OH)4-")
addOBIGT("Sn++++")
# Now the file can be processed to include As- and Sn-bearing
# species (see the thermo.Rmd vignette)
#logKcalc("thermo.tdat")
}

\references{
Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \doi{10.1093/petrology/29.2.445}

Dick, J. M. (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{http://pid.geoscience.gov.au/dataset/ga/25473}

Zavarin, M., Wolery, T. J. and Atkins-Duffin, C. (2016) M4FT-16LL080302052-Update to Thermodynamic Database Development and Sorption Database Integration. Lawrence Livermore National Laboratory, LLNL-TR-701407. \doi{10.2172/1325873}
}
