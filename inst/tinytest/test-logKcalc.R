library(logKcalc)

# Function extracted from individual test_that calls 20200625
getlines <- function(file) {
  lines <- readLines(file)
  # Normalize “ and ” to "
  # (straight quote is created in R CMD check -- locale setting??)
  lines <- gsub('“', '"', lines)
  lines <- gsub('”', '"', lines)
  rmspecies <- function(pattern, lines) {
    iname <- grep(pattern, lines)
    if(length(iname)==1) {
      iallrefs <- grep("^\\*\\ \\[", lines)
      irefline <- min(iallrefs[iallrefs > iname])
      lines <- lines[-(iname:(irefline + 1))]
    }
    lines
  }
  if(utils::packageVersion("CHNOSZ") < "1.4.0") {
    # References block isn't available in CHNOSZ < 1.4.0 20200625
    iReferences <- match("* References", lines)
    if(!is.na(iReferences)) lines <- lines[-(iReferences:length(lines))]
    # List unavailable minerals and aqueous speices; H2O, H+, e- because they had no references in OBIGT
    patterns <- c("^Wustite", "^Dickite", "^Arsenopyrite", "^FeCO3", "^FeHCO3\\+", "^FeSO4", "^NaCO3\\-", "NaHCO3", "^H2O", "^H\\+", "^e\\-")
    for(pattern in patterns) lines <- rmspecies(pattern, lines)
  }
  if(utils::packageVersion("CHNOSZ") < "1.3.4") {
    # Changed dawsonite to use Joules
    lines <- rmspecies("^Dawsonite", lines)
  }
  # Exclude header lines (timestamp, package version and unavailable species might change)
  lines <- lines[-(7:12)]
  # exclude mineral and aqueous species counts
  lines <- lines[!grepl("minerals$", lines)]
  lines[!grepl("aqueous\\ species$", lines)]
}

info <- "Adding duplicate species produces an error"
infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
iAuCl2 <- info("AuCl2-")
expect_error(logKcalc(infile, ispecies = iAuCl2), "duplicated aqueous species")

info <- "Modifying the database and adding species to the output work as expected"
# Process the thermo_12elements.tdat file
infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
outfile <- file.path(tempdir(), "thermo_12OBIGT.tdat")
#outfile <- "thermo_12OBIGT.tdat"
reset()
modOBIGT(c("addSUPCRT", "steam"))
addOBIGT("AuCl4-")
# Add a selection of 5 aqueous, 1 mineral and 1 gas species
ispecies <- info(c("glycinate", "C2H4", "AuCl", "Au(OH)2-", "Au(HS)2-", "pyrrhotite", "ammonia"))
# Set a non-default ion size parameter for the ions
a0_ion <- 3.5
# Use particular ion size parameters for each neutral aqueous species
a0_neutral <- c(NA, -0.5, 1, NA, NA, NA, NA)
logKcalc(infile, outfile, ispecies = ispecies, a0_ion = a0_ion, a0_neutral = a0_neutral)
reffile <- system.file("extdata/tests/thermo_12OBIGT.tdat", package = "logKcalc")
reflines <- getlines(reffile)
outlines <- getlines(outfile)
expect_identical(outlines, reflines)

info <- "Changing the temperature, water model, and Debye-Hückel method work as expected"
# Process the thermo_12elements.tdat file
infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
outfile <- file.path(tempdir(), "thermo_12OBIGT_bgamma.tdat")
#outfile <- "thermo_12OBIGT_bgamma.tdat"
reset()
modOBIGT(c("addSUPCRT", "steam"))
addOBIGT("AuCl4-")
# Add a selection of aqueous, mineral and gas species
ispecies <- info(c("glycinate", "C2H4", "AuCl", "Au(OH)2-", "Au(HS)2-", "pyrrhotite", "ammonia"))
# Set T and P (1000 bar is the minimum pressure for the DEW model in CHNOSZ)
water("DEW")
T <- seq(300, 650, 50)
P <- 1000
logKcalc(infile, outfile, T, P, ispecies = ispecies, DH.method = "bgamma")
reffile <- system.file("extdata/tests/thermo_12OBIGT_bgamma.tdat", package = "logKcalc")
reflines <- getlines(reffile)
outlines <- getlines(outfile)
# Note: 650 degC, 1000 bar is out of the applicable range of HKF (low-density region),
# so this is also a test that we get "500" values for the last T,P pair
expect_identical(outlines, reflines)

# The next test depends on CHNOSZ >= 1.4.0 (has As(OH)3 from PPB+08) 20201012
info <- "Processing a K2GWB file works as expected"
# Added this test to make sure the Methane(g) reaction is correct 20200625
infile <- system.file("extdata/ThermoGWB_15_6_2020.tdat", package = "logKcalc")
outfile <- file.path(tempdir(), "ThermoGWB_OBIGT.tdat")
#outfile <- "ThermoGWB_OBIGT.tdat"
reset()
# Setup as in vig1.Rmd
modOBIGT(c("addSUPCRT", "steam", "realgar*4"))
logKcalc(infile, outfile)
reffile <- system.file("extdata/tests/ThermoGWB_OBIGT.tdat", package = "logKcalc")
reflines <- getlines(reffile)
outlines <- getlines(outfile)
expect_identical(outlines, reflines)
