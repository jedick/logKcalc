library(logKcalc)

info <- "Adding duplicate species produces an error"
infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
iAuCl2 <- info("AuCl2-")
expect_error(logKcalc(infile, ispecies = iAuCl2), info = info)

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
calcdat <- unlist(readlogK(outfile))
refdat <- unlist(readlogK(reffile))
expect_equal(calcdat, refdat, tolerance = 1e-4, info = info)

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
calcdat <- unlist(readlogK(outfile))
refdat <- unlist(readlogK(reffile))
expect_equal(calcdat, refdat, tolerance = 1e-4, info = info)

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
calcdat <- unlist(readlogK(outfile))
refdat <- unlist(readlogK(reffile))
expect_equal(calcdat, refdat, tolerance = 1e-4, info = info)
