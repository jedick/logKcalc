context("logKcalc")

test_that("Adding duplicate species produces an error", {
  infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
  iAuCl2 <- info("AuCl2-")
  expect_error(logKcalc(infile, ispecies = iAuCl2), "duplicated aqueous species")
})

test_that("Modifying the database and adding species to the output work as expected", {
  # process the thermo_12elements.tdat file
  infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
  outfile <- file.path(tempdir(), "thermo_12OBIGT.tdat")
  reset()
  modOBIGT(c("addSUPCRT", "steam"))
  addOBIGT("AuCl4-")
  # add a selection of 2 aqueous, 1 mineral and 1 gas species
  ispecies <- info(c("C2H4", "AuCl", "Au(OH)2-", "Au(HS)2-", "pyrrhotite", "ammonia"))
  # set a non-default ion size parameter for the ions
  a0_ion <- 3.5
  # use particular ion size parameters for each neutral aqueous species
  a0_neutral <- c(-0.5, 1, NA, NA, NA, NA)
  logKcalc(infile, outfile, ispecies = ispecies, a0_ion = a0_ion, a0_neutral = a0_neutral)
  reffile <- system.file("extdata/tests/thermo_12OBIGT.tdat", package = "logKcalc")
  reflines <- readLines(reffile)
  outlines <- readLines(outfile)
  # normalize “ and ” to "
  # (straight quote is created in R CMD check -- locale setting??)
  reflines <- gsub('“', '"', reflines)
  reflines <- gsub('”', '"', reflines)
  outlines <- gsub('“', '"', outlines)
  outlines <- gsub('”', '"', outlines)
  # make the test excluding the lines with the timestamp and package versions
  expect_identical(outlines[-(6:7)], reflines[-(6:7)])
#  # Testing for identical files is problematic with ongoing updates to OBIGT,
#  # so just check that the values are close to each other
#  mineral <- logKcomp(reffile, outfile, type = "mineral", plot.it = FALSE)
#  expect_equal(mineral$logK1, mineral$logK2, tolerance = 0.1)
#  aqueous <- logKcomp(reffile, outfile, type = "aqueous", plot.it = FALSE)
#  expect_equal(aqueous$logK1, aqueous$logK2, tolerance = 0.1)
})
