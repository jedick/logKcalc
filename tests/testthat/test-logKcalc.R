context("logKcalc")

#test_that("Processing a K2GWB file works as expected", {
#  # process the K2GWB file
#  infile <- system.file("extdata/ThermoGWB_15_6_2020.tdat", package = "logKcalc")
#  outfile <- file.path(tempdir(), "ThermoGWB_OBIGT.tdat")
#  reset()
#  modOBIGT(c("addSUPCRT", "steam", "realgar*4"))
#  logKcalc(infile, outfile)
#  # the reference file
#  reffile <- system.file("extdata/tests/ThermoGWB_OBIGT.tdat", package = "logKcalc")
#  # instead of testing for identical files,
#  # just check that the values are tolerably close to each other
#  LOGK <- logKcomp(reffile, outfile, type = "mineral", plot.it = FALSE)
#  expect_equal(LOGK$logK1, LOGK$logK2, tolerance = 0.1)
#})

test_that("Modifying the database and processing a GWB file work as expected", {
  # process the thermo_24elements.tdat file
  infile <- system.file("extdata/thermo_24elements.tdat", package = "logKcalc")
  outfile <- file.path(tempdir(), "thermo_OBIGT.tdat")
  reset()
  modOBIGT(c("addSUPCRT", "steam", "antigorite/2"))
  addOBIGT("As(OH)4-")
  addOBIGT("Sn++++")
  logKcalc(infile, outfile)
  reffile <- system.file("extdata/tests/thermo_24OBIGT.tdat", package = "logKcalc")
  # Testing for identical files is problematic with ongoing updates to OBIGT,
  # so just check that the values are close to each other
  mineral <- logKcomp(reffile, outfile, type = "mineral", plot.it = FALSE)
  # The formula for scorodite was corrected in CHNOSZ_1.3.6-9
  mineral <- mineral[rownames(mineral)!="Scorodite", ]
  expect_equal(mineral$logK1, mineral$logK2, tolerance = 0.1)
  aqueous <- logKcomp(reffile, outfile, type = "aqueous", plot.it = FALSE)
  expect_equal(aqueous$logK1, aqueous$logK2, tolerance = 0.1)
})

test_that("Adding duplicate species produces an error", {
  infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
  iAuCl2 <- info("AuCl2-")
  expect_error(logKcalc(infile, ispecies = iAuCl2), "duplicated aqueous species")
})

test_that("Adding species to the output works as expected", {
  # process the thermo_12elements.tdat file
  infile <- system.file("extdata/thermo_12elements.tdat", package = "logKcalc")
  outfile <- file.path(tempdir(), "thermo_12OBIGT.tdat")
  reset()
  modOBIGT(c("addSUPCRT", "steam"))
  # add a selection of aqueous, mineral and gas species
  ispecies <- info(c("AuCl", "pyrrhotite", "ammonia"))
  logKcalc(infile, outfile, ispecies = ispecies)
  reffile <- system.file("extdata/tests/thermo_12OBIGT.tdat", package = "logKcalc")
  reflines <- readLines(reffile)
  outlines <- readLines(outfile)
  # normalize “ and ” to "
  # (straight quote is created in R CMD check -- locale setting??)
  reflines <- gsub('“', '"', reflines)
  reflines <- gsub('”', '"', reflines)
  outlines <- gsub('“', '"', outlines)
  outlines <- gsub('”', '"', outlines)
  expect_identical(outlines, reflines)
})
