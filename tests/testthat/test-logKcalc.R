context("logKcalc")

test_that("Processing a K2GWB file works as expected", {
  # process the K2GWB file
  infile <- system.file("extdata/ThermoGWB_15_6_2020.tdat", package = "logKcalc")
  outfile <- file.path(tempdir(), "ThermoGWB_OBIGT.tdat")
  reset()
  modOBIGT(c("addSUPCRT", "steam", "realgar*4"))
  logKcalc(infile, outfile)
  # the reference file
  reffile <- system.file("extdata/tests/ThermoGWB_OBIGT.tdat", package = "logKcalc")
  # instead of testing for identical files,
  # just check that the values are tolerably close to each other
  LOGK <- logKcomp(reffile, outfile, type = "mineral", plot.it = FALSE)
  expect_equal(LOGK$logK1, LOGK$logK2, tolerance = 0.1)
})

test_that("Processing a GWB file works as expected", {
  # process the thermo_24elements.tdat file
  infile <- system.file("extdata/thermo_24elements.tdat", package = "logKcalc")
  outfile <- file.path(tempdir(), "thermo_OBIGT.tdat")
  reset()
  modOBIGT(c("addSUPCRT", "steam", "antigorite/2"))
  addOBIGT("As(OH)4-")
  addOBIGT("Sn++++")
  logKcalc(infile, outfile)
  reffile <- system.file("extdata/tests/thermo_OBIGT.tdat", package = "logKcalc")
  LOGK <- logKcomp(reffile, outfile, type = "mineral", plot.it = FALSE)
  expect_equal(LOGK$logK1, LOGK$logK2, tolerance = 0.1)
})
