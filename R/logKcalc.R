# logKcalc/logKcalc.R
# - find reactions in a GWB file
# - calculate logK at specified T and P using OBIGT database
# - save calculated values of logK in modified GWB file
# 20200422 jmd first version - calculate logK for redox species
# 20200429 output logK for redox species
# 20200503 calculate and output logK for aqueous species
# 20200523 split out functions readhead, calclogK
# 20200524 split out functions mapnames, writedat
# 20200525 new functions readlogK and logKcomp

## for debugging:
#source("mapnames.R")
#source("readhead.R")
#source("calclogK.R")
#source("writedat.R")

logKcalc <- function(infile = "thermo.tdat", outfile = "thermo_OBIGT.tdat", T = NULL, P = "Psat", maxprint = Inf) {
## for debugging:
#  infile <- "thermo.tdat"
#  outfile <- "thermo_OBIGT.tdat"
#  T <- NULL
#  P <- "Psat"
  maxprint <- Inf
  if(!file.exists(infile)) stop("file specified by 'infile' doesn't exist")
  # read the GWB data file
  LINES <- readLines(infile)
  # get the header data
  HEAD <- readhead(LINES)
  # calculate the logK values for available species
  LOGK <- calclogK(LINES, HEAD, T = T, P = P, maxprint = maxprint)
  # write the new file
  writedat(outfile, LINES, HEAD, LOGK)
}
