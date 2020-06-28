# logKcalc/logKcalc.R
# - find reactions in a GWB file
# - calculate logK at specified T and P using OBIGT database
# - save calculated values of logK in modified GWB file
# 20200422 jmd first version - calculate logK for redox species
# 20200429 output logK for redox species
# 20200503 calculate and output logK for aqueous species
# 2020052u3 split out functions readhead, calclogK
# 20200524 split out functions mapnames, writedat
# 20200525 new functions readlogK and logKcomp

## for debugging:
#source("utils.R")
#source("mapnames.R")
#source("readhead.R")
#source("calclogK.R")
#source("addspecies.R")
#source("nonideal.R")
#source("writedat.R")

logKcalc <- function(infile = "thermo.tdat", outfile = "thermo_OBIGT.tdat",
  T = NULL, P = "Psat", ispecies = NULL, a0_ion = NULL, a0_neutral = 0,
  DH.method = "bdot", maxprint = Inf) {
  # set defaults for a0_ion
  if(missing(a0_ion)) {
    if(identical(DH.method, "bdot")) a0_ion <- 4.5
    else if(identical(DH.method, "bgamma")) a0_ion <- 3.72
    else stop("DH.method should be 'bdot' or 'bgamma'")
  }
  if(!file.exists(infile)) stop("file specified by 'infile' doesn't exist")
  # read the GWB data file
  LINES <- readLines(infile)
  LINES <- cleanUTF8(LINES, infile)
  # get the header data
  HEAD <- readhead(LINES)
  # calculate the logK values for available species
  LOGK <- calclogK(LINES, HEAD, T, P, maxprint)
  # get lines for added species
  ADDS <- addspecies(LOGK, ispecies, a0_ion, a0_neutral, DH.method)
  # write the new file
  if(DH.method=="bdot") writedat(outfile, LINES, HEAD, LOGK, ADDS, infile, DH.method)
  if(DH.method=="bgamma") writedat(outfile, LINES, HEAD, LOGK, ADDS, infile, DH.method, a0_ion)
}
