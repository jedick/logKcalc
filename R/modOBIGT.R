# logKcalc/modOBIGT.R

# modify OBIGT for converting GWB files 20200614
modOBIGT <- function(mod) {
  # 'mod' is one or more keywords to identify a specific modification
  if(length(mod) == 0) {
    stop("missing 'mod' argument")
  } else if(length(mod) > 1) {
    invisible(lapply(mod, modOBIGT))
  } else if(mod=="addSUPCRT") {
    # adds minerals from SUPCRT92 that are not in OBIGT
    #CHNOSZ::add.obigt("SUPCRT92", force = FALSE)
    # 20200614 that should work, but doesn't because of addition of some cr2 phases to existing Berman cr entries, so subcrt("albite") produces
    #Error in berman(PAR$name, T = T, P = P, thisinfo = PAR) : Data for ALBITE not available.
    supfile <- system.file("extdata/OBIGT/SUPCRT92.csv", package = "CHNOSZ")
    supdat <- utils::read.csv(supfile, as.is = TRUE)
    newnames <- setdiff(supdat$name, CHNOSZ::thermo()$obigt$name)
    CHNOSZ::add.obigt("SUPCRT92", newnames)
  } else if(mod=="allSUPCRT") {
    # adds all minerals from SUPCRT92, possibly replacing Berman ones
    CHNOSZ::add.obigt("SUPCRT92", force = TRUE)
  } else if(mod=="noBerman") {
    # removes minerals that use the Berman equations
    message("modOBIGT: removing Berman minerals")
    OBIGT <- CHNOSZ::thermo()$obigt
    # adjust for earlier versions of CHNOSZ that didn't have the "units" column in OBIGT
    if(utils::packageVersion("CHNOSZ") < "1.3.3") iBerman <- OBIGT$state=="cr" & apply(is.na(OBIGT[, 8:20]), 1, all)
    else iBerman <- OBIGT$state=="cr" & apply(is.na(OBIGT[, 9:21]), 1, all)
    OBIGT <- OBIGT[!iBerman, ]
    CHNOSZ::thermo(obigt = OBIGT)
  } else if(mod=="steam") {
    steamfile <- system.file("extdata/steam.csv", package = "logKcalc")
    CHNOSZ::add.obigt(steamfile)
  } else if(mod=="antigorite/2") {
    # divide the properties of antigorite by 2 (for thermo.tdat) 20200529
    CHNOSZ::mod.obigt("antigorite", formula = "Mg24Si17O49.5H31O24")
    name <- NULL
    bdat <- subset(CHNOSZ::berman(), name == "antigorite")
    bdat[, 2:30] <- bdat[, 2:30] / 2
    bfile <- tempfile()
    utils::write.csv(bdat, bfile, row.names = FALSE, quote = FALSE)
    CHNOSZ::thermo("opt$Berman" = bfile)
  } else if(mod=="realgar*4") {
    # multiply the properties of realgar by 4 (for Unitherm) 20200614
    irealgar <- CHNOSZ::info("realgar,alpha")
    obigt <- CHNOSZ::thermo()$obigt
    obigt$formula[irealgar] <- "As4S4"
    if(utils::packageVersion("CHNOSZ") < "1.3.3") obigt[irealgar, 8:20] <- obigt[irealgar, 8:20] * 4
    else obigt[irealgar, 9:21] <- obigt[irealgar, 9:21] * 4
    CHNOSZ::thermo(obigt = obigt)
  } else stop("unrecognized 'mod' argument: ", mod)
}

