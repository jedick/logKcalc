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
    #add.OBIGT("SUPCRT92", force = FALSE)
    # 20200614 that should work, but doesn't because of addition of some cr2 phases to existing Berman cr entries, so subcrt("albite") produces
    #Error in Berman(PAR$name, T = T, P = P, thisinfo = PAR) : Data for ALBITE not available.
    supfile <- system.file("extdata/OBIGT/SUPCRT92.csv", package = "CHNOSZ")
    supdat <- utils::read.csv(supfile, as.is = TRUE)
    newnames <- setdiff(supdat$name, CHNOSZ::thermo()$OBIGT$name)
    add.OBIGT("SUPCRT92", newnames)
  } else if(mod=="allSUPCRT") {
    # adds all minerals from SUPCRT92, possibly replacing Berman ones
    add.OBIGT("SUPCRT92", force = TRUE)
  } else if(mod=="noBerman") {
    # removes minerals that use the Berman equations
    message("modOBIGT: removing Berman minerals")
    OBIGT <- CHNOSZ::thermo()$OBIGT
    # Adjust for addition of columns to OBIGT in versions of CHNOSZ
    # ("units" column in 1.3.3, "model" column in 2.0.0)
    iBerman <- OBIGT$state=="cr" & apply(is.na(OBIGT[, 10:22]), 1, all)
    OBIGT <- OBIGT[!iBerman, ]
    CHNOSZ::thermo(OBIGT = OBIGT)
  } else if(mod=="steam") {
    steamfile <- system.file("extdata/steam.csv", package = "logKcalc")
    add.OBIGT(steamfile)
  } else if(mod=="antigorite/2") {
    # divide the properties of antigorite by 2 (for thermo.tdat) 20200529
    mod.OBIGT("antigorite", formula = "Mg24Si17O49.5H31O24")
    name <- NULL
    bdat <- subset(Berman(), name == "antigorite")
    bdat[, 2:30] <- bdat[, 2:30] / 2
    bfile <- tempfile()
    utils::write.csv(bdat, bfile, row.names = FALSE, quote = FALSE)
    CHNOSZ::thermo("opt$Berman" = bfile)
  } else if(mod=="realgar*4") {
    # multiply the properties of realgar by 4 (for Unitherm) 20200614
    irealgar <- CHNOSZ::info("realgar,alpha")
    OBIGT <- CHNOSZ::thermo()$OBIGT
    OBIGT$formula[irealgar] <- "As4S4"
    OBIGT[irealgar, 10:22] <- OBIGT[irealgar, 10:22] * 4
    CHNOSZ::thermo(OBIGT = OBIGT)
  } else stop("unrecognized 'mod' argument: ", mod)
}

