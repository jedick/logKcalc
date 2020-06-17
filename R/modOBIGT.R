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

# (v1) Calculate thermodynamic parameters of redox or aqueous species in a GWB file 20200615
# (v0) Calculate thermodynamic parameters of As(OH)4- from logK values 20200614
addOBIGT <- function(species, formula = NULL, file = system.file("extdata/thermo_24swapped.tdat", package = "logKcalc"), tolerance = 0.05) {

  ## Figure out the chemical formula of the species
  if(is.null(formula)) {
    formula <- mapnames(species, return.processed.name = TRUE)
  }

  ## Get T, logK and reaction coefficients from file
  LINES <- readLines(file)
  HEAD <- readhead(LINES, quiet = TRUE)
  T <- HEAD$T
  LOGK <- readlogK(file, quiet = TRUE)
  # Look for the species in the redox or aqueous species
  iredox <- match(species, LINES[HEAD$ispecies$redox])
  iaqueous <- match(species, LINES[HEAD$ispecies$aqueous])
  if(!is.na(iredox)) {
    logK <- LOGK$redox[[iredox]]
    rxn <- readrxn(LINES, HEAD$ihead, HEAD$ispecies$redox[iredox])
  } else if(!is.na(iaqueous)) {
    logK <- LOGK$aqueous[[iaqueous]]
    rxn <- readrxn(LINES, HEAD$ihead, HEAD$ispecies$aqueous[iaqueous])
  } else stop("species ", species, " is not a redox or aqueous species in ", file)
  # Change "500" to NA
  logK[logK == 500] <- NA

  ## Get Gibbs energy of species from logK of reaction
  # Calculate T in Kelvin
  TK <- CHNOSZ::convert(T, "K")
  # logK gives input values for ΔG°r of the reaction
  Grin <- CHNOSZ::convert(logK, "G", TK)
  # Calculate ΔG°r of the dissociation reaction without the reactant species
  rxnspecies <- mapnames(rxn$species)$CHNOSZ
  rxncoeff <- rxn$coeff
  Grout <- suppressWarnings(suppressMessages(CHNOSZ::subcrt(rxnspecies, rxncoeff, T = T)$out$G))
  # Calculate ΔG°f of the reactant species
  Gf <- Grout - Grin

  ## Solve for G, S, and Cp
  # Make an 'lm' model object for given Cp
  Gfun <- function(Cp = 0) {
    Tr <- 298.15
    TTr <- TK - Tr
    # Subtract Cp term from Gf
    GCp <- Cp * (TK - Tr - TK * log(TK / Tr))
    GCp[is.na(GCp)] <- 0
    GfCp <- Gf - GCp
    # Write linear model in Ttr -- slope is -S
    stats::lm(GfCp ~ TTr)
  }
  # Calculate the sum of squares of residuals for given Cp
  Sqfun <- function(Cp) sum(Gfun(Cp)$residuals ^ 2)
  # Find the Cp with minimum sum of squares of residuals
  Cp <- stats::optimize(Sqfun, c(-100, 100))$minimum
  # Calculate the fitted G and S for this Cp
  G <- Gfun(Cp)$coefficients[[1]]
  S <- - Gfun(Cp)$coefficients[[2]]

  ## Make a new species
  # Use ref1 to store the Tmax
  Tmax <- max(T[!is.na(logK)])
  CHNOSZ::mod.obigt(species, formula = formula, state = "aq", ref1 = Tmax, G = G, S = S, c1 = Cp)
  # We need to call mod.obigt() a second time to set Z = 0 (to avoid triggering HKF omega derivatives)
  suppressMessages(CHNOSZ::mod.obigt(species, state = "aq", z = 0))
  # Calculate ΔG°r of the reaction with the reactant species
  rxnspecies <- c(species, rxnspecies)
  rxncoeff <- c(-1, rxncoeff)
  logKcalc <- suppressMessages(CHNOSZ::subcrt(rxnspecies, rxncoeff, T = T)$out$logK)
  # Re-read Tmax from the database and set the appropriate values of logK to NA
  Tmax <- as.numeric(suppressMessages(CHNOSZ::info(CHNOSZ::info(species))$ref1))
  logKcalc[T > Tmax] <- NA
  # Check that calculated values are close to input values
  stopifnot(all.equal(logK, logKcalc, tolerance = tolerance, scale = 1))

}
