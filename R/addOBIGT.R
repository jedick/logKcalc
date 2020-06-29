# logKcalc/addOBIGT.R
# add species to OBIGT using G, S, and Cp
# fit to logK values from a GWB file 20200615

### Taken from R: src/library/base/demo/error.catching.R 20200624
##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings (with value) and errors
##' @param expr an \R expression to evaluate
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler;
##' Copyright (C) 2010-2012  The R Core Team
tryCatch.W.E <- function(expr)
{
    W <- NULL
    w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
				     warning = w.handler),
	 warning = W)
}

# (v2) Add minerals 20200624
# (v1) Calculate thermodynamic parameters of redox or aqueous species in a GWB file 20200615
# (v0) Calculate thermodynamic parameters of As(OH)4- from logK values 20200614
addOBIGT <- function(species, formula = NULL, file = system.file("extdata/thermo_24swapped.tdat", package = "logKcalc"), tolerance = 0.1) {

  ## Figure out the chemical formula using heuristics for aqueous species formulas
  if(is.null(formula)) {
    # FIXME: this doesn't work for named species (e.g. minerals)
    formula <- mapnames(species, return.processed.name = TRUE)
  }
  # Check if the formula is valid
  trymakeup <- tryCatch.W.E(CHNOSZ::makeup(formula))
  if(!is.null(trymakeup$warning) | inherits(trymakeup$value, "error"))
    stop(paste0("'", formula, "' isn't a valid chemical formula; please provide one in the 'formula' argument"))

  ## Get T, logK and reaction coefficients from file
  LINES <- readLines(file)
  HEAD <- readhead(LINES, quiet = TRUE)
  T <- HEAD$T
  LOGK <- readlogK(file, quiet = TRUE)
  # Look for the species in the redox or aqueous species
  iredox <- match(species, line2name(LINES[HEAD$ispecies$redox]))
  iaqueous <- match(species, line2name(LINES[HEAD$ispecies$aqueous]))
  imineral <- match(species, line2name(LINES[HEAD$ispecies$mineral]))
  # A flag to indicate if the species is aqueous (in redox or aqueous species blocks)
  isaq <- TRUE
  # Get the logK values and reaction coefficients
  if(!is.na(iredox)) {
    logK <- LOGK$redox[[iredox]]
    rxn <- readrxn(LINES, HEAD$ihead, HEAD$ispecies$redox[iredox])
  } else if(!is.na(iaqueous)) {
    logK <- LOGK$aqueous[[iaqueous]]
    rxn <- readrxn(LINES, HEAD$ihead, HEAD$ispecies$aqueous[iaqueous])
  } else if(!is.na(imineral)) {
    logK <- LOGK$mineral[[imineral]]
    rxn <- readrxn(LINES, HEAD$ihead, HEAD$ispecies$mineral[imineral])
    isaq <- FALSE
  } else stop(species, " is not a redox, aqueous or mineral species in ", file)
  # Change "500" to NA
  logK[logK == 500] <- NA

  ## Get Gibbs energy of species from logK of reaction
  # Calculate T in Kelvin
  TK <- CHNOSZ::convert(T, "K")
  # logK gives input values for ΔG°r of the reaction
  Grin <- CHNOSZ::convert(logK, "G", TK)
  # Calculate ΔG°r of the dissociation reaction without the reactant species
  rxnspecies <- mapnames(rxn$species)$OBIGT
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
  # Get the highest temperature for the fitted logK values
  Tmax <- max(T[!is.na(logK)])
  if(isaq) {
    # Use the formula as the name (because CHNOSZ::expr.species("AsO4---") produces an error) 20200624
    species <- formula
    # For aqueous species, use abbrv to store the Tmax
    # Explicitly set all parameters in case this species already exists in OBIGT 20200624
    moargs <- list(name = species, formula = formula, state = "aq", ref1 = "logK_fit", ref2 = basename(file),
      G = G, S = S, c1 = Cp, abbrv = Tmax,
      H = NA, V = NA, a1 = NA, a2 = NA, a3 = NA, a4 = NA, c2 = NA, omega = NA
    )
    if(utils::packageVersion("CHNOSZ") >= "1.3.3") moargs <- c(moargs, list(E_units = "cal"))
    do.call(CHNOSZ::mod.obigt, moargs)
    # We need to call mod.obigt() a second time to set Z = 0 (to avoid triggering HKF omega derivatives)
    suppressMessages(CHNOSZ::mod.obigt(species, state = "aq", z = 0))
  } else {
    # Make a mineral 20200624
    # Nudge the Tmax to allow calculation at exactly that temperture 20200625
    if(utils::packageVersion("CHNOSZ") < "1.3.4") Tmax <- Tmax + 0.01
    moargs <- list(species, formula = formula, state = "cr", ref1 = "logK_fit", ref2 = basename(file),
      G = G, S = S, Cp = Cp, T = CHNOSZ::convert(Tmax, "K"),
      H = NA, V = NA, a = NA, b = NA, c = NA, d = NA, e = NA, f = NA, lambda = NA)
    if(utils::packageVersion("CHNOSZ") >= "1.3.3") moargs <- c(moargs, list(E_units = "cal"))
    do.call(CHNOSZ::mod.obigt, moargs)
  }
  # Calculate ΔG°r of the reaction with the reactant species
  rxnspecies <- c(species, rxnspecies)
  rxncoeff <- c(-1, rxncoeff)
  logKcalc <- suppressMessages(CHNOSZ::subcrt(rxnspecies, rxncoeff, T = T)$out$logK)
  if(isaq) {
    # Re-read Tmax from the database and set the appropriate values of logK to NA
    Tmax <- as.numeric(suppressMessages(CHNOSZ::info(CHNOSZ::info(species))$abbrv))
    logKcalc[T > Tmax] <- NA
  }
  # Check that calculated values are close to input values
  stopifnot(all.equal(logK, logKcalc, tolerance = tolerance, scale = 1))
}
