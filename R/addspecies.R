# logKcalc/addspecies.R
# Create entries for one or more new species 20200619

# Function to create reaction lines 20200620
rxnlines <- function(stoich) {
  # remove species with 0 coefficient
  stoich <- stoich[, stoich != 0]
  # create reaction header line
  rhead <- paste0("     ", ncol(stoich), " species in reaction")
  # intialize reaction species lines
  rlines <- character()
  k <- 0
  # loop over species
  for(j in 1:ncol(stoich)) {
    # start a newline for every three species
    if((j-1) %% 3 == 0) {
      k <- k + 1
      rlines[k] <- "     "
    }
    # add coefficient and species name
    coeff <- sprintf("%8.3f", stoich[, j])
    sname <- sprintf("%-12s", names(stoich)[j])
    rlines[k] <- paste0(rlines[k], coeff, " ", sname, " ")
  }
  # return header and species lines
  c(rhead, rlines)
}

# ispecies: species index in thermo()$obigt
addspecies <- function(LOGK, ispecies, a0_ion, a0_neutral, DH.method) {
  # set defaults for a0_ion
  if(is.null(a0_ion)) {
    if(identical(DH.method, "bdot")) a0_ion <- 4.5
    else if(identical(DH.method, "bgamma")) a0_ion <- 3.72
    else stop("DH.method should be 'bdot' or 'bgamma'")
  }
  # Initialize output
  init <- list(n = 0, lines = character(), refs = character())
  # aqueous, mineral, gas: the ones we can add
  # redox, electron, oxide: will not be changed
  ADDS <- list(redox = init, aqueous = init, electron = init, mineral = init, gas = init, oxide = init)
  if(is.null(ispecies)) return(ADDS)
  if(length(ispecies) == 0) return(ADDS)
  # Set basis species (default logact is 0)
  CHNOSZ::basis(LOGK$basis$map$OBIGT)
  on.exit(CHNOSZ::basis(delete = TRUE))
  # Load species (set logact to 0!)
  CHNOSZ::species(ispecies, 0)
  # Calculate affinity for all species
  a <- suppressMessages(CHNOSZ::affinity(T = LOGK$T, P = LOGK$P))
  # Make a0_ion and a0_species the same length as the number of species
  a0_ion <- rep(a0_ion, length.out = length(ispecies))
  a0_neutral <- rep(a0_neutral, length.out = length(ispecies))
  # Loop over species to add them to the output
  extramsg <- character()
  for(i in 1:length(ispecies)) {
    # Get species name and state
    name <- a$species$name[i]
    state <- a$species$state[i]
    # Create header lines depending on state
    # TODO: use element mole wt. from GWB file
    mw <- CHNOSZ::mass(ispecies[i])
    mw <- sprintf("%10.4f", mw)
    if(state == "aq") {
      type <- "aqueous"
      # first header line: charge, ion size, and molecular weight
      #Z <- suppressMessages(CHNOSZ::info(ispecies[i], check.it = FALSE)$Z)
      # use the charge from the formula, not the "Z" parameter, because
      # the latter is set to zero to disable the g-function for some organic ions 20200622
      Z <- CHNOSZ::makeup(ispecies[i])["Z"]
      if(is.na(Z)) Z <- 0
      # get ion size parameter for ions or aqueous species
      if(Z == 0) a0 <- a0_neutral[i] else a0 <- a0_ion[i]
      extramsg <- paste(" with a0 =", a0)
      Z <- sprintf("%3.0f", Z)
      # add a decimal place for a0 in bgamma (default 3.72) 20200623
      if(DH.method=="bdot") a0 <- sprintf("%5.1f", a0)
      if(DH.method=="bgamma") a0 <- sprintf("%5.2f", a0)
      head1 <- paste0("     charge=", Z, "      ion size=", a0, " A      mole wt.=", mw, " g")
      # there's no second header line
      head2 <- character()
    } else if(state == "cr") {
      type <- "mineral"
      # first header line: formula
      formula <- CHNOSZ::info(ispecies[i], check.it = FALSE)$formula
      head1 <- paste("     formula=", formula)
      # second header line: volume and molecular weight
      V <- CHNOSZ::info(ispecies[i], check.it = FALSE)$V
      V <- sprintf("%10.4f", V)
      head2 <- paste0("     mole vol.=", V, " cc      mole wt.=", mw, " g")
    } else if(state == "gas") {
      type <- "gas"
      # first header line: molecular weight
      head1 <- paste0("     mole wt.=", mw, " g")
      head2 <- character()
    } else stop(paste("can't handle", state, "species"))
    # Don't allow duplicated species -- check both existing GWB and OBIGT names
    existing <- c(names(LOGK[[type]]$logKs), LOGK[[type]]$speciesOBIGT)
    if(name %in% existing) stop(paste(name, "is a duplicated", type, "species"))
    message("Adding ", type, " species ", name, extramsg)
    # Get reaction stoichiometry
    stoich <- a$species[i, 1:nrow(a$basis)]
    # Use GWB names in reaction
    colnames(stoich) <- LOGK$basis$map$GWB
    # Get reaction lines
    rxn <- rxnlines(stoich)
    # Get logK lines
    # NOTE: logK of dissociation reaction is opposite that of formation reaction
    logK1 <- formatline(-a$values[[i]], 1, na.500 = TRUE)
    logK2 <- formatline(-a$values[[i]], 2, na.500 = TRUE)
    # Make a reference line 20200623
    refline <- "* no references available"
    r1 <- CHNOSZ::info(ispecies[i], check.it = FALSE)$ref1
    # remove suffixes
    r1 <- sapply(strsplit(sapply(strsplit(r1, "\\."), "[", 1), " "), "[", 1)
    if(!is.na(r1)) {
      ADDS[[type]]$refs <- c(ADDS[[type]]$refs, r1)
      refline <- paste0("* [", r1, "]")
      r2 <- CHNOSZ::info(ispecies[i], check.it = FALSE)$ref2
      r2 <- sapply(strsplit(sapply(strsplit(r2, "\\."), "[", 1), " "), "[", 1)
      if(!is.na(r2)) {
        ADDS[[type]]$refs <- c(ADDS[[type]]$refs, r2)
        refline <- paste0("* [", r1, ", ", r2, "]")
      }
    }
    # Put together the lines for this species (including a blank line at the end)
    alllines <- c(name, head1, head2, rxn, logK1, logK2, refline, "")
    # Add the lines to the output
    ADDS[[type]]$n <- ADDS[[type]]$n + 1
    ADDS[[type]]$lines <- c(ADDS[[type]]$lines, alllines)
  }
  ADDS
}
