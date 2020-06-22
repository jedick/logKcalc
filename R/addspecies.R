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
addspecies <- function(LOGK, ispecies) {
  # Initialize output
  init <- list(n = 0, lines = character())
  # aqueous, mineral, gas: the ones we can add
  # redox, electron, oxide: will not be changed
  ADDS <- list(redox = init, aqueous = init, electron = init, mineral = init, gas = init, oxide = init)
  if(is.null(ispecies)) return(ADDS)
  if(length(ispecies) == 0) return(ADDS)
  # Set basis species
  CHNOSZ::basis(LOGK$basis$map$CHNOSZ)
  on.exit(CHNOSZ::basis(delete = TRUE))
  # Load species
  CHNOSZ::species(ispecies)
  # Calculate affinity for all species
  a <- suppressMessages(CHNOSZ::affinity(T = LOGK$T, P = LOGK$P))
  # Loop over species to add them to the output
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
      Z <- suppressMessages(CHNOSZ::info(ispecies[i], check.it = FALSE)$Z)
      Z <- sprintf("%3.0f", Z)
      # default ion size from UNITHERM
      # TODO: allow user to change ion size
      r <- 4.5
      r <- sprintf("%5.1f", r)
      head1 <- paste0("     charge=", Z, "      ion size=", r, " A      mole wt.=", mw, " g")
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
    message("Adding ", type, " species: ", name)
    # Get reaction stoichiometry
    stoich <- a$species[i, 1:nrow(a$basis)]
    # Use GWB names in reaction
    colnames(stoich) <- LOGK$basis$map$GWB
    # Get reaction lines
    rxn <- rxnlines(stoich)
    # Get logK lines
    # NOTE: logK of dissociation reaction is opposite that of formation reaction
    logK1 <- formatline(-a$values[[i]], 1)
    logK2 <- formatline(-a$values[[i]], 2)
    # Put together the lines for this species (including a blank line at the end)
    alllines <- c(name, head1, head2, rxn, logK1, logK2, "")
    # Add the lines to the output
    ADDS[[type]]$n <- ADDS[[type]]$n + 1
    ADDS[[type]]$lines <- c(ADDS[[type]]$lines, alllines)
  }
  ADDS
}
