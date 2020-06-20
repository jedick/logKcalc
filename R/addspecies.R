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
  ADDS <- list(aqueous = init, mineral = init, gas = init)
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
    if(state == "aq") {
      # get and format charge, ion size, and molecular weight
      Z <- suppressMessages(CHNOSZ::info(ispecies[i])$Z)
      Z <- sprintf("%3.0f", Z)
      # TODO: lookup or have user specify ion size
      r <- 4.0
      r <- sprintf("%5.1f", r)
      # TODO: use element mole wt. from GWB file
      mw <- CHNOSZ::mass(ispecies[i])
      mw <- sprintf("%10.4f", mw)
      head <- paste0("     charge=", Z, "      ion size=", r, " A      mole wt.=", mw, " g")
    }
    # Get reaction stoichiometry
    stoich <- a$species[i, 1:nrow(a$basis)]
    # Use GWB names in reaction
    colnames(stoich) <- LOGK$basis$map$GWB
    # Get reaction lines
    rxn <- rxnlines(stoich)
    # Get logK lines
    logK1 <- formatline(a$values[[i]], 1)
    logK2 <- formatline(a$values[[i]], 2)
    # Put together the lines for this species (including a blank line at the end)
    alllines <- c(name, head, rxn, logK1, logK2, "")
    # Add the lines to the output
    if(state=="aq") type <- "aqueous"
    else if(state=="cr") type <- "mineral"
    else if(state=="gas") type <- "gas"
    else stop(paste("can't handle", state, "species"))
    ADDS[[type]]$n <- ADDS[[type]]$n + 1
    ADDS[[type]]$lines <- c(ADDS[[type]]$lines, alllines)
  }
  ADDS
}
