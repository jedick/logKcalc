# calclogK.R
# calculate logK values for inclusion in data file for GWB
# 20200525 jmd initial version

# a function to read the *reaction* coefficients and species
readrxn <- function(LINES, ihead, i) {
  # find the next header line after this one
  # skip at least one line for "*    formula=" after the name 20200611
  j <- min(ihead[ihead > i + 1])
  # get all lines between lines i and j
  k <- (i + 1):(j - 1)
  thisdat <- LINES[k]
  # find the line describing the reaction
  jrxn <- grep("in\\ reaction", thisdat)
  # get the number of species in the reaction
  desc <- strsplit(thisdat[jrxn], " ")[[1]]
  desc <- desc[! desc == ""]
  nspecies <- as.numeric(desc[1])
  # calculate the lines of data for the reaction
  nlines <- ceiling(nspecies / 3)
  irxn <- jrxn + 1:nlines
  # parse the reaction data
  rxn <- thisdat[irxn]
  rxn <- strsplit(rxn, " ")
  rxn <- lapply(rxn, function(x) x[nzchar(x)])
  rxn <- unlist(rxn)
  # odd values are coefficients, even are species
  coeff <- as.numeric(rxn[c(TRUE, FALSE)])
  species <- rxn[c(FALSE, TRUE)]
  list(coeff = coeff, species = species)
}

# a function to print messages about unavailable species,
# and possibly the species names themselves 20200615
printNA <- function(start, end, type, species, maxprint) {
  if(length(species) > maxprint) {
    message(start, " ", length(species), " ", type, " ", end, ".")
  } else {
    message(start, " these ", type, " ", end, ":")
    print(paste(species, collapse = ", "))
  }
}

calclogK <- function(LINES, HEAD, T = NULL, P = "Psat", maxprint = Inf) {
  # process T and P arguments
  if(is.null(T)) T <- HEAD$T
  message("The temperatures for the logK calculation are:")
  print(paste(T, collapse = ", "))
  if(is.null(P)) P <- HEAD$P else {
    if(identical(P, "Psat")) {
      # get rounded values and nudge last decimal place to ensure we are in the liquid region
      P <- round(CHNOSZ::water("Psat", T = CHNOSZ::convert(T, "K"))$Psat, 4)
      P[P > 1] <- P[P > 1] + 0.0001
    }
  }
  message("The pressures for the logK calculation are:")
  print(paste(P, collapse = ", "))
  # make a list to hold information on the different types of species
  OUT <- list(T = T, P = P, basismap = NA)
  # get information on the different types of species
  for(type in c("redox", "aqueous", "electron", "mineral", "gas", "oxide")) {
    # get the line numbers with the species names
    inames <- HEAD$ispecies[[type]]
    # read the species names
    # and exclude type= ... at end 20200526
    speciesGWB <- line2name(LINES[inames])
    # get the reactions for each species
    rxnGWB <- lapply(inames, function(i) readrxn(LINES, HEAD$ihead, i))
    # test if all basis species or redox species are available
    if(type=="redox") {
      # set the basis species
      OUT$basismap <- mapnames(LINES[HEAD$ispecies$basis], type = "basis", na.omit = TRUE)
      basis <- CHNOSZ::basis(OUT$basismap$CHNOSZ)
      message(paste("The basis has", nrow(basis), "elements:"))
      print(paste(colnames(basis)[1:nrow(basis)], collapse=", "))
      inbasis <- sapply(rxnGWB, function(x) all(x$species %in% OUT$basismap$GWB))
      if(!all(inbasis)) {
        printNA("The basis species for", "species are not available", type, speciesGWB[!inbasis], maxprint)
      }
      # remove the basis species so that incorrect reactions for other
      # types of species are not automatically balanced 20200611
      CHNOSZ::basis(delete = TRUE)
    } else {
      # add O2(g) here, needed for the free electron in thermo.tdat 20200526
      inbasis <- sapply(rxnGWB, function(x) all(x$species %in% c(OUT$basismap$GWB, names(OUT$redox$logKs), "O2(g)")))
      if(!all(inbasis)) {
        printNA("The basis or redox species for", "species are not available", type, speciesGWB[!inbasis], maxprint)
      }
    }
    # skip types with no species 20200609
    if(length(inames)==0) {
      # don't forget to remove basis species 20200614
      if(type == "redox") basis(delete = TRUE)
      next
    }
    # remove species with missing basis species
    speciesGWB <- speciesGWB[inbasis]
    rxnGWB <- rxnGWB[inbasis]
    # map the names from GWB to CHNOSZ
    # use the 'type' argument here to restrict matches to a subset of OBIGT
    # so e.g. the mineral BaCrO4 doesn't match aqueous BaCrO4
    speciesmap <- lapply(speciesGWB, mapnames, type = type)
    speciesCHNOSZ <- unlist(lapply(speciesmap, "[", 2))
    speciesNA <- speciesGWB[is.na(speciesCHNOSZ)]
    if(length(speciesNA) > 0) {
      printNA("Can't map", "species to the OBIGT database", type, speciesNA, maxprint)
    }
    speciesGWB <- speciesGWB[!is.na(speciesCHNOSZ)]
    rxnGWB <- rxnGWB[!is.na(speciesCHNOSZ)]
    speciesCHNOSZ <- stats::na.omit(speciesCHNOSZ)
    if(length(speciesCHNOSZ) > 0) {
      print(paste("calculating logKs for", length(speciesCHNOSZ), type, "species ..."))
      logKs <- lapply(1:length(speciesCHNOSZ), function(i) {
        rxnCHNOSZ <- mapnames(rxnGWB[[i]]$species)$CHNOSZ
        coeff <- c(-1, rxnGWB[[i]]$coeff)
        species <- c(speciesCHNOSZ[i], rxnCHNOSZ)
        # print NA species names for debugging mapping problems
        if(any(is.na(species))) print(species)
        sres <- suppressMessages(CHNOSZ::subcrt(species, coeff, T = T, P = P, property = "logK"))
        logK <- sres$out$logK
        # get Tmax from ref1 (for species added by addOBIGT) 20200615
        Tmax <- stats::na.omit(suppressWarnings(as.numeric(suppressMessages(CHNOSZ::info(CHNOSZ::info(species))$ref1))))
        if(length(Tmax) > 0) logK[T > min(Tmax)] <- NA
        # if there are any warnings for unbalanced reactions, set logK to NA 20200614
        if(!is.null(sres$warnings)) {
          iswarn <- grepl("unbalanced", sres$warnings)
          if(any(iswarn)) {
            # parse the warning message to get the stoichiometry of the missing composition
            warn <- sres$warnings[iswarn]
            missing <- strsplit(warn, "missing ")[[1]][2]
            # tolerate a difference of up to 0.001 (for Realgar in thermo.dat)
            if(any(abs(round(CHNOSZ::makeup(missing), 3)) > 0.001)) logK[] <- NA
          }
        }
        logK
      })
      # get references 20200617
      ispecies <- suppressMessages(CHNOSZ::info(speciesCHNOSZ, check.it = FALSE))
      iinfo <- suppressMessages(CHNOSZ::info(ispecies, check.it = FALSE))
      ref1 <- iinfo$ref1
      ref2 <- iinfo$ref2
      # remove suffixes
      ref1 <- sapply(strsplit(sapply(strsplit(ref1, "\\."), "[", 1), " "), "[", 1)
      ref2 <- sapply(strsplit(sapply(strsplit(ref2, "\\."), "[", 1), " "), "[", 1)
      names(logKs) <- speciesGWB
    } else ref1 <- ref2 <- logKs <- NULL
    # save the information
    OUT[[type]] <- list(logKs = logKs, ref1 = ref1, ref2 = ref2)
  }
  OUT
}

