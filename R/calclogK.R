# calclogK.R
# calculate logK values for inclusion in data file for GWB
# 20200525 jmd initial version

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
  if(length(T) == 1) T <- rep(T, 8)
  if(length(T) != 8) stop("please supply 1 or 8 values for T")
  message("The temperatures for the logK calculation are:")
  print(paste(T, collapse = ", "))
  if(is.null(P)) P <- HEAD$P else {
    if(identical(P, "Psat")) {
      # get rounded values and nudge last decimal place to ensure we are in the liquid region
      P <- round(CHNOSZ::water("Psat", T = CHNOSZ::convert(T, "K"))$Psat, 4)
      P[P > 1] <- P[P > 1] + 0.0001
    }
  }
  if(length(P) == 1) P <- rep(P, 8)
  if(length(P) != 8) stop("please supply 1 or 8 values for P")
  message("The pressures for the logK calculation are:")
  print(paste(P, collapse = ", "))
  # make a list to hold information on the different types of species
  OUT <- list(T = T, P = P, basis = list())
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
      # first set the basis species
      OUT$basis$map <- mapnames(LINES[HEAD$ispecies$basis], type = "basis", na.omit = TRUE)
      # print list of elements without using CHNOSZ::basis() 20200701
      ispecies <- suppressMessages(CHNOSZ::info(OUT$basis$map$OBIGT))
      elements <- names(CHNOSZ::makeup(ispecies, sum = TRUE))
      message(paste("The", length(ispecies), "available basis species have", length(elements), "elements:"))
      print(paste(elements, collapse=", "))
      OUT$basis$nNA <- length(HEAD$ispecies$basis) - length(ispecies)
      # check if basis species are available for dissociation reactions for the redox species
      inbasis <- sapply(rxnGWB, function(x) all(x$species %in% OUT$basis$map$GWB))
      if(!all(inbasis)) {
        printNA("The basis species for", "species are not available", type, speciesGWB[!inbasis], maxprint)
      }
      # get references 20200618
      iinfo <- suppressMessages(CHNOSZ::info(ispecies, check.it = FALSE))
      ref1 <- iinfo$ref1
      ref2 <- iinfo$ref2
      # remove suffixes
      OUT$basis$ref1 <- sapply(strsplit(sapply(strsplit(ref1, "\\.[0-9]+"), "[", 1), " "), "[", 1)
      OUT$basis$ref2 <- sapply(strsplit(sapply(strsplit(ref2, "\\.[0-9]+"), "[", 1), " "), "[", 1)
      # get formulas 20200701
      OUT$basis$formula <- iinfo$formula
    } else {
      # add O2(g) here, needed for the free electron in thermo.tdat 20200526
      inbasis <- sapply(rxnGWB, function(x) all(x$species %in% c(OUT$basis$map$GWB, names(OUT$redox$logKs), "O2(g)")))
      if(!all(inbasis)) {
        printNA("The basis or redox species for", "species are not available", type, speciesGWB[!inbasis], maxprint)
      }
    }
    # skip types with no species 20200609
    if(length(inames)==0) next
    # remove species with missing basis species
    speciesGWB <- speciesGWB[inbasis]
    rxnGWB <- rxnGWB[inbasis]
    # skip calculating logKs for oxide species 20200622
    if(type == "oxide") {
      # https://stackoverflow.com/questions/58081243/turning-a-character-vector-into-an-empty-list-with-names-from-the-character-vect
      logKs <- Map(function(x) NULL, speciesGWB)
      nNA <- sum(!inbasis)
      OUT[[type]] <- list(logKs = logKs, nNA = nNA)
      next
    }
    # map the names from GWB to CHNOSZ
    # use the 'type' argument here to restrict matches to a subset of OBIGT
    # so e.g. the mineral BaCrO4 doesn't match aqueous BaCrO4
    speciesmap <- lapply(speciesGWB, mapnames, type = type)
    speciesOBIGT <- unlist(lapply(speciesmap, "[", 2))
    speciesNA <- speciesGWB[is.na(speciesOBIGT)]
    if(length(speciesNA) > 0) {
      printNA("Can't map", "species to the OBIGT database", type, speciesNA, maxprint)
    }
    speciesGWB <- speciesGWB[!is.na(speciesOBIGT)]
    rxnGWB <- rxnGWB[!is.na(speciesOBIGT)]
    speciesOBIGT <- as.character(stats::na.omit(speciesOBIGT))
    if(length(speciesOBIGT) > 0) {
      #print(paste("calculating logKs for", length(speciesOBIGT), type, "species ..."))
      logKs <- lapply(1:length(speciesOBIGT), function(i) {
        rxnCHNOSZ <- mapnames(rxnGWB[[i]]$species)$OBIGT
        coeff <- c(-1, rxnGWB[[i]]$coeff)
        species <- c(speciesOBIGT[i], rxnCHNOSZ)
        # print NA species names for debugging mapping problems
        if(any(is.na(species))) {
          message("Some species in the reaction are NA:")
          print(species)
        }
        sargs <- list(species = species, coeff, T = T, P = P, property = "logK")
        # Force "methane" to be gas in previous versions of CHNOSZ 20200625
        if(utils::packageVersion("CHNOSZ") <= "1.3.6") {
          if(species[1] == "methane") {
            # use species indices instead of names
            infsp <- suppressMessages(CHNOSZ::info(species))
            infsp[1] <- suppressMessages(CHNOSZ::info("methane", "gas"))
            sargs <- list(species = infsp, coeff, T = T, P = P, property = "logK")
          }
        }
        sres <- suppressWarnings(suppressMessages(do.call(CHNOSZ::subcrt, sargs)))
        logK <- sres$out$logK
        # get Tmax from abbrv (for species added by addOBIGT) 20200615
        Tmax <- stats::na.omit(suppressWarnings(as.numeric(suppressMessages(CHNOSZ::info(sres$reaction$ispecies, check.it = FALSE)$abbrv))))
        if(length(Tmax) > 0) logK[T > min(Tmax)] <- NA
        # if there are any warnings for unbalanced reactions, set logK to NA 20200614
        if(!is.null(sres$warnings)) {
          iswarn <- grepl("unbalanced", sres$warnings)
          if(any(iswarn)) {
            # parse the warning message to get the stoichiometry of the missing composition
            warn <- sres$warnings[iswarn]
            missing <- strsplit(warn, "missing ")[[1]][2]
            # tolerate a difference of up to 0.001 (for Realgar in thermo.dat)
            if(any(abs(round(CHNOSZ::makeup(missing), 3)) > 0.001)) {
              logK[] <- NA
              # show the warning as a message (capitalize first letter)
              message(gsub("^r", "R", warn))
            }
          }
        }
        logK
      })
      # remove species with all NA values (produced by unbalanced reactions) 20200621
      allisna <- sapply(logKs, function(x) all(is.na(x)))
      if(any(allisna)) {
        removed <- speciesGWB[allisna]
        message("Removing species with unavailable logKs:")
        print(paste(removed, collapse = ", "))
        speciesGWB <- speciesGWB[!allisna]
        speciesOBIGT <- speciesOBIGT[!allisna]
        logKs <- logKs[!allisna]
      }
      # get references 20200617
      ispecies <- suppressMessages(CHNOSZ::info(speciesOBIGT, check.it = FALSE))
      # force "methane" to be a gas in previous versions of CHNOSZ 20200625
      if(!utils::packageVersion("CHNOSZ") > "1.3.6") {
        imethane <- speciesOBIGT == "methane"
        if(any(imethane)) {
          ispecies[imethane] <- suppressMessages(CHNOSZ::info("methane", "gas", check.it = FALSE))
        }
      }
      iinfo <- suppressMessages(CHNOSZ::info(ispecies, check.it = FALSE))
      ref1 <- iinfo$ref1
      ref2 <- iinfo$ref2
      # remove suffixes
      ref1 <- sapply(strsplit(sapply(strsplit(ref1, "\\.[0-9]+"), "[", 1), " "), "[", 1)
      ref2 <- sapply(strsplit(sapply(strsplit(ref2, "\\.[0-9]+"), "[", 1), " "), "[", 1)
      names(logKs) <- speciesGWB
      # get formula 20200701
      formula <- iinfo$formula
    } else {
      ref1 <- ref2 <- logKs <- NULL
      allisna <- logical()
    }
    # get number of unavailable species
    nNA <- sum(!inbasis) + length(speciesNA) + sum(allisna)
    # save the information
    OUT[[type]] <- list(logKs = logKs, ref1 = ref1, ref2 = ref2, speciesOBIGT = speciesOBIGT, formula = formula, nNA = nNA)
  }
  OUT
}

