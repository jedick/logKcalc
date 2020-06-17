# writedat.R
# write a GWB data file with updated logK values
# 20200524

# a function to create one line of logK values
logKline <- function(LINES, iname, logKs, iline) {
  GWBnames <- names(logKs)
  ispecies <- match(line2name(LINES[iname]), GWBnames)
  if(iline==1) thislogK <- logKs[[ispecies]][1:4]
  if(iline==2) thislogK <- logKs[[ispecies]][5:8]
  # use 500 for NA 20200526
  thislogK[is.na(thislogK)] <- 500
  paste0("   ", paste(sprintf("%12.4f", thislogK), collapse = ""))
}

# a function to create one line of T or P values 20200610
TPline <- function(x, iline) {
  if(iline==1) x <- x[1:4]
  if(iline==2) x <- x[5:8]
  paste0("   ", paste(sprintf("%12.4f", x), collapse = ""))
}

writedat <- function(outfile, LINES, HEAD, LOGK) {
  # put together the output
  # start with empty lines
  out <- rep(NA, length(LINES))
  # initialize a counter for the line of the output file
  j <- 1
  # flag to indicate whether to include the current basis species
  okbasis <- FALSE
  icurrent <- NA
  # loop over lines of the input file
  for(i in 1:length(LINES)) {

    # we won't add anything if this is NA
    outline <- NA

    # check if we're in the top header block
    if(i < HEAD$ibasis) {
      # just copy the line from the input file
      outline <- LINES[i]
      # change T and P 20200610
      if(i == HEAD$iT + 1) outline <- TPline(LOGK$T, 1)
      if(i == HEAD$iT + 2) outline <- TPline(LOGK$T, 2)
      if(i == HEAD$iP + 1) outline <- TPline(LOGK$P, 1)
      if(i == HEAD$iP + 2) outline <- TPline(LOGK$P, 2)
    }

    # check if we're in the basis species header
    if(i == HEAD$ibasis) {
      # update the number of basis species
      outline <- gsub(length(HEAD$ispecies$basis), length(LOGK$basismap$GWB), LINES[HEAD$ibasis])
      message(paste("Writing output to", outfile, "..."))
      print(paste("outputting", length(LOGK$basismap$GWB), "basis species"))
    }

    # check if we're in the basis species block
    if(i > HEAD$ibasis & i < HEAD$iredox) {
      # check if we're reading the name of a basis species
      if(i %in% HEAD$ispecies$basis) {
        # update the flag
        if(LINES[i] %in% LOGK$basismap$GWB) okbasis <- TRUE else okbasis <- FALSE
        # include -end- marker
        if(i %in% HEAD$iend) okbasis <- TRUE
      }
      # copy the line from the input file
      if(okbasis) outline <- LINES[i]
    }

    # check if we're in a species type header
    if(i %in% c(HEAD$iredox, HEAD$iaqueous, HEAD$ielectron, HEAD$imineral, HEAD$igas, HEAD$ioxide)) {
      if(identical(i, HEAD$iredox)) type <- "redox"
      if(identical(i, HEAD$iaqueous)) type <- "aqueous"
      if(identical(i, HEAD$ielectron)) type <- "electron"
      if(identical(i, HEAD$imineral)) type <- "mineral"
      if(identical(i, HEAD$igas)) type <- "gas"
      if(identical(i, HEAD$ioxide)) type <- "oxide"
      # get the logK values for this type of species
      logKs <- LOGK[[type]]
      # update the number of this type of species
      inames <- HEAD$ispecies[[type]]
      outline <- gsub(length(inames), length(logKs), LINES[i])
      print(paste("outputting", length(logKs), type, "species"))
      # initialize the species inclusion flag
      okspecies <- FALSE
    }

    if(i > HEAD$iredox) {

      # check if we're reading the name of a species
      if(i %in% inames) {
        # test if the species is in our database
        if(!line2name(LINES[i]) %in% names(logKs)) okspecies <- FALSE else {
          okspecies <- TRUE
          # get the number of species in the reaction
          # account for an extra line for minerals and some gases 20200526
          for(nf in 0:1) {
            nspecies <- suppressWarnings(as.numeric(gsub(" ", "", gsub("species in reaction", "", LINES[i + 2 + nf]))))
            if(!is.na(nspecies)) break
          }
          # calculate the number of reaction lines
          nrxnlines <- ceiling(nspecies/3)
          # identify the lines with the logK values
          ilogK1 <- i + nrxnlines + 3 + nf
          ilogK2 <- i + nrxnlines + 4 + nf
          # keep the line number with the name of the current species
          icurrent <- i
        }
      }

      # test to detect end of species type block
      if(i %in% HEAD$iend) {
        # include -end- marker and reset icurrent to NA
        outline <- LINES[i]
        icurrent <- NA
      } else if(okspecies) {
        if(!i %in% c(ilogK1, ilogK2)) {
          # if we're inside a species block, copy the line from the input file
          if(!is.na(icurrent)) outline <- LINES[i]
        } else {
          # don't attempt to insert logK for oxide species 20200528
          if(i < HEAD$ioxide) {
            # insert calculated logK values
            if(i == ilogK1) iline <- 1
            if(i == ilogK2) iline <- 2
            outline <- logKline(LINES, icurrent, logKs, iline)
          } else outline <- ""
        }
      }

    }

    # keep a blank line after the headers
    if(i %in% (c(HEAD$ibasis, HEAD$iredox, HEAD$iaqueous, HEAD$ielectron, HEAD$imineral, HEAD$igas, HEAD$ioxide) + 1)) outline <- ""
    # keep a blank line after the -end- markers
    if(i %in% (HEAD$iend + 1)) outline <- ""
    # add the line to the output
    if(!is.na(outline)) {
      out[j] <- outline
      j <- j + 1
    }

  }

  # cleanup and write the output to file
  out <- stats::na.omit(out)
  writeLines(out, outfile)
}
