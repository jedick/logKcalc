# readlogK.R
# read logK values from data file in GWB format
# 20200523 jmd initial version
# 20200525 revision to move not-reading logic to calclogK()

# a function to exclude type= ... at end of line with names 20200526
# also remove formula= (occurs in thermo.com.V8.R6+.tdat) 20200623
line2name <- function(LINE) {
  LINE <- gsub("\ *type=.*", "", LINE)
  LINE <- gsub("\ *formula=.*", "", LINE)
  LINE
}

readlogK <- function(file, quiet = FALSE) {
  # read the GWB data file
  message(paste("Reading file", file))
  LINES <- readLines(file)
  # get the header data
  HEAD <- readhead(LINES, quiet = quiet)
  # make a list to hold information on the different types of species
  LOGK <- list(T = HEAD$T, P = HEAD$P, redox = NA, aqueous = NA)
  # get logK for each types of species
  for(type in c("redox", "aqueous", "electron", "mineral", "gas")) {
    # get the line numbers with the species names
    inames <- HEAD$ispecies[[type]]
    # read the species names
    speciesGWB <- line2name(LINES[inames])
    logKs <- numeric()
    # loop over species
    logKs <- lapply(inames, function(i) {
      # get the number of species in the reaction
      # account for an extra line for minerals and some gases 20200526
      for(nf in 0:2) {
        nspecies <- suppressWarnings(as.numeric(gsub(" ", "", gsub("species in reaction", "", LINES[i + 2 + nf]))))
        if(!is.na(nspecies)) break
      }
      if(is.na(nspecies)) stop("can't find number of species in reaction for ", LINES[i])
      # calculate the number of reaction lines
      nrxnlines <- ceiling(nspecies/3)
      # identify the lines with the logK values
      ilogK1 <- i + nrxnlines + 3 + nf
      ilogK2 <- i + nrxnlines + 4 + nf
      # parse the logK values
      logK1 <- LINES[ilogK1]
      logK1 <- strsplit(logK1, " ")[[1]]
      logK1 <- as.numeric(logK1[nzchar(logK1)])
      logK2 <- LINES[ilogK2]
      logK2 <- strsplit(logK2, " ")[[1]]
      logK2 <- as.numeric(logK2[nzchar(logK2)])
      # return the values
      c(logK1, logK2)
    })
    names(logKs) <- speciesGWB
    # save the information
    LOGK[[type]] <- logKs
  }
  LOGK
}
