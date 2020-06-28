# logKcalc/utils.R
# utility functions 20200628

# Function to create reaction lines 20200620
# used in addspecies, readlogK, writedat
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

# a function to read the *reaction* coefficients and species
# used in calclogK, logKcomp, modOBIGT
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

# a function to exclude type= ... at end of line with names 20200526
# also remove formula= (occurs in thermo.com.V8.R6+.tdat) 20200623
# used in calclogK, modOBIGT, readlogK, writedat
line2name <- function(LINE) {
  LINE <- gsub("\ *type=.*", "", LINE)
  LINE <- gsub("\ *formula=.*", "", LINE)
  LINE
}

# a function to format values on separate lines 20200610
# used in addspecies, writedat
formatline <- function(values, iline, na.500 = FALSE, ndec = 4) {
  if(iline==1) values <- values[1:4]
  if(iline==2) values <- values[5:8]
  # use 500 for NA 20200526
  if(na.500) values[is.na(values)] <- 500
  if(ndec == 4) line <- paste0("   ", paste(sprintf("%12.4f", values), collapse = ""))
  if(ndec == 6) line <- paste0("   ", paste(sprintf("%12.6f", values), collapse = ""))
  line
}

