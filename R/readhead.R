# readhead.R
# read header information from a GWB data file 20200523

# a function to read the *numeric* values in a data block
# 'i' is the line number for the header line in the file
readblock <- function(LINES, ihead, i) {
  # find the next header line after this one
  j <- min(ihead[ihead > i])
  # get all lines between lines i and j
  k <- (i + 1):(j - 1)
  thisdat <- LINES[k]
  # use the line numbers as names for the list
  names(thisdat) <- k
  # separate the content and remove blanks
  thisdat <- strsplit(thisdat, " ")
  thisdat <- lapply(thisdat, function(x) x[nzchar(x)])
  # also remove empty lines
  thisdat <- thisdat[! sapply(thisdat, length) == 0]
  # convert values to numeric
  thisdat <- suppressWarnings(lapply(thisdat, as.numeric))
  # remove any lines that have non-numeric values
  hasNA <- sapply(thisdat, function(x) any(is.na(x)))
  thisdat <- thisdat[!hasNA]
  thisdat
}

readhead <- function(LINES, quiet = FALSE) {
  # find header lines before data blocks
  iT <- grep("^\\*\ temperatures\\s*$", LINES)
  iP <- grep("^\\*\ pressures\\s*$", LINES)
  iadh <- grep("^\\*\ debye.*adh", LINES)
  ibdh <- grep("^\\*\ debye.*bdh", LINES)
  ibdot <- grep("^\\*\ bdot\\s*$", LINES)
  ico2_1 <- grep("* c co2 1", LINES, fixed = TRUE)
  ico2_2 <- grep("* c co2 2", LINES, fixed = TRUE)
  ico2_3 <- grep("* c co2 3", LINES, fixed = TRUE)
  ico2_4 <- grep("* c co2 4", LINES, fixed = TRUE)
  ih2o_1 <- grep("* c h2o 1", LINES, fixed = TRUE)
  ih2o_2 <- grep("* c h2o 2", LINES, fixed = TRUE)
  ih2o_3 <- grep("* c h2o 3", LINES, fixed = TRUE)
  ih2o_4 <- grep("* c h2o 4", LINES, fixed = TRUE)
  # find header lines before sections
  jbasis <- grep("basis species$", LINES)
  jredox <- grep("redox couples$", LINES)
  jaqueous <- grep("aqueous species$", LINES)
  jelectron <- grep("free electron$", LINES)
  jmineral <- grep("minerals$", LINES)
  jgas <- grep("gases$", LINES)
  joxide <- grep("oxides$", LINES)
  # find -end- markers
  iend <- grep("^-end-$", LINES)
  # --- find header lines ---
  # non-empty lines that start with a non-blank character
  ihead <- grep("^\\ ", LINES, invert = TRUE)
  ihead <- ihead[! LINES[ihead] == ""]
  # also exclude lines beginning with "-" but not "-end-"
  # (some reactions output by DBCreate have negative sign in first column)
  ihead <- ihead[! substr(LINES[ihead], 1, 1) == "-" | LINES[ihead] == "-end-" ]
  # --- find lines that have names of species ---
  # exclude -end- markers and lines beginning with "*" (comment lines)
  iname <- ihead[! ihead %in% iend & ! substr(LINES[ihead], 1, 1) == "*"]
  # find lines that have the names of species
  ibasis <- iname[iname > jbasis & iname < jredox]
  iredox <- iname[iname > jredox & iname < jaqueous]
  iaqueous <- iname[iname > jaqueous & iname < jelectron]
  # in case there is no electron block (as in K2GWB output - dataset format: oct94) 20200610
  if(length(iaqueous)==0) iaqueous <- iname[iname > jaqueous & iname < jmineral]
  ielectron <- iname[iname > jelectron & iname < jmineral]
  imineral <- iname[iname > jmineral & iname < jgas]
  igas <- iname[iname > jgas & iname < joxide]
  ioxide <- iname[iname > joxide & iname < utils::tail(iend, 1)]

  # get the temperature and pressure
  T <- as.numeric(unlist(readblock(LINES, ihead, iT)))
  if(!quiet) {
    message("The temperatures in the input file are:")
    print(paste(T, collapse = ", "))
  }
  P <- as.numeric(unlist(readblock(LINES, ihead, iP)))
  if(!quiet) {
    message("The pressures in the input file are:")
    print(paste(P, collapse = ", "))
  }

  # put together output
  out <- list(
    iT = iT,
    iP = iP,
    iadh = iadh,
    ibdh = ibdh,
    ibdot = ibdot,
    ico2_1 = ico2_1,
    ico2_2 = ico2_2,
    ico2_3 = ico2_3,
    ico2_4 = ico2_4,
    ih2o_1 = ih2o_1,
    ih2o_2 = ih2o_2,
    ih2o_3 = ih2o_3,
    ih2o_4 = ih2o_4,
    ibasis = jbasis,
    iredox = jredox,
    iaqueous = jaqueous,
    ielectron = jelectron,
    imineral = jmineral,
    igas = jgas,
    ioxide = joxide,
    iend = iend,
    ihead = ihead,
    iname = iname,
    ispecies = list(
      basis = ibasis,
      redox = iredox,
      aqueous = iaqueous,
      electron = ielectron,
      mineral = imineral,
      gas = igas,
      oxide = ioxide
    ),
    T = T,
    P = P
  )
  out
}
