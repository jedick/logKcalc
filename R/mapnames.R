# mapnames.R
# 20200505 map GWB to OBIGT species names (only those listed in map_names.csv)
# 20200525 also map names that match names or formulas in OBIGT
#   - possibly useful trick: set NA for a CHNOSZ name in map_names.csv to remove that species

mapnames <- function(GWBnames, type = NULL, na.omit = FALSE, return.processed.name = FALSE) {

  # remove (aq) and ,aq suffixes 20200610
  Gnames <- gsub("\\(aq\\)$", "", GWBnames)
  Gnames <- gsub(",aq$", "", Gnames)
  # replace underscore with space 20200611
  Gnames <- gsub("_", " ", Gnames)
  # remove leading "1-" (e.g. 1-Pentanol)
  Gnames <- gsub("^1-", "", Gnames)
  # replace repeated "-" and "+" with number 20200611
  Gnames <- gsub("\\-\\-\\-\\-\\-\\-\\-\\-\\-\\-$", "-10", Gnames)
  Gnames <- gsub("\\-\\-\\-\\-\\-\\-\\-\\-\\-$", "-9", Gnames)
  Gnames <- gsub("\\-\\-\\-\\-\\-\\-\\-\\-$", "-8", Gnames)
  Gnames <- gsub("\\-\\-\\-\\-\\-\\-\\-$", "-7", Gnames)
  Gnames <- gsub("\\-\\-\\-\\-\\-\\-$", "-6", Gnames)
  Gnames <- gsub("\\-\\-\\-\\-\\-$", "-5", Gnames)
  Gnames <- gsub("\\-\\-\\-\\-$", "-4", Gnames)
  Gnames <- gsub("\\-\\-\\-$", "-3", Gnames)
  Gnames <- gsub("\\-\\-$", "-2", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+\\+\\+\\+\\+\\+\\+$", "+10", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+\\+\\+\\+\\+\\+$", "+9", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+\\+\\+\\+\\+$", "+8", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+\\+\\+\\+$", "+7", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+\\+\\+$", "+6", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+\\+$", "+5", Gnames)
  Gnames <- gsub("\\+\\+\\+\\+$", "+4", Gnames)
  Gnames <- gsub("\\+\\+\\+$", "+3", Gnames)
  Gnames <- gsub("\\+\\+$", "+2", Gnames)
  if(return.processed.name) return(Gnames)

  # first look in OBIGT
  OBIGT <- get("thermo", CHNOSZ)$obigt
  if(!is.null(type)) {
    if(type %in% c("basis", "redox", "aqueous", "electron")) OBIGT <- OBIGT[OBIGT$state=="aq", ]
    if(type == "mineral") OBIGT <- OBIGT[OBIGT$state=="cr", ]
    if(type == "gas") OBIGT <- OBIGT[OBIGT$state=="gas", ]
  }
  # H2O is a liquid, so doesn't fall in the other categories
  allOBIGTnames <- unique(c("H2O", OBIGT$name, OBIGT$formula))
  iOBIGT <- match(Gnames, allOBIGTnames)
  OBIGTnames <- allOBIGTnames[iOBIGT]
  # also try lower-case names
  iOBIGT <- match(tolower(Gnames), allOBIGTnames)
  ina <- is.na(OBIGTnames)
  OBIGTnames[ina] <- allOBIGTnames[iOBIGT[ina]]

  # then look in map_names.csv
  mapfile <- system.file("extdata/map_names.csv", package = "logKcalc")
  map <- utils::read.csv(mapfile, as.is = TRUE)
  # earlier versions of CHNOSZ use lowercase "a" for Angstrom abbreviation in mineral names 20200609
  if(!utils::packageVersion("CHNOSZ") > "1.3.6") {
    map$OBIGT <- gsub("7A$", "7a", map$OBIGT)
    map$OBIGT <- gsub("14A$", "14a", map$OBIGT)
  }
  imap <- match(GWBnames, map$GWB)
  # get the names of the mapped species
  iimap <- !is.na(imap)
  mappednames <- map$OBIGT[imap[iimap]]
  # check that they're actually present in OBIGT 20200614
  infonames <- suppressMessages(CHNOSZ::info(mappednames))
  iimap[iimap][is.na(infonames)] <- FALSE
  OBIGTnames[iimap] <- map$OBIGT[imap[iimap]]

  if(na.omit) {
    # remove species that are NA in CHNOSZ
    ina <- is.na(OBIGTnames)
    if(any(ina)) {
      message("Can't map these ", type, " species to the OBIGT database:")
      print(paste(GWBnames[ina], collapse = ", "))
      GWBnames <- GWBnames[!ina]
      OBIGTnames <- OBIGTnames[!ina]
    }
  } else {
    # identify species that are NA in CHNOSZ
    GWBnames[is.na(OBIGTnames)] <- NA
  }

  # return the mapping
  list(GWB = GWBnames, OBIGT = OBIGTnames)
}
