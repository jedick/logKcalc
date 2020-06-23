# writedat.R
# write a GWB data file with updated logK values
# 20200524

# a function to format values on separate lines 20200610
formatline <- function(values, iline, na.500 = FALSE, ndec = 4) {
  if(iline==1) values <- values[1:4]
  if(iline==2) values <- values[5:8]
  # use 500 for NA 20200526
  if(na.500) values[is.na(values)] <- 500
  if(ndec == 4) line <- paste0("   ", paste(sprintf("%12.4f", values), collapse = ""))
  if(ndec == 6) line <- paste0("   ", paste(sprintf("%12.6f", values), collapse = ""))
  line
}

writedat <- function(outfile, LINES, HEAD, LOGK, ADDS, infile, DH.method, a0_ion = NULL) {
  # put together the output
  # start with empty lines
  out <- rep(NA, length(LINES))
  # initialize a counter for the line of the output file
  j <- 1
  # flag to indicate whether to include the current basis species
  okbasis <- FALSE
  # line number with the name of the current species
  icurrent <- NA
  # content for reference line
  refline <- NA
  # calculate Debye-Hückel coefficients 20200622
  DH <- dh(LOGK$T, LOGK$P, DH.method)
  if(!is.null(a0_ion)) a0 <- sprintf("%5.2f", a0_ion)
  # calculate coefficients for CO2 and H2O activity 20200623
  CO2 <- co2(LOGK$T)
  H2O <- h2o(LOGK$T)
  # loop over lines of the input file
  for(i in 1:length(LINES)) {

    # we won't add anything if this is NA
    outline <- NA

    # check if we're in the top header block
    if(i < HEAD$ibasis) {
      # just copy the line from the input file
      outline <- LINES[i]
      # exclude comment lines 20200621
      if(grepl("^\\*", LINES[i]) & i < HEAD$iT) {
        outline <- NA
      }
      # output T and P 20200610
      if(i == HEAD$iT + 1) outline <- formatline(LOGK$T, 1)
      if(i == HEAD$iT + 2) outline <- formatline(LOGK$T, 2)
      if(i == HEAD$iP + 1) outline <- formatline(LOGK$P, 1)
      if(i == HEAD$iP + 2) outline <- formatline(LOGK$P, 2)
      # output Debye-Hückel parameters 20200622
      if(i == HEAD$iadh + 1) outline <- formatline(DH$adh, 1)
      if(i == HEAD$iadh + 2) outline <- formatline(DH$adh, 2)
      if(i == HEAD$ibdh + 1) outline <- formatline(DH$bdh * 1e-8, 1)
      if(i == HEAD$ibdh + 2) outline <- formatline(DH$bdh * 1e-8, 2)
      if(i == HEAD$ibdot + 1) outline <- formatline(DH$bdot, 1)
      if(i == HEAD$ibdot + 2) outline <- formatline(DH$bdot, 2)
      # output CO2 parameters 20200623
      if(i == HEAD$ico2_1 + 1) outline <- formatline(CO2$co2_1, 1, ndec = 6)
      if(i == HEAD$ico2_1 + 2) outline <- formatline(CO2$co2_1, 2, ndec = 6)
      if(i == HEAD$ico2_2 + 1) outline <- formatline(CO2$co2_2, 1, ndec = 6)
      if(i == HEAD$ico2_2 + 2) outline <- formatline(CO2$co2_2, 2, ndec = 6)
      if(i == HEAD$ico2_3 + 1) outline <- formatline(CO2$co2_3, 1, ndec = 6)
      if(i == HEAD$ico2_3 + 2) outline <- formatline(CO2$co2_3, 2, ndec = 6)
      if(i == HEAD$ico2_4 + 1) outline <- formatline(CO2$co2_4, 1, ndec = 6)
      if(i == HEAD$ico2_4 + 2) outline <- formatline(CO2$co2_4, 2, ndec = 6)
      # output H2O parameters 20200623
      if(i == HEAD$ih2o_1 + 1) outline <- formatline(H2O$h2o_1, 1, ndec = 6)
      if(i == HEAD$ih2o_1 + 2) outline <- formatline(H2O$h2o_1, 2, ndec = 6)
      if(i == HEAD$ih2o_2 + 1) outline <- formatline(H2O$h2o_2, 1, ndec = 6)
      if(i == HEAD$ih2o_2 + 2) outline <- formatline(H2O$h2o_2, 2, ndec = 6)
      if(i == HEAD$ih2o_3 + 1) outline <- formatline(H2O$h2o_3, 1, ndec = 6)
      if(i == HEAD$ih2o_3 + 2) outline <- formatline(H2O$h2o_3, 2, ndec = 6)
      if(i == HEAD$ih2o_4 + 1) outline <- formatline(H2O$h2o_4, 1, ndec = 6)
      if(i == HEAD$ih2o_4 + 2) outline <- formatline(H2O$h2o_4, 2, ndec = 6)
    }

    # check if we're in the basis species header
    if(i == HEAD$ibasis) {
      # update the number of basis species
      outline <- gsub(length(HEAD$ispecies$basis), length(LOGK$basis$map$GWB), LINES[HEAD$ibasis])
      message(paste("Writing output to", outfile, "..."))
      print(paste(length(LOGK$basis$map$GWB), "basis species"))
    }

    # check if we're in the basis species block
    if(i > HEAD$ibasis & i < HEAD$iredox) {
      # check if we're reading the name of a basis species
      if(i %in% HEAD$ispecies$basis) {
        # update the output flag
        if(LINES[i] %in% LOGK$basis$map$GWB) okbasis <- TRUE else okbasis <- FALSE
        icurrent <- i
        # include -end- marker
        if(i %in% HEAD$iend) okbasis <- TRUE
        # get the number of elements in the species 20200617
        for(nf in 0:1) {
          nspecies <- suppressWarnings(as.numeric(gsub(" ", "", gsub("elements in species", "", LINES[i + 2 + nf]))))
          if(!is.na(nspecies)) break
        }
        if(is.na(nspecies)) stop("can't find number of elements in species for ", LINES[i])
        # calculate the line number *after* all the element lines
        iafter <- i + 2 + ceiling(nspecies/3) + nf
      }
      # copy the line from the input file
      if(okbasis) {
        # don't copy comment lines from input file 20200617
        if(!grepl("^\\*", LINES[i])) {
          outline <- LINES[i]
          # adjust a0 parameter for ions with bgamma DH method 20200623
          if(!is.null(a0_ion)) {
            if(grepl("ion size", outline, fixed = TRUE)) {
              if(!grepl("charge=\\s*0", outline)) {
                start <- paste0(strsplit(outline, "ion size=", fixed = TRUE)[[1]][1], "ion size=")
                end <- paste0(" A", strsplit(outline, "A", fixed = TRUE)[[1]][2])
                outline <- paste0(start, a0, end)
              }
            }
          }
        }
        # add reference line 20200618
        if(i == iafter) {
          ispecies <- match(line2name(LINES[icurrent]), LOGK$basis$map$GWB)
          refline <- "* no references available"
          r1 <- LOGK$basis$ref1[ispecies]
          if(!is.na(r1)) {
            refline <- paste("* reference:", r1)
            r2 <- LOGK$basis$ref2[ispecies]
            if(!is.na(r2)) refline <- paste("* references: ", r1, ", ", r2, sep = "")
          }
        }
      }
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
      logKs <- LOGK[[type]]$logKs
      # also get the references 20200617
      ref1 <- LOGK[[type]]$ref1
      ref2 <- LOGK[[type]]$ref2
      # update the number of this type of species
      inames <- HEAD$ispecies[[type]]
      oldn <- length(inames)
      newn <- length(logKs)
      if(i %in% c(HEAD$iaqueous, HEAD$imineral, HEAD$igas)) newn <- newn + ADDS[[type]]$n
      outline <- gsub(oldn, newn, LINES[i])
      print(paste(newn, type, "species"))
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
          # account for extra intervening lines for minerals and some gases 20200526
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
          icurrent <- i
        }
      }

      # test to detect end of species type block
      if(i %in% HEAD$iend) {
        # include -end- marker and reset icurrent to NA
        outline <- LINES[i]
        icurrent <- NA
        # get lines for added species 20200621
        addlines <- ADDS[[type]]$lines
      } else if(okspecies) {
        # copy non-logK lines from the input file
        if(!i %in% c(ilogK1, ilogK2)) {
          # are we inside a species block?
          if(!is.na(icurrent)) {
            # don't copy comment lines from input file 20200617
            if(!grepl("^\\*", LINES[i])) {
              outline <- LINES[i]
              # adjust a0 parameter for bgamma DH method 20200623
              if(!is.null(a0_ion)) {
                if(grepl("ion size", outline, fixed = TRUE)) {
                  if(!grepl("charge=\\s*0", outline)) {
                    start <- paste0(strsplit(outline, "ion size=", fixed = TRUE)[[1]][1], "ion size=")
                    end <- paste0(" A", strsplit(outline, "A", fixed = TRUE)[[1]][2])
                    outline <- paste0(start, a0, end)
                  }
                }
              }
            }
          }
        } else {
          # write the logK or other data for a species
          # don't attempt to insert logK for oxide species 20200528
          if(i < HEAD$ioxide) {
            # get calculated logK values
            ispecies <- match(line2name(LINES[icurrent]), names(logKs))
            values <- logKs[[ispecies]]
            # insert calculated logK values
            if(i == ilogK1) iline <- 1
            if(i == ilogK2) iline <- 2
            outline <- formatline(values, iline, na.500 = TRUE)
            # add reference line 20200617
            if(i == ilogK2) {
              refline <- "* no references available"
              r1 <- ref1[ispecies]
              if(!is.na(r1)) {
                refline <- paste("* reference:", r1)
                r2 <- ref2[ispecies]
                if(!is.na(r2)) refline <- paste("* references: ", r1, ", ", r2, sep = "")
              }
            }
          } else outline <- ""
        }
      }

    }

    # keep a blank line after the headers
    if(i %in% (c(HEAD$ibasis, HEAD$iredox, HEAD$iaqueous, HEAD$ielectron, HEAD$imineral, HEAD$igas, HEAD$ioxide) + 1)) outline <- ""
    # keep a blank line after the -end- markers
    if(i %in% (HEAD$iend + 1)) outline <- ""
    # insert the lines for added species before the -end- markers 20200621
    if(i %in% HEAD$iend & i > HEAD$iredox) {
      nlines <- length(addlines)
      if(nlines > 0) out[j : (j + nlines - 1)] <- addlines
      j <- j + nlines
    }
    # add the line to the output
    if(!is.na(outline)) {
      out[j] <- outline
      j <- j + 1
    }
    # add the reference line to the output 20200617
    if(!is.na(refline)) {
      out[j] <- refline
      j <- j + 1
      refline <- NA
    }
    # insert comment block at top of file 20200621
    if(i == HEAD$iT - 2) {
      lver <- utils::packageDescription("logKcalc")$Version
      cver <- utils::packageDescription("CHNOSZ")$Version
      nadd <- sum(sapply(ADDS, "[[", "n"))
      nNA <- sum(unlist(sapply(LOGK, "[", "nNA")), na.rm = TRUE)
      clines <- c(
        paste0("* Thermodynamic database: OBIGT in CHNOSZ"),
        paste0("* File generated at ", date()),
        paste0("* by logKcalc ", lver, " with CHNOSZ ", cver, "."),
        paste0("* (https://github.com/jedick/logKcalc)"),
        paste0("* System based on ", basename(infile)),
        paste0("* with ", nadd, " added and ", nNA, " unavailable species.")
      )
      nlines <- length(clines)
      out[j : (j + nlines - 1)] <- clines
      j <- j + nlines
    }

  }

  # cleanup
  out <- stats::na.omit(out)
  # add reference block 20200618
  if(utils::packageVersion("CHNOSZ") > "1.3.6") {
    # get all reference keys used in output GWB file
    basisrefs <- c(LOGK$basis$ref1, LOGK$basis$ref2)
    speciestypes <- c("redox", "aqueous", "electron", "mineral", "gas")
    speciesrefs <- lapply(speciestypes, function(x) c(LOGK[[x]]$ref1, LOGK[[x]]$ref2))
    addrefs <- sapply(ADDS, "[[", "refs")
    allrefs <- c(basisrefs, speciesrefs, addrefs)
    keys <- sort(unique(stats::na.omit(unlist(allrefs))))
    # read the bibtex files from CHNOSZ and logKcalc
    bibfile1 <- system.file("doc/obigt.bib", package = "CHNOSZ")
    bibentry1 <- bibtex::read.bib(bibfile1)
    bibfile2 <- system.file("extdata/logKcalc.bib", package = "logKcalc")
    bibentry2 <- bibtex::read.bib(bibfile2)
    # use the current year for the logK_fit entry 20200623
    bibentry2["logK_fit"]$year <- substr(date(), 21, 24)
    bibentry <- c(bibentry1, bibentry2)
    # check for missing entries
    inbib <- keys %in% names(bibentry)
    if(any(!inbib)) warning(paste("reference(s) not found in bibtex file:", paste(keys[!inbib], collapse = ", ")))
    keys <- keys[inbib]
    # format the printed references
    op <- options(width = 90)
    on.exit(options(op))
    reftext <- utils::capture.output(print(bibentry[keys]))
    # insert the reference keys before the references
    keys <- paste0("[", keys, "]")
    refout <- character()
    j <- 1
    for(i in 1:length(reftext)) {
      if(i==1) { refout <- c(refout, keys[j]); j <- j + 1 }
      refout <- c(refout, reftext[i])
      if(reftext[i]=="") { refout <- c(refout, keys[j]); j <- j + 1 }
    }
    out <- c(out, refout)
  }

  # write the output to file
  # use CRLF lineendings (better for Windows users) 20200623
  writeLines(out, outfile, sep = "\r\n")
}
