# logKcomp.R
# compare logK values from two GWB files 20200525
# use ggrepel to label outliers 20200526

## for testing
#source("mapnames.R")
#source("utils.R")
#source("readhead.R")

logKcomp <- function(file1, file2, type = c("redox", "aqueous"), iTP = 2,
  lab1 = NULL, lab2 = NULL, xlim = NULL, ylim = NULL, plot.it = TRUE) {

  # for missing file2, just get values from file1 20200702
  missfile2 <- FALSE
  if(missing(file2)) {
    message("file2 was missing; running file2 <- file1")
    file2 <- file1
    missfile2 <- TRUE
  }

  # get data titles from file names
  if(is.null(lab1)) lab1 <- basename(file1)
  if(is.null(lab2)) lab2 <- basename(file2)
  # read logK values
  LOGK1 <- readlogK(file1)
  LOGK2 <- readlogK(file2)

  # loop over types 20200702
  logKs <- lapply(type, function(type) {
    # get logKs for the indicated species type
    logK1 <- LOGK1[[type]]
    logK2 <- LOGK2[[type]]
    # check for no species of this type
    if(length(logK1)==0 | length(logK2)==0) {
      if(length(logK1)==0 & length(logK2)==0) message(paste("No logKs for", type, "species in file1 or file2"))
      else if(length(logK1)==0) message(paste("No logKs for", type, "species in file1"))
      else if(length(logK2)==0) message(paste("No logKs for", type, "species in file2"))
      return(data.frame(logK1 = numeric(), logK2 = numeric()))
    }
    # remove (aq) suffix 20200616
    names(logK1) <- gsub("\\(aq\\)$", "", names(logK1))
    names(logK2) <- gsub("\\(aq\\)$", "", names(logK2))
    # get name of species that are present in both files
    names <- intersect(names(logK1), names(logK2))
    # get the matching indexes for the names
    i1 <- match(names, names(logK1))
    i2 <- match(names, names(logK2))
    # get *mapped* names of species that are present in both files 20200627
    map1 <- mapnames(names(logK1))$OBIGT
    map2 <- mapnames(names(logK2))$OBIGT
    mappednames <- stats::na.omit(intersect(map1, map2))
    j1 <- match(mappednames, map1)
    j2 <- match(mappednames, map2)
    # remove matching mapped names that already have matching *un*mapped names
    itxt <- paste(i1,i2)
    jtxt <- paste(j1,j2)
    inotmatched <- !jtxt %in% itxt
    if(any(inotmatched)) {
      j1 <- j1[inotmatched]
      j2 <- j2[inotmatched]
      # add these to the matching index and print update
      i1 <- c(i1, j1)
      i2 <- c(i2, j2)
      for(k in 1:length(j1)) {
        message(paste(names(logK1)[j1][k], "in file1 matches", names(logK2)[j2][k], "in file2 via mapped name", map1[j1][k]))
      }
    }
    # print names of unmatched species 20200616
    not1 <- names(logK1)[-i1]
    if(length(not1) > 0) {
      message(paste("Didn't match", length(not1), type, "species in file 1:"))
      print(paste(not1, collapse = ", "))
    }
    not2 <- names(logK2)[-i2]
    if(length(not2) > 0) {
      message(paste("Didn't match", length(not2), type, "species in file 2:"))
      print(paste(not2, collapse = ", "))
    }
    # extract the values for plotting
    logK1 <- logK1[i1]
    logK2 <- logK2[i2]

    # check for reaction balance 20200616
    LINES1 <- cleanUTF8(readLines(file1), file1)
    LINES2 <- cleanUTF8(readLines(file2), file2)
    HEAD1 <- readhead(LINES1, quiet = TRUE)
    HEAD2 <- readhead(LINES2, quiet = TRUE)
    ispecies1 <- HEAD1$ispecies[[type]]
    ispecies2 <- HEAD2$ispecies[[type]]
    iunequal <- numeric()
    for(i in 1:length(names)) {
      j1 <- i1[i]
      j2 <- i2[i]
      rxn1 <- readrxn(LINES1, HEAD1$ihead, ispecies1[j1])
      rxn2 <- readrxn(LINES2, HEAD2$ihead, ispecies2[j2])
      # sort the species in the reactions
      rxn1$coeff <- rxn1$coeff[order(rxn1$species)]
      rxn1$species <- rxn1$species[order(rxn1$species)]
      rxn2$coeff <- rxn2$coeff[order(rxn2$species)]
      rxn2$species <- rxn2$species[order(rxn2$species)]
      if(!identical(rxn1, rxn2)) iunequal <- c(iunequal, i)
    }
    if(length(iunequal) > 0) {
      message(paste("Excluding", length(iunequal), type, "species with unequal reactions in the files:"))
      print(paste(names[iunequal], collapse = ", "))
      logK1 <- logK1[-iunequal]
      logK2 <- logK2[-iunequal]
    }
    message(paste("Found", length(logK1), type, "species with same names and reactions in both files"))

    # get logK values for all species at this T, P
    logK_1 <- sapply(logK1, "[", iTP)
    logK_2 <- sapply(logK2, "[", iTP)
    # change "500" values to NA
    logK_1[logK_1==500] <- NA
    logK_2[logK_2==500] <- NA
    data.frame(logK1 = logK_1, logK2 = logK_2)
  })
  # put together data for each type of species 20200702
  logKs <- do.call(rbind, logKs)

  # if file2 was missing, return the values from file1 20200702
  if(missfile2) {
    logK1 <- logKs[, 1]
    names(logK1) <- row.names(logKs)
    return(logK1)
  }
  # if we're not making a plot, return the values 20200616
  if(!plot.it | nrow(logKs)==0) return(logKs)

  # get x, y values and identify outliers
  xy <- data.frame(
    x = logKs$logK1,
    y = logKs$logK2 - logKs$logK1,
    row.names = row.names(logKs)
  )
  ylab = paste(lab2, "-", lab1)
  iout <- abs(xy$y) > 1
  # get labels for outliers and plot title
  labels <- rownames(xy)
  labels[!iout] <- ""
  title <- bquote(log~italic(K)~of~.(type)~species~at~.(LOGK1$T[iTP])~"\u00BAC"~and~.(LOGK1$P[iTP])~bar)
  subtitle <- paste0("NA values: ", sum(is.na(logKs$logK1)), " in ", lab1, ", ", sum(is.na(logKs$logK2)), " in ", lab2)
  x <- y <- NULL
  ggplot2::ggplot(xy, ggplot2::aes(x, y, label = labels)) +
    ggrepel::geom_text_repel(na.rm = TRUE, max.overlaps = 20) +
    ggplot2::geom_point(color = 'red', na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::xlab(lab1) +
    ggplot2::ylab(ylab) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "gray30") +
    ggplot2::ggtitle(title, subtitle) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
}
