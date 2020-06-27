# logKcalc/vignettes/As_mosaic.R
# R code for As-S-O-H mosaic diagram
# 20200624 jmd first version

# Use default OBIGT database
obigt()

# Get As-bearing species in the As-S-O-H system
iaq <- retrieve("As", c("S", "O", "H"), "aq")
icr <- retrieve("As", c("S", "O", "H"), "cr")

# Remove species that don't have the parameters needed for high-T calculations
sout.aq <- subcrt(iaq, T = 100, property = "G")$out
na.aq <- is.na(unlist(sout.aq))
iaq <- iaq[!na.aq]
sout.cr <- subcrt(icr, T = 100, property = "G")$out
na.cr <- is.na(unlist(sout.cr))
icr <- icr[!na.cr]

# Set up basis species and system
basis(c("H3AsO3", "SO4-2", "H2O", "O2", "H+"))
basis("SO4-2", -4)
species(iaq, -3)
species(icr, 0)

# Set pH and loga(O2) range and grid resolution
pH <- c(0, 14, 500)
O2 <- c(-70, -40, 500)

# Calculate affinities for the mosaic diagram
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
m <- mosaic(bases, pH = pH, O2 = O2, T = 100)

# Make the diagram
par(cex = 2.5, mar = c(2.8, 3, 1.4, 1))
fill <- c(rep("aliceblue", length(iaq)), rep("antiquewhite1", length(icr)))
d <- diagram(m$A.species, fill = fill, lwd = 2)
diagram(m$A.bases, add = TRUE, col = "blue", lty = 3, lwd = 2, col.names = "blue", cex.names = 0.7, italic = TRUE)
water.lines(d)
legend("bottomleft", describe.property("T", 100), bty = "n")
title(main = "Default OBIGT database", font.main = 1)
