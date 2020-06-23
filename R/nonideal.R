# logKcalc/nonideal.R
# Calculate parameters in equations for activity and osmotic coefficients
# First version 20200622

# Parameters for the extended Debye-Hückel equation
# as a function of T (degrees C) and P (bars)
dh <- function(T, P, method = "Bdot") {
  TK <- CHNOSZ::convert(T, "K")
  adh <- CHNOSZ::water("A_DH", T = TK, P = P)[, 1]
  bdh <- CHNOSZ::water("B_DH", T = TK, P = P)[, 1]
  if(identical(method, "bdot")) {
    # taken from CHNOSZ:::Bdot()
    bdot <- stats::splinefun(c(25, 50, 100, 150, 200, 250, 300), c(0.0418, 0.0439, 0.0468, 0.0479, 0.0456, 0.0348, 0))(T)
    bdot[T > 300] <- 0
  } else if(identical(method, "bgamma")) {
    bdot <- CHNOSZ::bgamma(T, P)
  } else stop("DH.method should be 'bdot' or 'bgamma'")
  list(adh = adh, bdh = bdh, bdot = bdot)
}

# Parameters for the activity coefficient of CO2 in NaCl solutions 20200623
# - equation from K2GWB (Table 1 of Cleverley and Bastrakov, 2005)
co2 <- function(T) {
  aterms <- c(-1.2089359E-15, 3.8431738E-12, -4.8786404E-09, 3.2625527E-06, -3.6089246E-04, 1.0678128E-01)
  bterms <- c(-3.1815653E-29, 1.0770230E-25, -1.4354227E-22, 9.8571174E-20, -1.0735937E-04, 3.8913425E-04)
  cterms <- c(3.3202383E-30, -1.2334246E-26, 1.7781153E-23, -1.3083434E-20, 9.2100868E-06, -3.3382837E-05)
  a <- b <- c <- d <- rep(0, length(T))
  exp <- 5:0
  for(i in 1:6) {
    a <- a + aterms[i] * T^exp[i]
    b <- b + bterms[i] * T^exp[i]
    c <- c + cterms[i] * T^exp[i]
  }
  list(co2_1 = a, co2_2 = b, co2_3 = c, co2_4 = d)
}

# Parameters for the osmotic coefficient of H2O 20200623
# - equation from K2GWB (Table 2 of Cleverley and Bastrakov, 2005)
h2o <- function(T) {
  aterms <- c(1.348253E-03, 1.420267E+00)
  bterms <- c(1.882726E-04, 1.765034E-02)
  cterms <- c(-3.889674E-05, 1.034523E-02)
  dterms <- c(-2.360000E-06, -4.772000E-04)
  a <- b <- c <- d <- rep(0, length(T))
  exp <- 1:0
  for(i in 1:2) {
    a <- a + aterms[i] * T^exp[i]
    b <- b + bterms[i] * T^exp[i]
    c <- c + cterms[i] * T^exp[i]
    d <- d + dterms[i] * T^exp[i]
  }
  list(h2o_1 = a, h2o_2 = b, h2o_3 = c, h2o_4 = d)
}
