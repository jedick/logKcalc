# logKcalc/nonideal.R
# calculate coefficients for activity and fugacity models
# first version 20200622

# Coefficients of the extended Debye-Hückel equation
# as a function of T (degrees C) and P (bars)
DH <- function(T, P, method = "Bdot") {
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
