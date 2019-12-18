
## from the first to the sixth orders
fun.taylor <- function(aa, x) {

  zero <- (aa < 1E-10); no_zero <- sum(zero)
  a <- aa[!zero]
  I1  <- (a + x) * (x + 1)
  I2  <- 2 * x + a + 1
  I3  <- log(a + x) - log(x + 1)

  I   <- I1 * I3
  II  <- I2 * I3 + 1 - a
  III <- 2 * I3 + I2 * (1 - a) / I1
  IV  <- - (1 - a)^3 / I1^2
  V   <- 2 * (1 - a)^3 * I2 / I1^3
  VI  <- 2 * (1 - a)^3 * (2 * I1 - 3 * I2^2) / I1^4
  VII <- 12 * (1 - a)^3 * (I2 * (2 * I2 - 3 * I1)) / I1^5


  f0 <- log(-I3)
  c <- - (1 - a)
  f1 <- - c / I
  f2 <- c * II / I^2
  f3 <- c * (III * I - 2 * II^2) / I^3
  f4 <- c * (IV * I^2 - 6 * I * II * III + 6 * II^3) / I^4
  f5 <- c * (I^3 * V - 8 * I^2 * II * IV - 6 * I^2 * III^2 + 36 * I * II^2 * III - 24 * II^4) / I^5
  f6 <- c * (-10 * I^3 * II * V + I^4 * VI - 20 * I^3 * III * IV + 60 * I^2 * II^2 * IV + 90 * I^2 * II * III^2
             - 240 * I * II^3 * III + 120 * II^5) / I^6

  if (no_zero > 0) {
      val <- rbind(f0, f1, f2, f3, f4, f5, f6)
      cbind(val, matrix(0, nrow(val), no_zero))
  } else {
      rbind(f0, f1, f2, f3, f4, f5, f6)
  }

}


## Sixth order
fun6 <- function(x, aa) {
  zero <- (aa < 1E-10); no_zero <- sum(zero)
  a <- aa[!zero]
  I1  <- (a + x) * (x + 1)
  I2  <- 2 * x + a + 1
  I3  <- log(a + x) - log(x + 1)

  I   <- I1 * I3
  II  <- I2 * I3 + 1 - a
  III <- 2 * I3 + I2 * (1 - a) / I1
  IV  <- - (1 - a)^3 / I1^2
  V   <- 2 * (1 - a)^3 * I2 / I1^3
  VI  <- 2 * (1 - a)^3 * (2 * I1 - 3 * I2^2) / I1^4
  VII <- 12 * (1 - a)^3 * (I2 * (2 * I2 - 3 * I1)) / I1^5

  c <- - (1 - a)
  f6 <- c * (-10 * I^3 * II * V + I^4 * VI - 20 * I^3 * III * IV + 60 * I^2 * II^2 * IV + 90 * I^2 * II * III^2
             - 240 * I * II^3 * III + 120 * II^5) / I^6
  if (no_zero > 0) {
    rep(0, length(x))
  } else {
    f6
  }
}


## Numerical approximation for high orders
f0.deriv6789 <- function(smix.chosen) {

  make.bfun <- function(p, x) {
    n <- length(x)
    x1 <- x[1:(n-1)]
    x2 <- x[2:n]
    if(all(x1 < x2) | all(x1 > x2)) {method <- "hyman"}
    else {method <- "fmm"}
    splinefun(p, x, method = method)
  }

  p1 <- seq(-0.001, 0.001, 0.00001)
  sapply(smix.chosen, function(a) {
    tFp <- fun6(p1, aa=a); id.nonNAN <- !is.na(tFp)
    Fp  <- tFp[id.nonNAN]; p1 <- p1[id.nonNAN]

    sfun6 <- make.bfun(p1, Fp)
    sfun7 <- function(p) {sfun6(p, deriv=1)}
    sfun8 <- function(p) {sfun6(p, deriv=2)}
    sfun9 <- function(p) {sfun6(p, deriv=3)}
    c(sfun6(0), sfun7(0), sfun8(0), sfun9(0))
  })

}
#f0.deriv6789(smix.chosen=c(0,0))


