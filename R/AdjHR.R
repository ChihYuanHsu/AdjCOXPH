#' An adjustment for Cox hazard ratio and an estimation for difference in proportions.
#' @description A novel adjustment method for Cox proportional hazards model in data with long-term survival.
#' @import stats



#' @param HR_cox Input. Hazard ratio obtained from Cox PH model.
#' @param HR_cox_CI Input. Confidence interval of the hazard ratio obtained from Cox PH model.
#' @param s1mix.chosen Input. Survival probabilities of Arm 1 estimated by KM at the chosen time.
#' @param s0mix.chosen Input. Survival probabilities of Arm 0 estimated by KM at the chosen time.
#' @param pi1.est Input. Proportion of poor-responders (uncured proportion) in Arm 1.
#' @param pi0.est Input. Proportion of poor-responders (uncured proportion) in Arm 0.
#' @param m Input (Optional). The polynomial order used in Taylor approximation.
#' The default value is NULL that means \emph{m} is selected automatically.

#' @return \item{HR_cox}{Cox hazard ratio you inputed.}
#' @return \item{HR_cox_CI}{Confidence interval of the Cox hazard ratio you inputed.}
#' @return \item{HR_adj}{Hazard ratio after adjustment.}
#' @return \item{HR_adj_CI}{Confidence interval of the hazard ratio after adjustment.}
#' @return \item{DP_adj}{Difference in proportions of the true responders.}
#' @return \item{DP_adj_CI}{Confidence interval of the difference in proportions of the true responders.}


#' @examples library(AdjCOXPH)
#' @examples s1mix.chosen <- c(0.48,0.39,0.35,0.35)
#' @examples s0mix.chosen <- c(0.36,0.28,0.26,0.25)
#' @examples pi1.est <- 0.65; pi0.est <- 0.75
#' @examples HR_cox <- 0.71; HR_cox_CI <- c(0.51, 0.91)
#' @examples adjustment(HR_cox, HR_cox_CI, s1mix.chosen, s0mix.chosen, pi1.est, pi0.est)
#'

#' @export

adjustment <- function(HR_cox, HR_cox_CI, s1mix.chosen, s0mix.chosen, pi1.est, pi0.est, m = NULL) {

  rule <- function(cure_rate) {
    if (cure_rate <= .3) {
      k <- 4
    } else if (cure_rate > 0.7) {
      k <- 9
    } else {
      k <- ceiling(10 * cure_rate + 1)
    }
    return(k)
  }
  #rule(0.5)

  if (is.null(m)) {
    #K1 <- rule(1 - pi1.est); K0 <- rule(1 - pi0.est)
    K1 <- K0 <- max(rule(1 - pi1.est), rule(1 - pi0.est))
  } else {
    #K1 <- m[1]; K0 <- m[2]
    K1 <- K0 <- m
  }

  #### Taylor expansion
  if (K1 <= 6) {
    s1_06  <- fun.taylor(aa = s1mix.chosen, x=0)
    T1 <- s1_06
  } else {
    s1_06  <- fun.taylor(aa = s1mix.chosen, x=0)
    s1_789 <- f0.deriv6789(s1mix.chosen)[-c(1),]
    T1 <- rbind(s1_06, s1_789)
  }
  if (K0 <= 6) {
    s0_06  <- fun.taylor(aa = s0mix.chosen, x=0)
    T0 <- s0_06
  } else {
    s0_06  <- fun.taylor(aa = s0mix.chosen, x=0)
    s0_789 <- f0.deriv6789(s0mix.chosen)[-c(1),]
    T0 <- rbind(s0_06, s0_789)
  }
  c1 <- (pi1.est - 1)^(0:(nrow(T1)-1))/gamma(1:nrow(T1))
  c0 <- (pi0.est - 1)^(0:(nrow(T0)-1))/gamma(1:nrow(T0))
  CF <- exp(mean( colSums( (T1 * c1)[2:(K1+1),] ) - colSums( (T0 * c0)[2:(K0+1),] ) ))
  #CF <- 7
  #exp(mean( sum( (T1 * c1)[1,] ) - sum( (T0 * c0)[1,] ) )); HR_cox
  #mean( colSums( (T1 * c1)[2:5,] ) - colSums( (T0 * c0)[2:5,] ) )

  HR_cal <- CF * HR_cox; HR_cal_CI <- CF * HR_cox_CI

  adj.factor <- function(pi1.est, pi0.est) {
    c1 <- (pi1.est - 1)^(0:(nrow(T1)-1))/gamma(1:nrow(T1))
    c0 <- (pi0.est - 1)^(0:(nrow(T0)-1))/gamma(1:nrow(T0))
    F_adj <- exp(mean( colSums( (T1 * c1)[2:(K1+1),] ) - colSums( (T0 * c0)[2:(K0+1),] ) ))
    F_adj
  }

  #### CI of DP
  solve1 <- function(pi1.limit, HRcox.limit) {
    adj.factor(pi1.limit, pi0.est) * HRcox.limit - HR_cal
  }

  CI.pi1 <- c(tryCatch({
    aa <- uniroot(f = solve1, interval = c(0.001, .999), HRcox.limit = HR_cox_CI[1])$root
  }, error = function(e) {
    if (abs(solve1(0.001, HR_cox_CI[1])) > abs(solve1(0.999, HR_cox_CI[1]))) {
      aa <- 1
    } else {
      aa <- 0
    }
    return(aa)
  }), tryCatch({
    aa <- uniroot(f = solve1, interval = c(0.001, .999), HRcox.limit = HR_cox_CI[2])$root
  }, error = function(e) {
    if (abs(solve1(0.001, HR_cox_CI[2])) > abs(solve1(0.999, HR_cox_CI[2]))) {
      aa <- 1
    } else {
      aa <- 0
    }
    return(aa)
  }) )


  solve0 <- function(pi0.limit, HRcox.limit) {
    adj.factor(pi1.est, pi0.limit) * HRcox.limit - HR_cal
  }

  CI.pi0 <- c(tryCatch({
    aa <- uniroot(f = solve0, interval = c(0.001, .999), HRcox.limit = HR_cox_CI[1])$root
  }, error = function(e) {
    if (abs(solve0(0.001, HR_cox_CI[1])) > abs(solve0(0.999, HR_cox_CI[1]))) {
      aa <- 1
    } else {
      aa <- 0
    }
    return(aa)
  }), tryCatch({
    aa <- uniroot(f = solve0, interval = c(0.001, .999), HRcox.limit = HR_cox_CI[2])$root
  }, error = function(e) {
    if (abs(solve0(0.001, HR_cox_CI[2])) > abs(solve0(0.999, HR_cox_CI[2]))) {
      aa <- 1
    } else {
      aa <- 0
    }
    return(aa)
  }))


  pair4 <- c(max(CI.pi0) - CI.pi1, min(CI.pi0) - CI.pi1)
  DP_CI <- c(min(pair4), max(pair4))

  round(c(HR_cox = HR_cox, HR_cox_CI = HR_cox_CI,
          HR_adj = HR_cal, HR_adj_CI = HR_cal_CI,
          DP_adj = pi0.est - pi1.est, DP_adj_CI = DP_CI), 4)


}
