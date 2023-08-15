#' Cox PH – Taylor expansion adjustment for Long-term survival data
#' @description An adjustment for Cox hazard ratio and an estimation for difference in proportions.
#' @import stats
#' @importFrom utils head tail


#' @param HR_cox Input. Hazard ratio obtained from Cox PH model.
#' @param HR_cox_CI Input. Confidence interval of the hazard ratio obtained from Cox PH model.
#' @param s1mix.chosen Input. Survival probabilities of Arm 1 estimated by KM at the chosen time.
#' @param s0mix.chosen Input. Survival probabilities of Arm 0 estimated by KM at the chosen time.
#' @param pi1.est Input. Proportion of short-term survivors (uncured proportion) in Arm 1.
#' @param pi0.est Input. Proportion of long-term survivors (uncured proportion) in Arm 0.
#' @param m Input (Optional). The polynomial order used in Taylor approximation.
#' \emph{m} is selected automatically if m = NULL (default).
#' @param cure.show Input (Optional). Indicate whether or not to show the estimated proportions of long-term survivors
#' The default is TRUE.
#'

#' @return \item{Adj.before.after}{Results before and after Cox-TEL adjustment: Cox HR, Cox-TEL HR, and Cox-TEL DP.}
#' @return \item{Proportion.LTS}{Estimated proportions of long-term survivors with confidence intervals.}
#'
#' @references
#' C.Y. Hsu, E.P.Y. Lin and Y Shyr. (2021). Development and Evaluation of a Method to
#' Correct Misinterpretation of Clinical Trial Results With Long-term Survival.
#' JAMA Oncol. doi:10.1001/jamaoncol.2021.0289.


#' @examples library(CoxTEL)
#' @examples s1mix.chosen <- c(0.48,0.39,0.35,0.35)
#' @examples s0mix.chosen <- c(0.36,0.28,0.26,0.25)
#' @examples pi1.est <- 0.65; pi0.est <- 0.75 # Proportion of short-term survivors
#' @examples HR_cox <- 0.71; HR_cox_CI <- c(0.51, 0.91)
#' @examples adjustment(HR_cox, HR_cox_CI, s1mix.chosen, s0mix.chosen, pi1.est, pi0.est)
#'

#' @export

adjustment <- function(HR_cox, HR_cox_CI, s1mix.chosen, s0mix.chosen, pi1.est, pi0.est, m = NULL,
                       cure.show = TRUE) {

  if (!isTRUE(all.equal(s0mix.chosen, sort(s0mix.chosen, decreasing = T))) |
      !isTRUE(all.equal(s1mix.chosen, sort(s1mix.chosen, decreasing = T)))) {
    stop("s1mix.chosen or/and s0mix.chosen inputed is NOT decreasing.")
  }
  if ( (tail(s0mix.chosen,1) < 1-pi0.est - 1E-8) |
       (tail(s1mix.chosen,1) < 1-pi1.est - 1E-8) ) {
    stop("All values in sjmix.chosen inputed must be larger than (1 - pij.est).")
  }

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


  HR_cal <- CF * HR_cox
  HR_cal_CI <- CF * HR_cox_CI

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

  RR <- round(c(HR_cox = HR_cox, HR_cox_CI = HR_cox_CI,
          HR_adj = HR_cal, HR_adj_CI = HR_cal_CI,
          DP_adj = pi0.est - pi1.est, DP_adj_CI = DP_CI), 3)
  names(RR) <- c("Cox_HR", "Cox_HR_CIL", "Cox_HR_CIU",
                 "CoxTEL_HR", "CoxTEL_HR_CIL", "CoxTEL_HR_CIU",
                 "CoxTEL_DP", "CoxTEL_DP_CIL", "CoxTEL_DP_CIU")

  if (cure.show) {
    prop.cure0 <- 1 - c(pi0.est, max(CI.pi0), min(CI.pi0))
    prop.cure1 <- 1 - c(pi1.est, max(CI.pi1), min(CI.pi1))
    prop.cure <- round(c(prop.cure0, prop.cure1), 3)
    names(prop.cure) <- c("Arm0", "Arm0_CIL", "Arm0_CIU", "Arm1", "Arm1_CIL", "Arm1_CIU")
    list(Adj.before.after = RR, Proportion.LTS = prop.cure, s1mix.chosen = s1mix.chosen, s0mix.chosen = s0mix.chosen)
  } else {
    list(Adj.before.after = RR, s1mix.chosen = s1mix.chosen, s0mix.chosen = s0mix.chosen)
  }

}
