

#' Piecewise regression criteria
#' @description A function to determine if a KM curve is matured.
#' @importFrom stats lm

#' @param smix.chosen Input. Extracted survival probabilities.
#' @param timepoints Input. Time points of the extracted survival probabilities.
#' @param slope.relative.change.threshold Input. A threshold for the relative change of the slopes of the last segment and the first segment.
#' @param ratio.long.short.threshold Input. A threshold for the ratio of the lengths of the last segment to the first two segments.
#'
#' @return \item{knots}{knots.}
#' @return \item{slop.flatline.with.95CI}{the slop of the KM tail with 95\% CI.}
#' @return \item{slope.relative.change}{the relative change of the slopes of the last segment and the first segment.}
#' @return \item{ratio.long.short}{the ratio of the lengths of the last segment to the first two segments.}
#' @return \item{sig}{TRUE if matured}
#' @return \item{PReg.fit}{the fitted results of piecewise regression.}

#' @references
#' E.P.Y. Lin, C.Y. Hsu, J.F. Chiou et al. (2022). Cox Proportional
#' Hazard Ratios Overestimate Survival Benefit of Immune Checkpoint Inhibitors (ICI):
#' Cox-TEL Adjustment and Meta-analyses of PD-L1 Expression and ICI Survival Benefit.
#' J Thorac Oncol. 17, 1365-1374. doi:10.1016/j.jtho.2022.08.010.
#'
#' E.P.Y. Lin, C.Y. Hsu, L. Berry, P. Bunn, Y. Shyr (2022).
#' Analysis of Cancer Survival associated with Immune Checkpoint Inhibitors after
#' Statistical Adjustment: A Systematic Review and Meta-analyses. JAMA Network Open.
#' 5(8):e2227211. doi:10.1001/jamanetworkopen.2022.27211.


#' @examples library(CoxTEL)
#' @examples s1mix.chosen <- c(39, 25, 18.2, 16, 14, 12.4, 11, 11)/100
#' @examples s0mix.chosen <- c(32, 18,  9.9,  8,  4,    3,  1,  1)/100
#' @examples PRC1 <- PRcriteria(smix.chosen=s1mix.chosen)
#' @examples PRC0 <- PRcriteria(smix.chosen=s0mix.chosen)
#' @examples PRC1$sig
#' @examples PRC0$sig
#'

#' @export


PRcriteria <- function(smix.chosen, timepoints = NULL,
                       slope.relative.change.threshold = 0.7, ratio.long.short.threshold = 1/3) {

  if(is.null(timepoints)) {
    x <- 0:length(smix.chosen)
  } else {
    x <- c(0, timepoints)
  }

  y <- c(1, smix.chosen)

  breaks <- x[-c(1,length(x))]
  mse <- matrix(NA, nrow=choose(length(breaks), 2), ncol=3)
  k <- 1
  for(i in 1:(length(breaks)-1) ){
    for (j in (i+1):length(breaks) ) {
      bp1 <- (x - breaks[i]) * (x>=breaks[i])
      bp2 <- (x - breaks[j]) * (x>=breaks[j])
      piecewise1 <- lm(y ~ x + bp1 + bp2)
      mse[k,] <- as.vector(c(i, j, unlist(summary(piecewise1)[6])))
      k <- k + 1
    }
  }
  mse.min.idx <- mse[which.min(mse[,3]),1:2]

  bp1 <- (x - breaks[mse.min.idx[1]]) * (x>=breaks[mse.min.idx[1]])
  bp2 <- (x - breaks[mse.min.idx[2]]) * (x>=breaks[mse.min.idx[2]])
  fit1 <- lm(y ~ x + bp1 + bp2)
  fit1.summary <- summary(fit1)

  slop.flatline <- sum(fit1.summary$coefficients[2:4,1])
  se.slop.flatline <- sqrt(c(t(c(1,1,1)) %*% fit1.summary$cov.unscaled[2:4,2:4] %*% c(1,1,1)))
  slop.flatline.with.CI <- round(slop.flatline + c(0, -1, 1) * 1.96 * se.slop.flatline, 3)

  cover.zero <- (min(slop.flatline.with.CI[2:3]) < 0 & max(slop.flatline.with.CI[2:3]) > 0)
  slope.relative.change <- as.vector(sum(fit1$coefficients[3:4])/abs(fit1$coefficients[2]))

  ratio <- round((tail(x,1) - x[mse.min.idx[2]+1])/(x[mse.min.idx[2]+1] - x[1]), 3)

  sig <- (cover.zero & (slope.relative.change > slope.relative.change.threshold) & ratio > ratio.long.short.threshold)

  list(knots = c(x[mse.min.idx[1]+1], x[mse.min.idx[2]+1]),
       slop.flatline.with.95CI=slop.flatline.with.CI,
       slope.relative.change = slope.relative.change,
       ratio.long.short = ratio,
       sig = sig,
       PReg.fit=fit1)
}

