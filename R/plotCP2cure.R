
#' Plot conditional probability of Pts belonging to long-term survivors
#' @description A figure function to show conditional probability of Pts belonging to long-term survivors.
#' @importFrom graphics polygon
#' @importFrom grDevices adjustcolor
#' @importFrom graphics legend lines par points axis

#' @param CT.object Input. Object of \emph{adjustment}.
#' @param timeunit Input. Time scale.
#' @param alpha.f Input. Alpha transparency in \emph{graphics::polygon}.
#' @param ... Input. Arguments in \emph{plot}.
#'
#' @details P(Pts belongs to long-term survivors | T>t) = (1 - pi)/smix(t).
#' The confidence interval are defined as ((1 - pi.L)/smix(t), (1 - pi.U)/smix(t)).


#' @examples library(CoxTEL)
#' @examples s1mix.chosen <- c(93, 83, 74, 66, 61, 56, 52, 49, 44, 42)/100
#' @examples s0mix.chosen <- c(86, 75, 62, 55, 48, 42, 40, 36, 33, 33)/100
#' @examples pi1.est <- 1 - 0.42; pi0.est <- 1 - 0.33 # Proportion of short-term survivors
#' @examples HR_cox <- 0.68; HR_cox_CI <- c(0.53, 0.87)
#' @examples CT.object <- adjustment(HR_cox, HR_cox_CI, s1mix.chosen, s0mix.chosen, pi1.est, pi0.est)
#' @examples plotCP2cure(CT.object, xlab="Months")
#'

#' @export


plotCP2cure <- function(CT.object, timeunit=6, alpha.f=0.3, ...) {

  CP2Cure <- function(cureprop, smix) {
    pi0 <- 1 - cureprop
    #pi0Pu <- c(1, smix) - (1 - pi0)
    #prob <- 1 - pi0Pu/(pi0Pu + (1 - pi0))
    prob <- (1 - pi0)/c(1, smix)
    prob <- ifelse(prob>1, 1, prob)
    prob <- ifelse(prob<0, 0, prob)
    prob
  }

  par(mfrow = c(1,1), mar = c(5,6.5,3,1))
  time <- c(0, 1:length(CT.object$s0mix.chosen))
  plot(time, CP2Cure(CT.object$Proportion.LTS[1], CT.object$s0mix.chosen), ylim=c(0,1),
       type="l",
       col="#00AFBB",
       ylab="Conditional probability of \n Pts belonging to long-term survivors",
       xaxt='n',
       cex.lab=1.5,
       cex.axis=1.5, ...)
  axis(1, at=time, labels=time*timeunit, cex.axis=1.5)
  points(time, CP2Cure(CT.object$Proportion.LTS[1], CT.object$s0mix.chosen), pch=16, col="#00AFBB")

  ld1 <- CP2Cure(CT.object$Proportion.LTS[2], CT.object$s0mix.chosen)
  ud1 <- CP2Cure(CT.object$Proportion.LTS[3], CT.object$s0mix.chosen)
  polygon(x = c(time, rev(time)), y = c(ld1, rev(ud1)),
          col = adjustcolor("#CCE5FF", alpha.f = alpha.f + 0.2),
          border = NA)

  time <- c(0, 1:length(CT.object$s1mix.chosen))
  lines(time, CP2Cure(CT.object$Proportion.LTS[4], CT.object$s1mix.chosen), col="#E7B800")
  points(time, CP2Cure(CT.object$Proportion.LTS[4], CT.object$s1mix.chosen), pch=16, col="#E7B800")

  ld1 <- CP2Cure(CT.object$Proportion.LTS[5], CT.object$s1mix.chosen)
  ud1 <- CP2Cure(CT.object$Proportion.LTS[6], CT.object$s1mix.chosen)
  polygon(x = c(time, rev(time)), y = c(ld1, rev(ud1)),
          col = adjustcolor("#FFE489", alpha.f = alpha.f),
          border = NA)

  legend("bottomright",
         legend = c("Arm_0", "Arm_1"),
         col=c("#00AFBB", "#E7B800"), pch=c(16,16), bty = "n", lwd=2, cex=1.5)


}
