#' Sample Size and Power Comparisons of Placebo Crossover and Standard Two-Arm Designs with Two Periods of Follow-up
#'
#' Considering two periods of participant follow-up with associated levels of vaccine efficacy \eqn{VE_1} and \eqn{VE_2}, \code{plaXpower} computes for both a placebo crossover and standard non-crossover design (i) power to reject the null hypothesis \eqn{H_0: VE_1 = VE_2} in favor of the one-sided alternative hypothesis of waning VE, i.e., \eqn{H_1: VE_1 > VE_2}, (ii) power to reject the null hypothesis \eqn{H_0: VE_2 \ge 0} in favor of the one-sided alternative hypothesis of vaccine harm, i.e., \eqn{H_1: VE_2 < 0}, and (iii) for each hypothesis test, the sample size ratio (placebo crossover/standard) required to achieve the same power. The event of interest is assumed to be rare, with the event count following a Poisson distribution in each arm and time period. The vaccine administered in period 2 is assumed to have equal VE as in period 1.
#'
#' @aliases print.plaXpower
#' @param theta1 a numeric value specifying the Poisson mean event count in the original placebo arm in period 1 in the standard design. This mean count is assumed in all scenarios described by \code{ve1} and \code{ve2}.
#' @param theta2 a numeric value specifying the Poisson mean event count in the original placebo arm in period 2 in the standard design. This mean count is assumed in all scenarios described by \code{ve1} and \code{ve2}.
#' @param ve1 a numeric vector specifying different scenarios of vaccine efficacy in period 1 defined as 1 minus the ratio (original vaccine/original placebo) of mean event counts
#' @param ve2 a numeric vector specifying different scenarios of vaccine efficacy in period 2 defined as 1 minus the ratio (original vaccine/original placebo) of mean event counts. Vectors \code{ve1} and \code{ve2} must be of equal length, with matching components of \code{ve1} and \code{ve2} treated as a single trial scenario.
#' @param method a character string specifying the computational method. The default option is \code{"analytical"}, which assumes that the mean event counts and VE levels in both arms and time periods and in both designs are known. The alternative option is \code{"simulation"}, which first randomly samples event counts from Poisson distributions defined by \code{theta1}, \code{theta2}, \code{ve1}, and \code{ve2}, and then uses the data to estimate VE and the variance. Only the first character of the method is necessary.
#' @param alpha one-sided nominal significance level (0.025 by default)
#' @param iter the number of Monte-Carlo iterations applied for the \code{"simulation"} method (\eqn{10^3} by default)
#' @param seed a seed of the random number generator supplied to \code{set.seed} for reproducibility
#'
#' @return
#' An object of class \code{plaXpower}, which is a data frame with the following columns:
#' \itemize{
#' \item \code{ve1}: repeats the input argument \code{ve1}
#' \item \code{ve2}: repeats the input argument \code{ve2}
#' \item \code{ssRatioWaneXstd}: the sample size ratio (placebo crossover/standard) that ensures the same power to detect waning VE
#' \item \code{ssRatioHarmXstd}: the sample size ratio (placebo crossover/standard) that ensures the same power to detect vaccine harm
#' \item \code{powerWaneX}: power to detect waning VE in the placebo crossover design
#' \item \code{powerWaneStd}: power to detect waning VE in the standard non-crossover design
#' \item \code{powerHarmX}: power to detect vaccine harm in the placebo crossover design
#' \item \code{powerHarmStd}: power to detect vaccine harm in the standard non-crossover design
#' }
#'
#' @examples
#' # analytical method
#' tab <- plaXpower(theta1=100,
#'                  theta2=100,
#'                  ve1=rep(0.8, 5),
#'                  ve2=c(0.8, 0.5, 0.2, -0.1, -0.4))
#'
#' # pretty printing of 'tab' with 2 decimal places (default)
#' tab
#'
#' # pretty printing of 'tab' with 3 decimal places
#' print(tab, digits=3)
#'
#' # unrounded real numbers in 'tab' for, e.g., plotting of a power curve
#' tab$powerWaneX
#'
#' # simulation method
#' tab <- plaXpower(theta1=100,
#'                  theta2=100,
#'                  ve1=rep(0.8, 5),
#'                  ve2=c(0.8, 0.5, 0.2, -0.1, -0.4),
#'                  method="s",
#'                  iter=100,
#'                  seed=19283)
#'
#' @export
plaXpower <- function(theta1, theta2, ve1, ve2, method=c("analytical", "simulation"), alpha=0.025, iter=1e3, seed=NULL){
  if (length(theta1)>1 || length(theta2)>1){ stop("Placebo mean event counts 'theta1' and 'theta2' must be scalars.") }
  if (length(ve1)!=length(ve2)){ stop("Vectors 've1' and 've2' must be of equal length.") }

  method <- match.arg(method)

  qZ<- qnorm(1 - alpha)
  out <- data.frame(ve1=ve1, ve2=ve2)

  if (method=="analytical"){
    # the sample size of the crossover design is greater by factor 'varRatioWaneXstd' or 'varRatioHarmXstd'
    # to achieve the same power as in the standard design to detect waning VE or harm, respectively

    # var(log[(1 - \hat{VEX}_2) / (1 - \hat{VEX}_1)])
    varWaneX <- 1/(theta2 * (1-ve1)) +  1/(theta2 * (1-ve2))
    # var(log[(1 - \hat{VES}_2) / (1 - \hat{VES}_1)])
    varWaneStd <- 1/theta1 + 1/(theta1 * (1-ve1)) + 1/theta2 + 1/(theta2 * (1-ve2))
    varRatioWaneXstd <- varWaneX / varWaneStd

    # var(log(1 - \hat{VEX}_2))
    varHarmX <- 1/theta1 + 1/(theta1 * (1-ve1)) + 1/(theta2 * (1-ve1)) + 1/(theta2 * (1-ve2))
    # var(log(1 - \hat{VES}_2))
    varHarmStd <- 1/theta2 + 1/(theta2 * (1-ve2))
    varRatioHarmXstd <- varHarmX / varHarmStd

    # mean of asymptotically normal Wald test statistics with variance 1
    # test of H0: VE1 = VE2
    muWaneX <- (log(1-ve2) - log(1-ve1)) / sqrt(varWaneX)
    muWaneStd <- (log(1-ve2) - log(1-ve1)) / sqrt(varWaneStd)

    # test of H0: VE2 >= 0
    muHarmX <- log(1-ve2) / sqrt(varHarmX)
    muHarmStd <- log(1-ve2) / sqrt(varHarmStd)

    # power to detect waning VE, i.e., VE2 < VE1 (1-sided test)
    powerWaneX <- pnorm(muWaneX - qZ)
    powerWaneStd <- pnorm(muWaneStd - qZ)

    # power to detect harm, i.e., VE2 < 0 (1-sided test)
    powerHarmX <- pnorm(muHarmX - qZ)
    powerHarmStd <- pnorm(muHarmStd - qZ)

    out <- data.frame(out, ssRatioWaneXstd=varRatioWaneXstd, ssRatioHarmXstd=varRatioHarmXstd,
                      powerWaneX, powerWaneStd, powerHarmX, powerHarmStd)

  } else {

    if (!is.null(seed)){ set.seed(seed) }

    out2 <- lapply(1:length(ve1), function(i, theta1, theta2, ve1, ve2, iter, qZ){
      # first, simulate case counts
      Y01 <- rpois(iter, theta1)
      Y11 <- rpois(iter, theta1 * (1-ve1[i]))
      Y12 <- rpois(iter, theta2 * (1-ve2[i]))
      Y02 <- rpois(iter, theta2)
      X02 <- rpois(iter, theta2 * (1-ve1[i]))

      # replace zeros with a small positive number
      # to avoid a zero in the denominator in the test statistics
      zeroFudge <- 0.05
      Y01[Y01==0] <- zeroFudge; Y11[Y11==0] <- zeroFudge
      Y12[Y12==0] <- zeroFudge; Y02[Y02==0] <- zeroFudge
      X02[X02==0] <- zeroFudge

      # the sample size of the crossover design is greater by factor 'meanVarRatioWaneXstd' or 'meanVarRatioHarmXstd'
      # to achieve the same power as in the standard design to detect waning VE or harm, respectively
      meanVarRatioWaneXstd <- mean((1/Y12 + 1/X02) / (1/Y12 + 1/Y02 + 1/Y11 + 1/Y01))
      meanVarRatioHarmXstd <- mean((1/Y11 + 1/Y01 + 1/Y12 + 1/X02) / (1/Y12 + 1/Y02))

      # Wald test statistics for the test of H0: VE1 = VE2
      waldWaneX <- log(Y12/X02) / sqrt(1/Y12 + 1/X02)
      waldWaneStd <- (log(Y12/Y02) - log(Y11/Y01)) / sqrt(1/Y12 + 1/Y02 + 1/Y11 + 1/Y01)

      # Wald test statistics for the test of H0: VE2 >= 0
      waldHarmX <- (log(Y11/Y01) + log(Y12/X02)) / sqrt(1/Y11 + 1/Y01 + 1/Y12 + 1/X02)
      waldHarmStd <- log(Y12/Y02) / sqrt(1/Y12 + 1/Y02)

      # power to detect waning VE, i.e., VE2 < VE1 (1-sided test)
      powerWaneX <- mean(waldWaneX>qZ)
      powerWaneStd <- mean(waldWaneStd>qZ)

      # power to detect harm, i.e., VE2 < 0 (1-sided test)
      powerHarmX <- mean(waldHarmX>qZ)
      powerHarmStd <- mean(waldHarmStd>qZ)

      return(c(meanVarRatioWaneXstd, meanVarRatioHarmXstd, powerWaneX, powerWaneStd, powerHarmX, powerHarmStd))
    }, theta1=theta1, theta2=theta2, ve1=ve1, ve2=ve2, iter=iter, qZ=qZ)

    out <- data.frame(out, do.call(rbind, out2))
    colnames(out)[-(1:2)] <- c("ssRatioWaneXstd", "ssRatioHarmXstd", "powerWaneX", "powerWaneStd", "powerHarmX", "powerHarmStd")
  }

  class(out) <- c("plaXpower", "data.frame")
  return(out)
}

#' @rdname plaXpower
#' @param x an object of class \code{plaXpower} from a call of \code{plaXpower}
#' @param digits the number of decimal places to use for printing (2 by default)
#' @param ... further arguments passed to or from other methods
#' @export
print.plaXpower <- function(x, digits=2, ...){
  class(x) <- "data.frame"
  x$ve1 <- paste0(x$ve1 * 100, "%")
  x$ve2 <- paste0(x$ve2 * 100, "%")
  x[, -(1:2)] <- format(round(x[, -(1:2)], digits=digits), nsmall=digits, scientific=FALSE)
  print(x)
}
