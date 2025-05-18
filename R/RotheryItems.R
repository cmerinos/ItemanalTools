#' @title RotheryItems: Nonparametric Concordance Index (Rothery, 1979)
#'
#' @description
#' Computes a nonparametric concordance index (\code{psi}) for repeated ordinal measures,
#' using the revised-beta approach (Rbeta) based on Rothery (1979). Confidence intervals
#' and p-values are estimated following the nopaco package method.
#'
#' @param data.items A numeric matrix or data frame. Rows are subjects; columns are items or repeated measures.
#' @param alpha Significance level for the lower confidence bound. Default is 0.05.
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{psi}: Concordance index.
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper bounds of the confidence interval.
#'   \item \code{p}: p-value for testing H0: psi = 2/3.
#'   \item \code{r}: Correlation-like transformation of psi.
#' }
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation. \emph{Applied Statistics}, 28(1), 104–107.
#'
#' @importFrom stats optimize pbeta qbeta var
#' @export
RotheryItems <- function(data.items, alpha = 0.05) {

  # ----- Aux: get concordance statistic -----
  getPsi <- function(x) {
    x <- as.matrix(x)
    ranked <- apply(x, 1, rank)
    row_sums <- rowSums(ranked)
    var(row_sums)
  }

  # ----- Aux: get variance under H0 (psi = 2/3) -----
  getVar <- function(x) {
    x <- as.matrix(x)
    bn <- rowSums(!is.na(x))
    t <- sum(bn)
    omega <- getOmega(bn)
    cij <- 0
    for (i in 1:(length(bn) - 1)) {
      for (j in (i + 1):length(bn)) {
        cij <- cij + 2 * (bn[i] - 1) * bn[i] * (bn[j] - 1) * bn[j]
      }
    }
    cii <- sum((bn - 1) * bn * (bn + 3) * (t - bn))
    (cii - cij) / (45 * omega^2) * (t + 1)
  }

  # ----- Aux: get Omega constant -----
  getOmega <- function(bn) {
    sum(bn * (bn - 1) * (sum(bn) - bn))
  }

  # ----- Aux: minimum possible value of psi under the design -----
  .minPsi <- function(bn) {
    omega <- getOmega(bn)
    sum(bn * (bn - 1)) / (3 * omega)
  }

  # ----- Aux: optimization target -----
  .confEstimatorBeta <- function(beta, mu, p, targetValue, lower.tail) {
    alpha <- mu * beta / (1 - mu)
    p_est <- pbeta(targetValue, shape1 = alpha, shape2 = beta, lower.tail = lower.tail)
    abs(p_est - p)
  }

  # ----- Aux: r transformation -----
  rfromPsi <- function(psi) {
    psi <- sapply(sapply(psi, max, 0.5), min, 1)
    2 * cos(pi * (1 - psi)) - 1
  }

  # ---------- Main body ----------
  x <- as.matrix(data.items)
  if (any(is.na(x))) stop("Missing values are not allowed.")

  psi <- getPsi(x)
  v <- getVar(x)
  bn <- rowSums(!is.na(x))
  meanB <- sum(bn * (bn * (bn - 1))) / sum(bn * (bn - 1))

  zeta <- 2/3 - sqrt(meanB + 1) / (9/2 * (meanB - 1)^1.5)
  iii <- 2 / getOmega(bn)

  alphaPar <- (4 * (meanB + 1) / (81 * (meanB - 1)^3) -
                 ((8/9) * (meanB + 1)^1.5) / (81 * (meanB - 1)^4.5)) / v -
    sqrt(4 * (meanB + 1) / (81 * (meanB - 1)^3))

  betaPar <- alphaPar * (9 * (meanB - 1)^1.5 / (2 * sqrt(meanB + 1)) - 1)

  # p-value
  p <- pbeta(psi - zeta - iii, shape1 = alphaPar, shape2 = betaPar, lower.tail = FALSE)

  # confidence lower bound
  min_psi <- .minPsi(bn)
  if (is.nan(alphaPar) || is.nan(betaPar) || alphaPar <= 0 || betaPar <= 0) {
    ci.lower <- NA
    p <- NA
  } else {
    est <- optimize(interval = c(1, 1e8), mu = psi - zeta - iii, p = p,
                    lower.tail = TRUE, targetValue = 2/3 - zeta - iii,
                    f = .confEstimatorBeta)
    betaEst <- est$minimum
    alphaEst <- (psi - zeta - iii) * betaEst / (1 - (psi - zeta - iii))
    ci.lower <- max(min_psi,
                    qbeta(alpha, alphaEst, betaEst, lower.tail = TRUE) + zeta + iii)
  }

  # Output
  out <- data.frame(
    psi = round(psi, 3),
    lwr.ci = round(ci.lower, 3),
    upr.ci = 1,
    p = signif(p, 3),
    r = round(rfromPsi(psi), 3)
  )

  return(out)
}
