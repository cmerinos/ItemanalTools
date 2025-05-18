#' @title RotheryItems: Nonparametric Concordance Index (Rothery, 1979)
#'
#' @description
#' Computes a nonparametric measure of intraclass concordance for ordinal item responses,
#' using the Rbeta method from Rothery (1979), adapted from the nopaco package.
#'
#' @param data.items A matrix or data frame with ordinal item responses (rows = subjects, columns = items).
#' @param alpha Significance level for the lower confidence bound. Default is 0.05.
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{psi}: Concordance index.
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper confidence interval.
#'   \item \code{p}: p-value for H0: psi = 2/3.
#'   \item \code{r}: Equivalent correlation coefficient.
#' }
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation.
#' \emph{Applied Statistics}, 28(1), 104–107.
#'
#' @importFrom stats optimize pbeta qbeta
#' @export
RotheryItems <- function(data.items, alpha = 0.05) {
  # ----- Internal functions -----
  getPsi <- function(x) {
    x <- as.matrix(x)
    ranks <- apply(x, 2, rank)
    rowSums <- rowSums(ranks)
    var(rowSums)
  }

  getOmega <- function(bn) {
    sum(bn * (bn - 1) * (sum(bn) - bn))
  }

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

  .minPsi <- function(bn) {
    omega <- getOmega(bn)
    sum(bn * (bn - 1)) / (3 * omega)
  }

  .confEstimatorBeta <- function(beta, mu, p, targetValue, lower.tail, ...) {
    alpha <- mu * beta / (1 - mu)
    p_est <- pbeta(targetValue, shape1 = alpha, shape2 = beta, lower.tail = lower.tail)
    abs(p_est - p)
  }

  rfromPsi <- function(psi) {
    psi <- sapply(sapply(psi, max, 0.5), min, 1)
    2 * cos(pi * (1 - psi)) - 1
  }

  # ----- Input checks -----
  x <- as.matrix(data.items)
  if (any(is.na(x))) stop("Missing values are not allowed.")

  psi <- getPsi(x)
  v <- getVar(x)
  bn <- rowSums(!is.na(x))
  meanB <- sum(bn * (bn * (bn - 1))) / sum(bn * (bn - 1))
  zeta <- (2/3 - sqrt(meanB + 1) / (9/2 * (meanB - 1)^1.5))
  iii <- 2 / getOmega(bn)

  alphaPar <- (4 * (meanB + 1) / (81 * (meanB - 1)^3) -
                 ((8/9) * (meanB + 1)^1.5) / (81 * (meanB - 1)^4.5)) / v -
    sqrt(4 * (meanB + 1) / (81 * (meanB - 1)^3))

  betaPar <- alphaPar * (9 * (meanB - 1)^1.5 / (2 * sqrt(meanB + 1)) - 1)

  p <- pbeta(psi - zeta - iii, shape1 = alphaPar, shape2 = betaPar, lower.tail = FALSE)

  est <- optimize(interval = c(1, 1e8), mu = psi - zeta - iii, p = p,
                  lower.tail = TRUE, targetValue = 2/3 - zeta - iii, f = .confEstimatorBeta)

  betaEst <- est$minimum
  alphaEst <- (psi - zeta - iii) * betaEst / (1 - (psi - zeta - iii))
  ci.lower <- max(.minPsi(bn), qbeta(alpha, alphaEst, betaEst, lower.tail = TRUE) + zeta + iii)

  out <- data.frame(
    psi = round(psi, 3),
    lwr.ci = round(ci.lower, 3),
    upr.ci = 1,
    p = signif(p, 3),
    r = round(rfromPsi(psi), 3)
  )

  return(out)
}
