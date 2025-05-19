#' @title RotheryItems: Nonparametric Concordance Index (Rothery, 1979, scaled)
#'
#' @description
#' Computes a nonparametric concordance index (\code{psi}) for repeated ordinal measures,
#' scaled to [0,1], using the revised-beta approach (Rbeta) based on Rothery (1979) and the
#' nopaco package. Confidence intervals and p-values are estimated using the Rbeta method.
#'
#' @param data.items A numeric matrix or data frame. Rows are subjects; columns are items or repeated measures.
#' @param alpha Significance level for the lower confidence bound. Default is 0.05.
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{psi}: Scaled concordance index ([0,1]).
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper bounds of the confidence interval for the scaled index.
#'   \item \code{p}: p-value for testing H0: psi = 2/3.
#'   \item \code{r}: Correlation-like transformation of psi (scaled).
#' }
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation. \emph{Applied Statistics}, 28(1), 104–107.
#'
#' @importFrom stats optimize pbeta qbeta rank var
#' @export
RotheryItems <- function(data.items, alpha = 0.05) {

  # ----- Aux: get concordance statistic -----
  getPsi <- function(x) {
    x <- as.matrix(x)
    ranked <- apply(x, 2, rank)
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

  # ----- Aux: maximum possible value of psi under the design -----
  .maxPsi <- function(bn) {
    # Perfect concordance: all items/subjects have identical ranks
    n <- length(bn)
    maxB <- max(bn)
    X <- matrix(rep(seq_len(maxB), each = n), nrow = n, byrow = TRUE)
    getPsi(X)
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

  psi_raw <- getPsi(x)
  v <- getVar(x)
  bn <- rowSums(!is.na(x))
  meanB <- sum(bn * (bn * (bn - 1))) / sum(bn * (bn - 1))

  min_psi <- .minPsi(bn)
  max_psi <- .maxPsi(bn)
  # Avoid zero division
  psi_scaled <- (psi_raw - min_psi) / (max_psi - min_psi)
  psi_scaled <- min(max(psi_scaled, 0), 1)

  zeta <- 2/3 - sqrt(meanB + 1) / (9/2 * (meanB - 1)^1.5)
  iii <- 2 / getOmega(bn)

  alphaPar <- (4 * (meanB + 1) / (81 * (meanB - 1)^3) -
                 ((8/9) * (meanB + 1)^1.5) / (81 * (meanB - 1)^4.5)) / v -
    sqrt(4 * (meanB + 1) / (81 * (meanB - 1)^3))

  betaPar <- alphaPar * (9 * (meanB - 1)^1.5 / (2 * sqrt(meanB + 1)) - 1)

  # p-value for unscaled psi
  p <- pbeta(psi_raw - zeta - iii, shape1 = alphaPar, shape2 = betaPar, lower.tail = FALSE)

  # confidence lower bound for unscaled psi
  if (is.nan(alphaPar) || is.nan(betaPar) || alphaPar <= 0 || betaPar <= 0) {
    ci.lower.raw <- NA
    p <- NA
  } else {
    est <- optimize(interval = c(1, 1e8), mu = psi_raw - zeta - iii, p = p,
                    lower.tail = TRUE, targetValue = 2/3 - zeta - iii,
                    f = .confEstimatorBeta)
    betaEst <- est$minimum
    alphaEst <- (psi_raw - zeta - iii) * betaEst / (1 - (psi_raw - zeta - iii))
    ci.lower.raw <- max(min_psi,
                        qbeta(alpha, alphaEst, betaEst, lower.tail = TRUE) + zeta + iii)
  }

  # Scale the lower bound as well
  lwr.ci <- (ci.lower.raw - min_psi) / (max_psi - min_psi)
  lwr.ci <- min(max(lwr.ci, 0), 1)

  # Output
  out <- data.frame(
    psi = round(psi_scaled, 3),
    lwr.ci = round(lwr.ci, 3),
    upr.ci = 1,
    p = signif(p, 3),
    r = round(rfromPsi(psi_scaled), 3)
  )

  return(out)
}
