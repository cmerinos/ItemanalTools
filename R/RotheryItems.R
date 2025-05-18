#' @title RotheryItems: Nonparametric Concordance Index (Rothery, 1979)
#'
#' @description
#' Computes the nonparametric concordance statistic \eqn{\psi} from Rothery (1979),
#' using the revised-beta method to estimate a lower confidence bound and p-value.
#'
#' @param data.items A numeric matrix or data.frame. Rows are subjects, columns are items or replications.
#' @param alpha Significance level for the confidence interval. Default is 0.05.
#'
#' @return A data.frame with:
#' \itemize{
#'   \item \code{psi}: Concordance statistic (unscaled).
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper bounds of the confidence interval.
#'   \item \code{p}: One-sided p-value for H0: psi = 2/3.
#'   \item \code{r}: Correlation-like transformation of \eqn{\psi}.
#' }
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation. *Applied Statistics*, 28(1), 104–107.
#'
#' @importFrom stats var qbeta pbeta optimize
#' @export
RotheryItems <- function(data.items, alpha = 0.05) {

  # --- Internal functions ---

  getPsi <- function(x) {
    x <- as.matrix(x)
    R <- apply(x, 2, rank)  # rankea por ítem (columna)
    S <- rowSums(R)
    var(S)
  }

  getOmega <- function(bn) {
    sum(bn * (bn - 1) * (sum(bn) - bn))
  }

  getVar <- function(x) {
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
    bn <- sort(bn, decreasing = TRUE)
    maxB <- max(bn)
    n <- length(bn)
    Q <- R <- matrix(seq_len(maxB * n), nrow = n)
    for (i in seq_along(bn)) {
      R[i, -seq_len(bn[i])] <- NA
    }
    Q[] <- rank(R)
    Q[is.na(R)] <- NA
    getPsi(Q)
  }

  .confEstimatorBeta <- function(x, mu, p, lower.tail, targetValue) {
    b <- x
    a <- mu * b / (1 - mu)
    (pbeta(targetValue, shape1 = a, shape2 = b, lower.tail = lower.tail) - p)^2
  }

  rfromPsi <- function(psi) {
    psi <- pmax(pmin(psi, 1), 0)
    2 * cos(pi * (1 - psi)) - 1
  }

  # --- Main computation ---

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

  if (is.nan(alphaPar) || is.nan(betaPar) || alphaPar <= 0 || betaPar <= 0) {
    ci.lower <- NA
    p <- NA
  } else {
    p <- pbeta(psi - zeta - iii, shape1 = alphaPar, shape2 = betaPar, lower.tail = FALSE)

    est <- optimize(interval = c(1, 1e8), mu = psi - zeta - iii, p = p,
                    lower.tail = TRUE, targetValue = 2/3 - zeta - iii,
                    f = .confEstimatorBeta)
    betaEst <- est$minimum
    alphaEst <- (psi - zeta - iii) * betaEst / (1 - (psi - zeta - iii))
    ci.lower <- max(.minPsi(bn), qbeta(alpha, alphaEst, betaEst, lower.tail = TRUE) + zeta + iii)
  }

  return(data.frame(
    psi = round(psi, 3),
    lwr.ci = round(ci.lower, 3),
    upr.ci = 1,
    p = signif(p, 3),
    r = round(rfromPsi(psi), 3)
  ))
}
