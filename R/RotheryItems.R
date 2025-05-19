#' @title RotheryItems: Nonparametric Concordance Index (Rothery, 1979, scaled, bootstrap IC)
#'
#' @description
#' Computes a nonparametric concordance index (\code{psi}) for repeated ordinal measures,
#' scaled to [0,1], using the revised-beta approach (Rbeta) for p-value and
#' bootstrap for the confidence interval, following the nopaco package convention.
#'
#' Empates (ties) en las columnas se resuelven usando rangos promedio (midranks).
#'
#' @param data.items A numeric matrix or data frame. Rows are subjects; columns are items or repeated measures.
#' @param alpha Significance level for the lower confidence bound. Default is 0.05.
#' @param B Number of bootstrap samples for confidence interval. Default is 1000.
#' @param type Type of confidence interval: "perc" (percentile, default) or "norm" (normal approx.).
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
#' @importFrom stats optimize pbeta qbeta rank var quantile sd
#' @export
RotheryItems <- function(data.items, alpha = 0.05, B = 1000, type = "perc") {

  # ----- Aux: get concordance statistic -----
  getPsi <- function(x) {
    x <- as.matrix(x)
    ranked <- apply(x, 2, rank, ties.method = "average")
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

  getOmega <- function(bn) {
    sum(bn * (bn - 1) * (sum(bn) - bn))
  }

  .minPsi <- function(bn) {
    omega <- getOmega(bn)
    sum(bn * (bn - 1)) / (3 * omega)
  }

  .maxPsi <- function(bn) {
    n <- length(bn)
    maxB <- max(bn)
    X <- matrix(rep(seq_len(maxB), each = n), nrow = n, byrow = TRUE)
    getPsi(X)
  }

  .confEstimatorBeta <- function(beta, mu, p, targetValue, lower.tail) {
    alpha <- mu * beta / (1 - mu)
    p_est <- pbeta(targetValue, shape1 = alpha, shape2 = beta, lower.tail = lower.tail)
    abs(p_est - p)
  }

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
  psi_scaled <- (psi_raw - min_psi) / (max_psi - min_psi)
  psi_scaled <- min(max(psi_scaled, 0), 1)

  zeta <- 2/3 - sqrt(meanB + 1) / (9/2 * (meanB - 1)^1.5)
  iii <- 2 / getOmega(bn)

  alphaPar <- (4 * (meanB + 1) / (81 * (meanB - 1)^3) -
                 ((8/9) * (meanB + 1)^1.5) / (81 * (meanB - 1)^4.5)) / v -
    sqrt(4 * (meanB + 1) / (81 * (meanB - 1)^3))

  betaPar <- alphaPar * (9 * (meanB - 1)^1.5 / (2 * sqrt(meanB + 1)) - 1)

  p <- pbeta(psi_raw - zeta - iii, shape1 = alphaPar, shape2 = betaPar, lower.tail = FALSE)

  # ---- Bootstrap confidence interval for psi_scaled ----
  n <- nrow(x)
  boot.psi <- replicate(B, {
    idx <- sample(seq_len(n), replace = TRUE)
    xb <- x[idx, ]
    psi_b <- getPsi(xb)
    # Use fixed min/max psi for scaling
    psi_scaled_b <- (psi_b - min_psi) / (max_psi - min_psi)
    min(max(psi_scaled_b, 0), 1)
  })

  alpha2 <- alpha / 2

  if (type == "norm") {
    lwr.ci <- mean(boot.psi) - qnorm(1 - alpha2) * sd(boot.psi)
    upr.ci <- mean(boot.psi) + qnorm(1 - alpha2) * sd(boot.psi)
  } else {  # percentile
    lwr.ci <- quantile(boot.psi, probs = alpha2, na.rm = TRUE)
    upr.ci <- quantile(boot.psi, probs = 1 - alpha2, na.rm = TRUE)
  }
  lwr.ci <- min(max(lwr.ci, 0), 1)
  upr.ci <- min(max(upr.ci, 0), 1)

  out <- data.frame(
    psi = round(psi_scaled, 3),
    lwr.ci = round(lwr.ci, 3),
    upr.ci = round(upr.ci, 3),
    p = formatC(p, digits = 3, format = "e"),
    r = round(rfromPsi(psi_scaled), 3)
  )

  return(out)
}
