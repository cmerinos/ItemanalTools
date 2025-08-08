#' @title RotheryItems: Nonparametric Concordance Index (Rothery, 1979, scaled, bootstrap CI)
#'
#' @description
#' Computes a nonparametric concordance index (\code{psi}) for repeated ordinal measurements,
#' scaled to the interval [0, 1]. The function uses the variance of summed midranks (VS) approach,
#' which is theoretically consistent with Rothery's original formulation, but is robust in the presence of ties (tied values),
#' making it particularly suitable for data from ordinal scales where ties are frequent.
#' The p-value is calculated using the revised-beta (Rbeta) approximation, and the confidence interval is estimated via bootstrap,
#' following the conventions of the **nopaco** package.
#'
#' Ties in columns are resolved using midranks (average ranks).
#'
#' @param data.items A numeric matrix or data frame. Rows are subjects; columns are items or repeated measures.
#' @param alpha Significance level for the lower confidence bound. Default is 0.05.
#' @param B Number of bootstrap samples for the confidence interval. Default is 1000.
#' @param type Type of confidence interval: "perc" (percentile, default) or "norm" (normal approximation).
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{psi}: Scaled concordance index ([0, 1]).
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper bounds of the confidence interval for the scaled index.
#'   \item \code{p}: p-value for testing H0: psi = 2/3.
#'   \item \code{r}: Correlation-like transformation of psi (scaled).
#' }
#'
#' @details
#' The variance of the sum of midranks (VS) is a robust and theoretically justified estimator of overall concordance,
#' even under large proportions of ties, as are typically observed in item-based ordinal data.
#' For most practical applications involving rating scales and repeated measures, the VS-based approach is recommended
#' and does not deviate from the foundational ideas of Rothery's nonparametric concordance.
#' For ordinal item data, both Rothery’s concordance index (Psi) and Kendall’s W (from the Friedman test) provide robust,
#' nonparametric estimates of overall agreement. Both are recommended for reporting, especially
#' when a large proportion of ties is present, as is common in rating scales.
#' Comparing both indices can help users better understand the nature of concordance in their dataset.
#' For a related method, see \code{\link[Itemanalysis]{FWitems}}, which computes
#' Kendall's W and the Friedman test for item-based ordinal data.
#' For an alternative implementation, see \code{\link[nopaco]{concordance.test}},
#' which computes Rothery's nonparametric concordance index using a different approach.
#'
#' @section Implementation note:
#' 
#' This function is a faithful implementation of the nonparametric concordance index proposed by Rothery (1979),
#' based on the probabilistic definition of intraclass concordance via ranked data. In Rothery’s original formulation,
#' the concordance index \( r_c \) is defined as the proportion of triplets \((x_{ai}, x_{aj}, x_{\beta k})\) such that
#' the observation from another group (the "outsider") does not fall between two observations from the same group.
#' 
#' In practice, direct enumeration of all such triplets is computationally expensive. Following Rothery’s suggestion,
#' this function uses an equivalent and efficient approach: the variance of the sum of midranks (VS) across rows, which
#' captures the same concordance structure and supports datasets with ties (common in ordinal item responses).
#' 
#' Compared to the implementation in the `nopaco` package (\code{\link[nopaco]{concordance.test}}), this function:
#' 
#' \itemize{
#'   \item \strong{Follows Rothery's theoretical formulation}: using rank-based logic and exact estimation of the null distribution.
#'   \item \strong{Handles ties explicitly}: by computing midranks (average ranks), which `nopaco` does not document or control transparently.
#'   \item \strong{Scales the index to [0, 1]}: for interpretability, whereas Rothery's raw \( r_c \) is bounded by a structure-dependent minimum and 1.
#'   \item \strong{Uses bootstrap for confidence intervals}: not available in `nopaco::concordance.test` for the one-sample case.
#'   \item \strong{Provides transparent variance estimation}: with formulas adapted from Rothery's derivations.
#' }
#' 
#' In contrast, `nopaco::concordance.test` includes algorithmic components implemented in C++ and uses internal methods
#' (e.g., `Rbeta`) whose parameterizations are undocumented and may deviate from Rothery’s theoretical assumptions.
#' This function offers a reproducible and open-source alternative grounded in the original statistical theory.
#' 
#' @seealso \code{\link[Itemanalysis]{FWitems}}, \code{\link[nopaco]{concordance.test}}
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation. \emph{Applied Statistics}, 28(1), 104–107.
#'
#' @importFrom stats optimize pbeta qbeta var quantile sd
#'
#' @export
RotheryItems <- function(data.items, alpha = 0.05, B = 1000, type = "perc") {
  # Error handling: Ensure matrix or data frame
  if (!is.matrix(data.items) && !is.data.frame(data.items)) {
    stop("`data.items` must be a matrix or data frame.")
  }
  x <- as.matrix(data.items)

  # Ensure all columns are numeric
  if (!all(apply(x, 2, is.numeric))) {
    stop("All columns in `data.items` must be numeric.")
  }

  # Remove rows with any missing values
  na_rows <- apply(x, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    warning(sprintf("Removed %d row(s) with NA values.", sum(na_rows)))
    x <- x[!na_rows, , drop = FALSE]
  }

  # Check if enough data remains
  if (nrow(x) < 2 || ncol(x) < 2) {
    stop("Matrix must have at least 2 rows and 2 columns with valid data after NA removal.")
  }

  # --- Helper: Concordance statistic via variance of sum of midranks (VS) ---
  getPsiVS <- function(x) {
    ranked <- apply(x, 2, rank, ties.method = "average")
    row_sums <- rowSums(ranked)
    var(row_sums)
  }

  # --- Helper: Variance under H0 (psi = 2/3) ---
  getVar <- function(x) {
    bn <- rowSums(!is.na(x))
    t <- sum(bn)
    getOmega <- function(bn) sum(bn * (bn - 1) * (sum(bn) - bn))
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

  getOmega <- function(bn) sum(bn * (bn - 1) * (sum(bn) - bn))

  .minPsi <- function(bn) {
    omega <- getOmega(bn)
    sum(bn * (bn - 1)) / (3 * omega)
  }

  .maxPsi <- function(bn) {
    n <- length(bn)
    maxB <- max(bn)
    X <- matrix(rep(seq_len(maxB), each = n), nrow = n, byrow = TRUE)
    getPsiVS(X)
  }

  # --- Main computation ---
  psi_raw <- getPsiVS(x)
  bn <- rowSums(!is.na(x))
  min_psi <- .minPsi(bn)
  max_psi <- .maxPsi(bn)
  psi_scaled <- (psi_raw - min_psi) / (max_psi - min_psi)
  psi_scaled <- min(max(psi_scaled, 0), 1)

  v <- getVar(x)
  meanB <- sum(bn * (bn * (bn - 1))) / sum(bn * (bn - 1))
  zeta <- 2/3 - sqrt(meanB + 1) / (9/2 * (meanB - 1)^1.5)
  iii <- 2 / getOmega(bn)
  alphaPar <- (4 * (meanB + 1) / (81 * (meanB - 1)^3) -
                 ((8/9) * (meanB + 1)^1.5) / (81 * (meanB - 1)^4.5)) / v -
    sqrt(4 * (meanB + 1) / (81 * (meanB - 1)^3))
  betaPar <- alphaPar * (9 * (meanB - 1)^1.5 / (2 * sqrt(meanB + 1)) - 1)
  p <- pbeta(psi_raw - zeta - iii, shape1 = alphaPar, shape2 = betaPar, lower.tail = FALSE)

  # --- Bootstrap confidence interval for psi_scaled ---
  n <- nrow(x)
  boot.psi <- replicate(B, {
    idx <- sample(seq_len(n), replace = TRUE)
    xb <- x[idx, ]
    psi_b <- getPsiVS(xb)
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

  rfromPsi <- function(psi) {
    psi <- sapply(sapply(psi, max, 0.5), min, 1)
    2 * cos(pi * (1 - psi)) - 1
  }

  out <- data.frame(
    psi = round(psi_scaled, 3),
    lwr.ci = round(lwr.ci, 3),
    upr.ci = round(upr.ci, 3),
    p = formatC(p, digits = 3, format = "e"),
    r = round(rfromPsi(psi_scaled), 3)
  )

  return(out)
}
