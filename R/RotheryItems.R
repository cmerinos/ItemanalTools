#' @title RotheryItems: Concordance Index for Repeated Measures (Rothery, 1979)
#'
#' @description
#' Computes a nonparametric measure of concordance across items (columns) for each subject (row),
#' following the method of Rothery (1979). Confidence intervals and p-values are estimated by bootstrap.
#'
#' @details
#' The Rothery index (\code{psi}) ranges from 0 (no agreement) to 1 (perfect concordance),
#' based on the variance of row-wise rank sums. This function reproduces the behavior of
#' \code{nopaco::concordance.test(x)} with a single matrix input.
#'
#' @param data.items A matrix or data frame where rows are subjects and columns are ordinal items (e.g., Likert-type).
#' @param ci Logical. If TRUE (default), computes a confidence interval via bootstrap.
#' @param B Number of bootstrap samples. Default is 1000.
#' @param conf.level Confidence level for the interval. Default is 0.95.
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{psi}: Nonparametric concordance index.
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper confidence bounds (if \code{ci = TRUE}).
#'   \item \code{p}: Bootstrap-based p-value for H0: psi = 0.
#' }
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation.
#' \emph{Applied Statistics}, 28(1), 104–107.
#'
#' @examples
#' set.seed(123)
#' dat <- matrix(sample(1:5, 100, replace = TRUE), ncol = 5)
#' RotheryItems(dat)
#'
#' @importFrom stats quantile sd
#' @export
RotheryItems <- function(data.items, ci = TRUE, B = 1000, conf.level = 0.95) {
  if (!is.matrix(data.items) && !is.data.frame(data.items)) {
    stop("`data.items` must be a matrix or data frame.")
  }

  data <- as.matrix(data.items)
  if (any(is.na(data))) stop("Missing values are not allowed.")

  n <- nrow(data)
  k <- ncol(data)

  # Step 1: rank each row (subject)
  ranked <- t(apply(data, 1, rank))  # Each row is ranked across items

  # Step 2: compute row sums of ranks
  R_sum <- rowSums(ranked)
  Rbar <- mean(R_sum)

  # Step 3: compute psi
  ST <- sum((R_sum - Rbar)^2)
  maxST <- (k^2 * (n^2 - 1)) / 12  # max variance of row sums under perfect agreement
  psi <- ST / maxST

  out <- data.frame(psi = round(psi, 3))

  if (ci) {
    boot.psi <- replicate(B, {
      idx <- sample(seq_len(n), replace = TRUE)
      boot_data <- data[idx, , drop = FALSE]
      boot_rank <- t(apply(boot_data, 1, rank))
      boot_R_sum <- rowSums(boot_rank)
      STb <- sum((boot_R_sum - mean(boot_R_sum))^2)
      STb / maxST
    })

    alpha <- 1 - conf.level
    ci.lwr <- quantile(boot.psi, probs = alpha / 2)
    ci.upr <- quantile(boot.psi, probs = 1 - alpha / 2)
    p <- mean(boot.psi <= 0)

    out$lwr.ci <- round(ci.lwr, 3)
    out$upr.ci <- round(ci.upr, 3)
    out$p <- signif(p, 3)
  }

  return(out)
}
