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

  # Rank each row (subject) across items
  ranked <- t(apply(data, 1, rank))  # row-wise ranking

  # Compute sum of ranks per subject
  R_sum <- rowSums(ranked)
  psi <- var(R_sum) / var(sort(R_sum))  # normalized variance

  out <- data.frame(psi = round(psi, 3))

  if (ci) {
    boot.psi <- replicate(B, {
      idx <- sample(seq_len(n), replace = TRUE)
      boot_data <- data[idx, , drop = FALSE]
      boot_rank <- t(apply(boot_data, 1, rank))
      boot_R_sum <- rowSums(boot_rank)
      var(boot_R_sum) / var(sort(boot_R_sum))
    })

    alpha <- 1 - conf.level
    out$lwr.ci <- round(quantile(boot.psi, alpha / 2), 3)
    out$upr.ci <- round(quantile(boot.psi, 1 - alpha / 2), 3)
    out$p <- signif(mean(boot.psi <= 0), 3)
  }

  return(out)
}
