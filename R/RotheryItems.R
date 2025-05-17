#' @title RotheryItems: Nonparametric Intraclass Correlation (Rothery, 1979)
#'
#' @description
#' Computes a nonparametric measure of intraclass correlation (Rothery, 1979) based on ranks,
#' to assess agreement among items. Confidence intervals are estimated via bootstrap.
#'
#' @details
#' The Rothery index ranges from 0 (no agreement) to 1 (perfect agreement),
#' and is analogous in interpretation to the intraclass correlation coefficient (ICC1).
#' It is computed from the ratio of between-column variance in ranks over total variance.
#'
#' @param data.items A data frame or matrix with ordinal item responses (rows = subjects, columns = items).
#' @param ci Logical. If TRUE, computes confidence interval via bootstrap. Default is TRUE.
#' @param B Integer. Number of bootstrap replications. Default is 1000.
#' @param conf.level Confidence level for the interval. Default is 0.95.
#' @param type Character. Type of confidence interval: \code{"perc"} (percentile; default) or \code{"norm"} (normal approximation).
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{psi}: Rothery intraclass correlation.
#'   \item \code{lwr.ci}, \code{upr.ci}: lower and upper confidence limits (if \code{ci = TRUE}).
#'   \item \code{p}: approximate p-value testing the null hypothesis of no agreement (psi = 0).
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
#'@importFrom stats quantile qnorm sd
#'
#' @export
RotheryItems <- function(data.items, ci = TRUE, B = 1000, conf.level = 0.95, type = "perc") {
  if (!is.matrix(data.items) && !is.data.frame(data.items)) {
    stop("`data.items` must be a matrix or data frame.")
  }

  data <- as.matrix(data.items)
  if (any(is.na(data))) stop("Missing values are not allowed.")
  if (!type %in% c("perc", "norm")) stop("`type` must be 'perc' or 'norm'.")

  n <- nrow(data)
  k <- ncol(data)
  R <- matrix(rank(data), nrow = n, ncol = k)

  Rbar <- mean(R)
  Rj.bar <- colMeans(R)
  ST <- sum((R - Rbar)^2)
  SB <- n * sum((Rj.bar - Rbar)^2)

  psi <- SB / ST
  out <- data.frame(psi = round(psi, 3))

  if (ci) {
    boot.psi <- replicate(B, {
      idx <- sample(seq_len(n), replace = TRUE)
      Rb <- R[idx, ]
      Rbbar <- mean(Rb)
      Rjb.bar <- colMeans(Rb)
      STb <- sum((Rb - Rbbar)^2)
      SBb <- n * sum((Rjb.bar - Rbbar)^2)
      SBb / STb
    })

    alpha <- 1 - conf.level

    if (type == "perc") {
      ci.lwr <- quantile(boot.psi, probs = alpha / 2)
      ci.upr <- quantile(boot.psi, probs = 1 - alpha / 2)
    } else if (type == "norm") {
      m <- mean(boot.psi)
      s <- sd(boot.psi)
      z <- qnorm(1 - alpha / 2)
      ci.lwr <- m - z * s
      ci.upr <- m + z * s
    }

    p <- mean(boot.psi <= 0)

    out$lwr.ci <- round(ci.lwr, 3)
    out$upr.ci <- round(ci.upr, 3)
    out$p <- signif(p, 3)
  }

  return(out)
}
