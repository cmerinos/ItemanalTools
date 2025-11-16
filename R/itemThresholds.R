#' @title Item Thresholds and equidistance (single item)
#'
#' @description
#' Estimates thresholds (cutpoints) for an ordinal Likert-type item
#' using a cumulative link model (CLM).
#' Provides estimates, confidence intervals, differences between
#' consecutive thresholds, z-tests for those differences, and
#' indices/tests of threshold equidistance (Spratto-type index).
#'
#' @details
#' This function fits an ordinal regression model without predictors
#' (\code{ordinal::clm}) and extracts threshold estimates.
#' The differences between consecutive thresholds are directly
#' interpretable as standardized effect sizes on the latent variable scale,
#' especially when the \code{link = "probit"} option is used
#' (default). Larger differences indicate more distinct response
#' categories, while small differences suggest redundancy.
#'
#' In addition, it evaluates the equidistance of thresholds. Let
#' \eqn{\Delta_i = \tau_{i+1} - \tau_i} be the difference between
#' consecutive thresholds. The Equidistance Index (Spratto) is
#' computed as
#' \deqn{
#'   EI = \sqrt{\frac{1}{K-2}\sum_{i=1}^{K-2}(\Delta_i - \bar{\Delta})^2},
#' }
#' where \eqn{K} is the number of response categories and
#' \eqn{\bar{\Delta}} is the mean of the \eqn{\Delta_i}.
#'
#' Local deviations from equidistance are assessed via second
#' differences
#' \deqn{
#'   d_j = \tau_j - 2\tau_{j+1} + \tau_{j+2},
#' }
#' and their Z-tests. A global Wald test of equidistance is also
#' computed using all second differences jointly (when there are
#' at least four response categories).
#'
#' The z-tests for threshold differences are computed as
#' \deqn{ z = \frac{ \tau_{i+1} - \tau_i }{ \sqrt{SE(\tau_{i+1})^2 + SE(\tau_i)^2} } }
#' where \eqn{\tau_i} is the i-th threshold and \eqn{SE} its standard error.
#'
#' Thresholds are estimated using the \pkg{ordinal} package, which
#' implements cumulative link models via maximum likelihood.
#'
#' @param x An ordered factor representing a Likert-type item.
#' @param link The link function to use in the CLM. Options include
#' "logit", "probit" (default), "cloglog", etc.
#' @param conf.level Confidence level for the threshold intervals.
#' Default is 0.95.
#' @param nd Number of digits used to round numeric outputs.
#' Default is 3.
#'
#' @return A list with five elements:
#' \item{thresholds}{Threshold estimates and their standard errors.}
#' \item{ci}{Confidence intervals for each threshold.}
#' \item{differences}{Differences between consecutive thresholds
#' (interpretable as standardized effect sizes on the latent scale).}
#' \item{diff_tests}{z-tests for the differences between consecutive
#' thresholds, including p-values.}
#' \item{equidistance}{A list with:
#'   \itemize{
#'     \item \code{second_diffs}: data frame with second differences
#'       \eqn{d_j}, their standard errors, Z-tests and p-values.
#'     \item \code{wald}: data frame with the global Wald statistic,
#'       degrees of freedom and p-value for equidistance (or NA if
#'       not applicable).
#'     \item \code{EI.spratto}: the Equidistance Index as defined
#'       above.
#'   }}
#'
#' @examples
#' set.seed(123)
#' y <- ordered(sample(1:5, 200, replace = TRUE,
#'                    prob = c(.1, .2, .3, .25, .15)))
#' res <- itemThresholds(y)
#' res$thresholds
#' res$ci
#' res$differences
#' res$diff_tests
#' res$equidistance
#'
#' @seealso \code{\link[ordinal]{clm}}
#'
#' @references
#' Christensen, R. H. B. (2019). *ordinal: Regression Models for
#' Ordinal Data*. R package version 2019.12-10.
#' https://CRAN.R-project.org/package=ordinal
#'
#' Spratto, E. M. (2018). \emph{In search of equality: Developing an equal interval Likert response scale
#' (Doctoral dissertations)}. https://commons.lib.jmu.edu/diss201019/172
#'
#' Sideridis, G., Tsaousis, I., & Ghamdi, H. (2022). Equidistant Response Options on Likert-Type
#' Instruments: Testing the Interval Scaling Assumption Using Mplus. \emph{Educational and Psychological Measurement},
#' 83(5), 885-906. https://doi.org/10.1177/00131644221130482
#'
#' @importFrom stats coef pnorm qnorm vcov pchisq
#' @importFrom ordinal clm
#'
#' @export
itemThresholds <- function(x, link = "probit", conf.level = 0.95, nd = 3) {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("Package 'ordinal' is required. Please install it.")
  }

  # Coerción automática a ordered
  if (!is.ordered(x)) {
    if (is.factor(x)) {
      x <- ordered(x, levels = levels(x))
    } else if (is.numeric(x) || is.integer(x)) {
      x <- ordered(x)
    } else {
      stop("The argument 'x' must be convertible to an ordered factor (factor, numeric, or integer).")
    }
  }


  # helper to round data.frames
  round_df <- function(df, digits) {
    if (nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  # Fit CLM without predictors
  fit <- ordinal::clm(x ~ 1, link = link)

  # Extract thresholds and SE (from summary)
  coefs <- coef(summary(fit))
  thr   <- coefs[grepl("\\|", rownames(coefs)), , drop = FALSE]
  est   <- thr[, "Estimate"]
  se    <- thr[, "Std. Error"]
  thr_names <- rownames(thr)
  m <- length(est)          # number of thresholds (K - 1)

  thresholds <- data.frame(
    threshold = thr_names,
    estimate  = est,
    se        = se,
    stringsAsFactors = FALSE
  )

  # Confidence intervals
  zval <- qnorm(1 - (1 - conf.level) / 2)
  ci <- data.frame(
    threshold = thr_names,
    lwr.ci = est - zval * se,
    upr.ci = est + zval * se,
    stringsAsFactors = FALSE
  )

  # Differences and z-tests between adjacent thresholds
  if (length(est) > 1) {
    diffs <- diff(est)
    diffs_df <- data.frame(
      contrast   = paste0(thr_names[-1], " - ", thr_names[-length(thr_names)]),
      difference = diffs,
      stringsAsFactors = FALSE
    )

    se_diffs <- sqrt(se[-1]^2 + se[-length(se)]^2)
    z_vals   <- diffs / se_diffs
    p_vals   <- 2 * (1 - pnorm(abs(z_vals)))

    diff_tests <- data.frame(
      contrast   = paste0(thr_names[-1], " - ", thr_names[-length(thr_names)]),
      difference = diffs,
      se.diff    = se_diffs,
      z          = z_vals,
      p.value    = p_vals,
      stringsAsFactors = FALSE
    )
  } else {
    diffs_df   <- data.frame()
    diff_tests <- data.frame()
  }

  ## ---------- Equidistance: Spratto index, Z-test(s) and Wald ----------
  # VCOV of thresholds
  V_all <- vcov(fit)
  par_names <- rownames(V_all)
  idx_thr <- grepl("\\|", par_names)
  V <- as.matrix(V_all[idx_thr, idx_thr, drop = FALSE])

  # Equidistance Index (Spratto) based on first differences
  if (length(est) > 1) {
    deltas <- diff(est)
    mean_delta <- mean(deltas)
    EI_spratto <- sqrt(mean((deltas - mean_delta)^2))
  } else {
    EI_spratto <- NA_real_
  }

  # Second differences (only if at least 3 thresholds -> K >= 4)
  if (m >= 3) {
    L <- m - 2L
    C <- matrix(0, nrow = L, ncol = m)
    for (j in 1:L) {
      C[j, j]     <- 1
      C[j, j + 1] <- -2
      C[j, j + 2] <- 1
    }

    sec_vals <- as.numeric(C %*% est)
    # names like "t1 - 2*t2 + t3"
    sec_names <- vapply(
      1:L,
      function(j) {
        paste0(thr_names[j], " - 2*", thr_names[j + 1], " + ", thr_names[j + 2])
      },
      FUN.VALUE = character(1)
    )

    Var_mat <- C %*% V %*% t(C)
    se_sec  <- sqrt(diag(Var_mat))
    Z_sec   <- sec_vals / se_sec
    p_sec   <- 2 * (1 - pnorm(abs(Z_sec)))

    second_diffs <- data.frame(
      contrast = sec_names,
      value    = sec_vals,
      SE       = se_sec,
      Z        = Z_sec,
      p.value  = p_sec,
      stringsAsFactors = FALSE
    )

    # Wald global
    wald_stat <- as.numeric(t(sec_vals) %*% solve(Var_mat, sec_vals))
    df_wald   <- L
    p_wald    <- pchisq(wald_stat, df = df_wald, lower.tail = FALSE)

    wald_df <- data.frame(
      stat    = wald_stat,
      df      = df_wald,
      p.value = p_wald,
      stringsAsFactors = FALSE
    )
  } else {
    second_diffs <- data.frame()
    wald_df <- data.frame(
      stat    = NA_real_,
      df      = 0L,
      p.value = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  equidistance <- list(
    second_diffs = round_df(second_diffs, nd),
    wald         = round_df(wald_df, nd),
    EI.spratto   = round(EI_spratto, nd)
  )

  # Apply rounding to main data.frames
  thresholds  <- round_df(thresholds, nd)
  ci          <- round_df(ci, nd)
  diffs_df    <- round_df(diffs_df, nd)
  diff_tests  <- round_df(diff_tests, nd)

  return(list(
    thresholds   = thresholds,
    ci           = ci,
    differences  = diffs_df,
    diff_tests   = diff_tests,
    equidistance = equidistance
  ))
}
