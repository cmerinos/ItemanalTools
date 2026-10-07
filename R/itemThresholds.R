#' @title Item Thresholds, Monotonicity and Equidistance (Single Item)
#'
#' @description
#' Estimates thresholds (cutpoints) for an ordinal Likert-type item using a
#' cumulative link model (CLM). Provides estimates, confidence intervals,
#' differences between consecutive thresholds, z-tests for those differences
#' (with standard errors that account for the covariance between thresholds),
#' a monotonicity diagnostic, standardized effect sizes (Cohen's d), and
#' indices/tests of threshold equidistance (Spratto-type index).
#'
#' @details
#' This function fits an ordinal regression model without predictors
#' (\code{ordinal::clm}) and extracts threshold estimates. Three properties
#' of the response categories are diagnosed: **distance**, **monotonicity**,
#' and **intervality (equidistance)**.
#'
#' ## Distance between consecutive thresholds
#'
#' Let \eqn{\Delta_k = \tau_{k+1} - \tau_k} be the difference between
#' consecutive thresholds. Its standard error is computed using the full
#' covariance matrix of the thresholds,
#' \deqn{
#'   SE(\Delta_k) = \sqrt{ SE_k^2 + SE_{k+1}^2 - 2\,Cov(\tau_k, \tau_{k+1}) },
#' }
#' which accounts for the correlation between threshold estimates. Ignoring
#' this covariance (as in \eqn{\sqrt{SE_k^2 + SE_{k+1}^2}}) produces biased
#' z-tests and p-values.
#'
#' ## Monotonicity diagnostic
#'
#' For a cumulative link model to be appropriate, thresholds must be
#' strictly increasing (\eqn{\tau_1 < \tau_2 < \dots < \tau_{K-1}}). For each
#' adjacent pair, the function reports:
#' \itemize{
#'   \item the raw difference \eqn{\Delta_k};
#'   \item a z statistic \eqn{z_k = \Delta_k / SE(\Delta_k)};
#'   \item a one-tailed p-value testing \eqn{H_0: \Delta_k \ge 0} versus
#'         \eqn{H_1: \Delta_k < 0} (i.e., disordered thresholds);
#'   \item a binary flag \code{disordered}, equal to 1 if the raw difference
#'         is negative and 0 otherwise.
#' }
#' The flag reflects the raw sign; users can apply the formal test using
#' \code{p.one.sided} with their preferred significance level. A convenient
#' one-tailed criterion at \eqn{\alpha = .05} is \eqn{z_k < -1.645}.
#'
#' ## Standardized effect size (Cohen's d)
#'
#' Raw differences are expressed as standardized effect sizes. With
#' \code{link = "probit"} (default), the latent scale has a standard
#' deviation of 1, so the raw differences are **already interpretable as
#' Cohen's d**. With \code{link = "logit"}, the latent scale has a standard
#' deviation of \eqn{\pi/\sqrt{3} \approx 1.8138}, so raw differences are
#' divided by this constant to obtain d. For other link functions, d is set
#' to \code{NA} and only raw differences are reported.
#'
#' Benchmarks (Cohen, 1988; adapted for threshold distances):
#' \itemize{
#'   \item \eqn{d < 0.20}: categories are functionally redundant;
#'   \item \eqn{0.20 \le d < 0.50}: weak distinction;
#'   \item \eqn{d \ge 0.50}: clear distinction.
#' }
#' The column \code{redundant} flags pairs with \eqn{d < 0.20}.
#'
#' ## Equidistance (intervality)
#'
#' Let \eqn{\Delta_k = \tau_{k+1} - \tau_k}. The Equidistance Index
#' (Spratto) is
#' \deqn{
#'   EI = \sqrt{\frac{1}{K-2}\sum_{k=1}^{K-2}(\Delta_k - \bar{\Delta})^2},
#' }
#' where \eqn{K} is the number of response categories and \eqn{\bar{\Delta}}
#' is the mean of the \eqn{\Delta_k}. Local deviations from equidistance are
#' assessed via second differences
#' \eqn{d_j = \tau_j - 2\tau_{j+1} + \tau_{j+2}} and their Z-tests. A global
#' Wald test of equidistance is computed when there are at least four
#' response categories.
#'
#' ## Scope
#'
#' This function diagnoses the **categorization of responses** at the item
#' level. It does **not** assume a latent dimension underlying multiple
#' items, and it does not estimate item discrimination. For latent-model
#' diagnostics of thresholds, see \code{\link{itemThresholdsLat}}.
#'
#' @param x An ordered factor representing a Likert-type item.
#' @param link The link function to use in the CLM. Options include
#'   \code{"logit"}, \code{"probit"} (default), \code{"cloglog"}, etc.
#' @param conf.level Confidence level for the threshold intervals.
#'   Default is 0.95.
#' @param nd Number of digits used to round numeric outputs. Default is 3.
#'
#' @return A list with the following elements:
#' \item{thresholds}{Threshold estimates and their standard errors.}
#' \item{ci}{Confidence intervals for each threshold.}
#' \item{differences}{Differences between consecutive thresholds
#'   (interpretable as standardized effect sizes on the latent scale).}
#' \item{diff_tests}{z-tests for the differences between consecutive
#'   thresholds, including two-tailed p-values, one-tailed p-values for the
#'   monotonicity test, Cohen's d, and flags for redundant and disordered
#'   pairs.}
#' \item{monotonicity}{A focused table with the monotonicity diagnostic:
#'   contrast, raw difference, standard error, z statistic, one-tailed
#'   p-value, and \code{disordered} flag.}
#' \item{any_disord}{Integer flag (0/1) indicating whether any adjacent pair
#'   of thresholds is disordered at the item level.}
#' \item{equidistance}{A list with:
#'   \itemize{
#'     \item \code{second_diffs}: second differences, their standard errors,
#'       Z-tests, and p-values;
#'     \item \code{wald}: global Wald statistic, degrees of freedom, and
#'       p-value for equidistance;
#'     \item \code{EI.spratto}: the Equidistance Index.
#'   }}
#'
#' @examples
#' set.seed(123)
#' y <- ordered(sample(1:5, 200, replace = TRUE,
#'                    prob = c(.1, .2, .3, .25, .15)))
#' res <- itemThresholds(y)
#' res$thresholds
#' res$diff_tests
#' res$monotonicity
#' res$any_disord
#' res$equidistance
#'
#' @seealso \code{\link[ordinal]{clm}}
#'
#' @references
#' Christensen, R. H. B. (2019). \emph{ordinal: Regression Models for
#' Ordinal Data}. R package version 2019.12-10.
#' https://CRAN.R-project.org/package=ordinal
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral
#' Sciences} (2nd ed.). Lawrence Erlbaum Associates.
#'
#' Spratto, E. M. (2018). \emph{In search of equality: Developing an equal
#' interval Likert response scale} (Doctoral dissertation).
#' https://commons.lib.jmu.edu/diss201019/172
#'
#' Sideridis, G., Tsaousis, I., & Ghamdi, H. (2022). Equidistant Response
#' Options on Likert-Type Instruments: Testing the Interval Scaling
#' Assumption Using Mplus. \emph{Educational and Psychological Measurement},
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

  # Coercion to ordered factor
  if (!is.ordered(x)) {
    if (is.factor(x)) {
      x <- ordered(x, levels = levels(x))
    } else if (is.numeric(x) || is.integer(x)) {
      x <- ordered(x)
    } else {
      stop("The argument 'x' must be convertible to an ordered factor ",
           "(factor, numeric, or integer).")
    }
  }

  # Helper to round numeric columns of a data.frame
  round_df <- function(df, digits) {
    if (nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  # Fit CLM without predictors
  fit <- ordinal::clm(x ~ 1, link = link)

  # Extract thresholds and SEs
  coefs <- coef(summary(fit))
  thr   <- coefs[grepl("\\|", rownames(coefs)), , drop = FALSE]
  est   <- thr[, "Estimate"]
  se    <- thr[, "Std. Error"]
  thr_names <- rownames(thr)
  m <- length(est)  # number of thresholds (K - 1)

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

  # Full covariance matrix of thresholds
  V_all <- vcov(fit)
  par_names <- rownames(V_all)
  idx_thr <- grepl("\\|", par_names)
  V <- as.matrix(V_all[idx_thr, idx_thr, drop = FALSE])

  # Latent SD depending on link (for Cohen's d)
  link_sd <- switch(
    tolower(link),
    "probit" = 1,
    "logit"  = pi / sqrt(3),
    NA_real_
  )

  # ---------------------------------------------------------------
  # Differences between consecutive thresholds + monotonicity
  # ---------------------------------------------------------------
  if (m > 1) {
    # Contrast matrix for adjacent differences: Delta_k = tau_{k+1} - tau_k
    C_diff <- matrix(0, nrow = m - 1, ncol = m)
    for (k in seq_len(m - 1)) {
      C_diff[k, k]     <- -1
      C_diff[k, k + 1] <-  1
    }
    diffs    <- as.numeric(C_diff %*% est)
    se_diffs <- sqrt(diag(C_diff %*% V %*% t(C_diff)))
    z_vals   <- diffs / se_diffs
    p_2sided <- 2 * (1 - pnorm(abs(z_vals)))
    p_1sided <- pnorm(z_vals)  # H1: Delta < 0 (disorder)

    contrast_names <- paste0(thr_names[-1], " - ", thr_names[-length(thr_names)])

    diffs_df <- data.frame(
      contrast   = contrast_names,
      difference = diffs,
      stringsAsFactors = FALSE
    )

    # Cohen's d and flags
    d_cohen    <- if (!is.na(link_sd)) diffs / link_sd else rep(NA_real_, length(diffs))
    redundant  <- as.integer(!is.na(d_cohen) & d_cohen < 0.20)
    disordered <- as.integer(diffs < 0)

    diff_tests <- data.frame(
      contrast    = contrast_names,
      difference  = diffs,
      se.diff     = se_diffs,
      z           = z_vals,
      p.value     = p_2sided,
      p.one.sided = p_1sided,
      d.cohen     = d_cohen,
      redundant   = redundant,
      disordered  = disordered,
      stringsAsFactors = FALSE
    )

    monotonicity <- data.frame(
      contrast    = contrast_names,
      difference  = diffs,
      se.diff     = se_diffs,
      z           = z_vals,
      p.one.sided = p_1sided,
      disordered  = disordered,
      stringsAsFactors = FALSE
    )
    any_disord <- as.integer(any(disordered == 1L))
  } else {
    diffs_df     <- data.frame()
    diff_tests   <- data.frame()
    monotonicity <- data.frame()
    any_disord   <- 0L
  }

  # ---------------------------------------------------------------
  # Equidistance: Spratto index, second differences, Wald test
  # ---------------------------------------------------------------
  if (m > 1) {
    deltas <- diff(est)
    mean_delta <- mean(deltas)
    EI_spratto <- sqrt(mean((deltas - mean_delta)^2))
  } else {
    EI_spratto <- NA_real_
  }

  if (m >= 3) {
    L <- m - 2L
    C <- matrix(0, nrow = L, ncol = m)
    for (j in 1:L) {
      C[j, j]     <- 1
      C[j, j + 1] <- -2
      C[j, j + 2] <- 1
    }
    sec_vals <- as.numeric(C %*% est)
    sec_names <- vapply(
      1:L,
      function(j) paste0(thr_names[j], " - 2*", thr_names[j + 1], " + ", thr_names[j + 2]),
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

  # Round main outputs
  thresholds   <- round_df(thresholds, nd)
  ci           <- round_df(ci, nd)
  diffs_df     <- round_df(diffs_df, nd)
  diff_tests   <- round_df(diff_tests, nd)
  monotonicity <- round_df(monotonicity, nd)

  return(list(
    thresholds   = thresholds,
    ci           = ci,
    differences  = diffs_df,
    diff_tests   = diff_tests,
    monotonicity = monotonicity,
    any_disord   = any_disord,
    equidistance = equidistance
  ))
}
