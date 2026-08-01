#' Ulrich-Wirtz Correlation between a Naturally and an Artificially Dichotomized Variable
#'
#' Calculates the correlation coefficient \eqn{\nu} proposed by Ulrich and Wirtz (2004)
#' for assessing the relationship between a naturally dichotomous variable (\eqn{X})
#' and an artificially dichotomized variable (\eqn{Y}) derived from an underlying
#' continuous latent variable \eqn{L}.
#'
#' @param data.items A data frame or matrix with dichotomous items (numeric 0/1)
#'   in columns and subjects in rows. These are treated as the naturally
#'   dichotomous variables \eqn{X}.
#' @param criterion A dichotomous vector (factor, logical, or numeric with two
#'   unique values). This is treated as the artificially dichotomized variable
#'   \eqn{Y}. The user may specify which level is coded as 0 via `lower.level`.
#' @param lower.level Specifies which level of the criterion is considered the
#'   "lower" group (coded as 0). For factors, the default is the first level;
#'   for numeric vectors, the minimum value; for logical, `FALSE`.
#' @param add A numeric constant added to all four cell counts of the 2x2 table
#'   before computing proportions. This is a continuity correction for tables
#'   with zero cells. Common values are `0.5` (Haldane correction) or smaller
#'   constants like `0.1`. Default is `NULL` (no correction).
#' @param ci Logical. If `TRUE`, bootstrap confidence intervals and a two-sided
#'   p-value are computed. Default `FALSE`.
#' @param conf.level Confidence level for the bootstrap intervals. Default 0.95.
#'   Used only if `ci = TRUE`.
#' @param B Number of bootstrap replicates (if `ci = TRUE`). Default 1000.
#' @param min.boot Minimum number of valid (non-NA) bootstrap replicates required to calculate confidence intervals and the p-value. If there are fewer, NA is returned along with a warning.
#' @param type Method for bootstrap confidence intervals: `"perc"`, `"norm"`, or
#'   `"bca"`. Default `"perc"`. Ignored if `ci = FALSE`.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Item}{Name of the item (column name from `data.items`).}
#'     \item{r.nu}{Estimated correlation coefficient \eqn{\nu}.}
#'   }
#'   If `ci = TRUE`, additionally:
#'   \describe{
#'     \item{lwr.ci}{Lower bound of the bootstrap confidence interval.}
#'     \item{upr.ci}{Upper bound of the bootstrap confidence interval.}
#'     \item{p.value}{Two-sided bootstrap p-value (proportion of bootstrap estimates
#'       with absolute value >= observed |\eqn{\nu}|).}
#'   }
#'
#' @details
#' The coefficient \eqn{\nu} (Equation 1 in Ulrich & Wirtz, 2004) recovers the
#' product-moment correlation between a naturally dichotomous variable \eqn{X} and
#' a continuous latent variable \eqn{L} underlying the artificial dichotomy \eqn{Y}.
#' Under the assumptions that:
#' \enumerate{
#'   \item The conditional distribution of \eqn{L} given \eqn{X = x} is normal.
#'   \item The conditional variance of \eqn{L} does not depend on \eqn{X} (homoscedasticity).
#' }
#' the coefficient is computed from the 2x2 table of \eqn{X} and \eqn{Y} as:
#' \deqn{
#'   \nu = \frac{\Delta}{\sqrt{\Delta^2 + \frac{1}{p_1(1-p_1)}}}, \quad
#'   \Delta = \Phi^{-1}(p_{00}/p_0) - \Phi^{-1}(p_{10}/p_1)
#' }
#' where \eqn{p_0 = \Pr(X=0)}, \eqn{p_1 = \Pr(X=1)},
#' \eqn{p_{00} = \Pr(X=0, Y=0)}, \eqn{p_{10} = \Pr(X=1, Y=0)}, and
#' \eqn{\Phi^{-1}} is the quantile function of the standard normal distribution.
#'
#' When a cell in the 2x2 table is zero, the coefficient \eqn{\nu} is undefined
#' (infinite quantiles). To avoid losing such items, the user may supply a small
#' constant via the argument `add`, which is added to all four cell counts before
#' computing proportions. This is a continuity correction analogous to Haldane's
#' correction for the log-odds ratio (Haldane, 1956). Brown and Benedetti (1977)
#' showed that adding 0.5 to each cell reduces bias in tetrachoric correlations
#' when minimum cell frequencies are less than 5. Values between 0.01 and 0.5 are
#' commonly used; 0.5 is the standard, while smaller values introduce less bias
#' but may still yield infinite estimates if the corrected cell remains zero.
#' The default `add = NULL` performs no correction and returns `NA` for undefined
#' items.
#'
#' Ulrich and Wirtz also proposed a logistic approximation (Equation 3) that
#' replaces the normal quantiles by the log-odds ratio, yielding virtually
#' identical results (\eqn{|\Psi(x) - \Phi(x)| < 0.01}). This approximation is
#' not implemented here because the exact normal-quantile formulation is
#' computationally trivial in R and is the theoretical foundation of the
#' coefficient. The logistic version is mentioned for completeness and historical
#' context.
#'
#' @references
#' Ulrich, R., & Wirtz, M. (2004). On the correlation of a naturally and an
#' artificially dichotomized variable. \emph{British Journal of Mathematical and
#' Statistical Psychology, 57}(2), 235-251. \doi{10.1348/0007110042307203}
#'
#' Brown, M. B., & Benedetti, J. K. (1977). On the mean and variance of the
#' tetrachoric correlation coefficient. \emph{Psychometrika, 42}(3), 347-355. \doi{10.1007/BF02293655}
#'
#' Haldane, J. B. S. (1956). The estimation and significance of the logarithm of
#' a ratio of frequencies. \emph{Annals of Human Genetics, 20}(4), 309-311. \doi{j.1469-1809.1955.tb01285.x}
#'
#' Islam, T. U., & Rizwan, M. (2020). Comparison of Correlation Measures for
#' Nominal Data. \emph{Preprints}, \doi{10.20944/preprints202004.0276.v1}.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' items <- data.frame(
#'   Item1 = sample(0:1, 100, replace = TRUE),
#'   Item2 = sample(0:1, 100, replace = TRUE)
#' )
#' criterion <- sample(0:1, 100, replace = TRUE)  # artificial dichotomy
#'
#' # Without bootstrap
#' rvUlrich(items, criterion, ci = FALSE)
#'
#' # With bootstrap (percentile CIs)
#' rvUlrich(items, criterion, ci = TRUE, B = 500, type = "perc")
#'
#' # With Haldane correction for zero cells
#' rvUlrich(items, criterion, add = 0.5, ci = FALSE)
#' }
#'
#' @importFrom stats qnorm  complete.cases quantile sd
#'
#' @export
rvUlrich <- function(data.items,
                     criterion,
                     lower.level = NULL,
                     add = NULL,
                     ci = FALSE,
                     B = 1000,
                     type = "perc",
                     conf.level = 0.95,
                     min.boot = 50) {

  # --- 1. Validations ---
  if (!is.data.frame(data.items) && !is.matrix(data.items)) {
    stop("'data.items' must be a data frame or matrix.")
  }
  data.items <- as.data.frame(data.items)
  if (!all(sapply(data.items, is.numeric))) {
    stop("All columns in 'data.items' must be numeric (preferably 0/1).")
  }
  if (length(criterion) != nrow(data.items)) {
    stop("'criterion' must have the same length as the number of rows in 'data.items'.")
  }
  if (!is.null(add) && (!is.numeric(add) || length(add) != 1 || add < 0)) {
    stop("'add' must be a single non-negative numeric value or NULL.")
  }
  if (ci && !type %in% c("perc", "norm")) {
    stop("Invalid 'type'. Use 'perc' or 'norm'.")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single numeric value between 0 and 1.")
  }

  # --- 2. Process criterion into binary 0/1 vector ---
  if (is.factor(criterion)) {
    if (nlevels(criterion) != 2) {
      stop("'criterion' must be a factor with exactly two levels.")
    }
    if (is.null(lower.level)) lower.level <- levels(criterion)[1]
    dic <- ifelse(criterion == lower.level, 0, 1)
  } else if (is.logical(criterion)) {
    if (is.null(lower.level)) lower.level <- FALSE
    dic <- ifelse(criterion == lower.level, 0, 1)
  } else if (is.numeric(criterion)) {
    if (length(unique(criterion)) != 2) {
      stop("Numeric 'criterion' must have exactly two unique values.")
    }
    if (is.null(lower.level)) lower.level <- min(criterion, na.rm = TRUE)
    dic <- ifelse(criterion == lower.level, 0, 1)
  } else {
    stop("'criterion' must be a factor, logical, or numeric vector with two values.")
  }

  if (!any(dic == 0) || !any(dic == 1)) {
    stop("The specified 'lower.level' does not produce both groups.")
  }

  # Remove rows with NAs in items or criterion
  complete <- complete.cases(data.items, criterion)
  if (!all(complete)) {
    warning(sprintf("Removed %d rows with missing values.", sum(!complete)))
    data.items <- data.items[complete, , drop = FALSE]
    dic <- dic[complete]
  }

  N <- nrow(data.items)
  if (N < 3) stop("At least 3 complete observations are required.")

  # --- 3. Internal function to compute ν for one item ---
  calc_nu <- function(x) {
    tab <- table(x, dic)
    if (any(dim(tab) < 2)) return(NA_real_)
    n00 <- tab[1, 1]; n01 <- tab[1, 2]
    n10 <- tab[2, 1]; n11 <- tab[2, 2]

    if (!is.null(add)) {
      n00 <- n00 + add; n01 <- n01 + add
      n10 <- n10 + add; n11 <- n11 + add
      N_corr <- n00 + n01 + n10 + n11
    } else {
      N_corr <- N
    }

    p0 <- (n00 + n01) / N_corr
    p1 <- (n10 + n11) / N_corr
    if (p0 == 0 || p1 == 0 || p0 == 1 || p1 == 1) return(NA_real_)

    p00 <- n00 / N_corr
    p10 <- n10 / N_corr
    if (p00 == 0 || p00 == p0 || p10 == 0 || p10 == p1) return(NA_real_)

    delta <- qnorm(p00 / p0) - qnorm(p10 / p1)
    nu <- delta / sqrt(delta^2 + 1 / (p1 * (1 - p1)))
    return(nu)
  }

  # --- 4. Compute ν for all items, handling constant items ---
  n_items <- ncol(data.items)
  nu_vals <- numeric(n_items)
  constant_items <- character(0)

  for (i in seq_len(n_items)) {
    x <- data.items[[i]]
    if (length(unique(x)) < 2) {
      constant_items <- c(constant_items, names(data.items)[i])
      nu_vals[i] <- NA
    } else {
      nu_vals[i] <- calc_nu(x)
    }
  }

  if (length(constant_items) > 0) {
    warning(sprintf(
      "The following items have no variability (all values equal) and were set to NA: %s",
      paste(constant_items, collapse = ", ")
    ))
  }

  df_out <- data.frame(Item = names(data.items),
                       r.nu = nu_vals,
                       stringsAsFactors = FALSE)

  # --- 5. Return if no bootstrap ---
  if (isFALSE(ci)) {
    num_cols <- sapply(df_out, is.numeric)
    df_out[num_cols] <- lapply(df_out[num_cols], round, 3)
    return(df_out)
  }

  # --- 6. Bootstrap CIs and p-value ---
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for bootstrap CIs. Please install it.")
  }

  # Bootstrap function: compute ν for all items on resampled data
  boot_nu <- function(data, indices, dic_orig, add_orig) {
    d <- data[indices, , drop = FALSE]
    dic_b <- dic_orig[indices]
    sapply(d, function(col) {
      # Check for constant item in bootstrap sample
      if (length(unique(col)) < 2) return(NA_real_)
      tab <- table(col, dic_b)
      if (any(dim(tab) < 2)) return(NA_real_)
      n00 <- tab[1, 1]; n01 <- tab[1, 2]
      n10 <- tab[2, 1]; n11 <- tab[2, 2]

      if (!is.null(add_orig)) {
        n00 <- n00 + add_orig; n01 <- n01 + add_orig
        n10 <- n10 + add_orig; n11 <- n11 + add_orig
        N_b <- n00 + n01 + n10 + n11
      } else {
        N_b <- length(dic_b)
      }

      p0 <- (n00 + n01) / N_b
      p1 <- (n10 + n11) / N_b
      if (p0 == 0 || p1 == 0 || p0 == 1 || p1 == 1) return(NA_real_)
      p00 <- n00 / N_b
      p10 <- n10 / N_b
      if (p00 == 0 || p00 == p0 || p10 == 0 || p10 == p1) return(NA_real_)
      delta <- qnorm(p00 / p0) - qnorm(p10 / p1)
      nu <- delta / sqrt(delta^2 + 1 / (p1 * (1 - p1)))
      return(nu)
    })
  }

  boot_obj <- boot::boot(data = data.items,
                         statistic = boot_nu,
                         R = B,
                         dic_orig = dic,
                         add_orig = add)

  # Extract bootstrap replicates matrix
  boot_coefs <- boot_obj$t

  # Compute percentiles based on conf.level
  alpha <- 1 - conf.level
  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2

  # For each item, compute CI and p-value manually
  lwr_ci <- upr_ci <- p_val <- rep(NA_real_, n_items)

  for (i in seq_len(n_items)) {
    obs <- df_out$r.nu[i]
    if (is.na(obs)) {
      lwr_ci[i] <- upr_ci[i] <- p_val[i] <- NA
      next
    }

    # Extract valid bootstrap replicates (non-NA)
    t_vals <- boot_coefs[, i]
    t_vals_ok <- t_vals[!is.na(t_vals)]
    n_ok <- length(t_vals_ok)

    if (n_ok < min.boot) {
      warning(sprintf(
        "Item '%s': only %d valid bootstrap replicates (min.boot = %d). CI and p-value set to NA.",
        names(data.items)[i], n_ok, min.boot
      ))
      lwr_ci[i] <- upr_ci[i] <- p_val[i] <- NA
      next
    }

    # Compute CI based on type
    if (type == "perc") {
      ci_perc <- quantile(t_vals_ok, probs = c(lower_prob, upper_prob), na.rm = TRUE)
      lwr_ci[i] <- ci_perc[1]
      upr_ci[i] <- ci_perc[2]
    } else { # type == "norm"
      m <- mean(t_vals_ok, na.rm = TRUE)
      s <- sd(t_vals_ok, na.rm = TRUE)
      z <- qnorm(1 - alpha / 2)
      lwr_ci[i] <- m - z * s
      upr_ci[i] <- m + z * s
    }

    # Compute two-sided p-value based on centered bootstrap distribution
    centered <- t_vals_ok - obs
    p_val[i] <- mean(abs(centered) >= abs(obs), na.rm = TRUE)
  }

  df_out$lwr.ci <- lwr_ci
  df_out$upr.ci <- upr_ci
  df_out$p.value <- p_val

  # Round numeric columns
  num_cols <- sapply(df_out, is.numeric)
  df_out[num_cols] <- lapply(df_out[num_cols], round, 3)

  return(df_out)
}
