#' Rothery's Nonparametric Intraclass Correlation Coefficient
#'
#' Computes Rothery's rank-based intraclass correlation coefficient (rho) using
#' the exact concordance estimator (Rothery, 1979). P-values can be obtained via
#' bootstrap or using the revised-beta approximation (as in nopaco). Confidence
#' intervals are computed by bootstrap, with options for percentile, normal, or
#' BCa methods.
#'
#' @param data A numeric matrix or data frame with subjects in rows and items
#'   (repeated measures) in columns. Missing values are not allowed; rows with
#'   NAs are removed.
#' @param ci Logical; if TRUE, bootstrap confidence intervals are computed.
#' @param B Integer; number of bootstrap replicates (ignored if ci = FALSE).
#' @param conf.level Confidence level for the intervals (default 0.95).
#' @param type Character; type of bootstrap confidence interval: \code{"perc"}
#'   (percentile), \code{"norm"} (normal approximation), or \code{"bca"}
#'   (bias-corrected and accelerated). Default is \code{"perc"}.
#' @param p.method Character; method for p-value: \code{"beta"} (revised-beta
#'   approximation, default) or \code{"bootstrap"} (based on the bootstrap
#'   distribution). The beta approximation is faster and generally more accurate.
#' @param alternative Character; alternative hypothesis: \code{"greater"}
#'   (default, tests if concordance is greater than chance) or \code{"two.sided"}.
#' @param seed Optional integer to set the random seed for reproducibility.
#'
#' @return A data frame with columns:
#'   \item{rho}{Point estimate of Rothery's rho (scaled: from -1 to 1).}
#'   \item{r_c}{Raw concordance estimate (on the native scale).}
#'   \item{ci.lower}{Lower confidence limit (NA if ci = FALSE).}
#'   \item{ci.upper}{Upper confidence limit (NA if ci = FALSE).}
#'   \item{p.value}{P-value (computed via the chosen method).}
#'   \item{method}{Character string indicating the p-value method used.}
#'   \item{alternative}{Character string indicating the alternative hypothesis.}
#'
#' @details
#' \subsection{Estimation}{
#'   For each subject (row), the sum of the within-subject absolute differences
#'   of ranks (mid‑ranks) is computed. The estimator \eqn{r_c} is the proportion
#'   of concordant triplets among all possible triplets formed by two
#'   observations from the same subject and one from a different subject.
#'   The formula used is exactly equation (2.1) in Rothery (1979).
#'
#'   The raw estimate \eqn{r_c} is then transformed to \eqn{\rho} (rho) using
#'   \deqn{\rho = 2 \sin(\pi (r_c - 1/2)) - 1,}
#'   which maps the range of \eqn{r_c} (typically between ~0.5 and 1) to
#'   \eqn{[-1, 1]}. Under the null hypothesis of no concordance (balanced
#'   design), \eqn{r_c = 2/3} and \eqn{\rho = 0}.
#' }
#'
#' \subsection{Bootstrap}{
#'   Bootstrap samples are obtained by resampling subjects (rows) with
#'   replacement, keeping all columns. For each bootstrap sample, \eqn{r_c} and
#'   \eqn{\rho} are recomputed. Confidence intervals are obtained from the
#'   bootstrap distribution of \eqn{\rho}. Three types are available:
#'   \itemize{
#'     \item \code{"perc"}: simple percentile interval.
#'     \item \code{"norm"}: interval based on the normal approximation with
#'           bootstrap standard error.
#'     \item \code{"bca"}: bias-corrected and accelerated interval (Efron, 1987),
#'           which adjusts for skewness and bias.
#'   }
#'
#'   The p-value can be computed either via the bootstrap distribution
#'   (\code{p.method = "bootstrap"}) or using the revised-beta approximation
#'   (\code{p.method = "beta"}) as implemented in the \code{nopaco} package.
#'   The beta approximation is faster and recommended for most applications.
#' }
#'
#' \subsection{Handling of ties}{
#'   Ranks are computed column-wise using \code{ties.method = "average"}
#'   (mid‑ranks), which is appropriate for ordinal data and follows Rothery's
#'   recommendation.
#' }
#'
#' \subsection{Note on balanced designs}{
#'   The function assumes that all subjects have the same number of items
#'   (columns) after removing rows with missing values. If your data are
#'   unbalanced (different numbers of items per subject), you need to adapt the
#'   calculation of \eqn{S_3} and \eqn{C_\alpha} using the per‑subject
#'   \eqn{n_\alpha} (the supplementary code shows how). This implementation
#'   focuses on the common balanced case. Here, Rocha et al. (2026) R code was adaptated.
#' }
#'
#' @references
#' Rothery, P. (1979). A nonparametric measure of intraclass correlation.
#' \emph{Biometrika}, 66(3), 629-639. \doi{10.1093/biomet/66.3.629}
#'
#' Rocha, B., Jonaitis, E. M., Hamwi, A., & Engelman, C. D. (2026). Longitudinal,
#' Intra-Individual Stability of Untargeted Plasma and Cerebrospinal Fluid
#' Metabolites. \emph{Metabolites}, 16(1), 35.
#' \doi{10.3390/metabo16010035}
#'
#' @examples
#' # Simulate balanced data: 30 subjects, 5 items with intra-subject correlation
#' set.seed(123)
#' dat <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
#' for (i in 1:30) dat[i, ] <- dat[i, ] + rnorm(1, sd = 1.5)
#' RotheryICC(dat, ci = TRUE, B = 500, type = "bca")
#'
#' @importFrom stats pbeta quantile sd
#'
#' @export
RotheryICC <- function(data,
                       ci = TRUE,
                       B = 1000,
                       conf.level = 0.95,
                       type = "perc",
                       p.method = "beta",
                       alternative = "greater",
                       seed = NULL) {

  # --- Helper functions (from supplemental_code_rotherys_icc.R) ---

  C_alpha <- function(N, n_alpha, R) {
    const <- (1/2) * (N - n_alpha) * n_alpha * (n_alpha - 1) +
      (1/6) * n_alpha * (n_alpha^2 - 1)
    sum_abs <- sum(abs(outer(R, R, "-")))
    out <- const - 0.5 * sum_abs
    return(out)
  }

  S3_alpha <- function(N, n_alpha) {
    0.5 * n_alpha * (n_alpha - 1) * (N - n_alpha)
  }

  # --- Main computation ---

  # Prepare data
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame.")
  }
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("All entries in 'data' must be numeric.")
  }
  na_rows <- apply(data, 1, function(x) any(is.na(x)))
  if (any(na_rows)) {
    warning(sprintf("Removed %d row(s) with missing values.", sum(na_rows)))
    data <- data[!na_rows, , drop = FALSE]
  }
  n_subj <- nrow(data)
  k_items <- ncol(data)
  if (n_subj < 2 || k_items < 2) {
    stop("Need at least 2 subjects and 2 items after removing NAs.")
  }

  # Compute r_c (exact formula)
  N <- n_subj * k_items
  R <- apply(data, 2, rank, ties.method = "average")

  C_sum <- 0
  S3_sum <- 0
  for (i in 1:n_subj) {
    row_ranks <- R[i, ]
    C_sum <- C_sum + C_alpha(N, k_items, row_ranks)
    S3_sum <- S3_sum + S3_alpha(N, k_items)
  }
  r_c <- C_sum / S3_sum
  rho <- 2 * sin(pi * (r_c - 0.5)) - 1

  # --- P-value ---

  if (p.method == "beta") {
    bn <- rep(k_items, n_subj)
    omega <- sum(bn * (bn - 1) * (sum(bn) - bn))
    v <- 4 * (N + 1) / (45 * n_subj^2 * k_items * (k_items - 1) * (n_subj - 1))
    meanB <- k_items
    zeta <- 2/3 - sqrt(meanB + 1) / (9/2 * (meanB - 1)^1.5)
    iii <- 2 / omega
    alphaPar <- (4*(meanB+1)/(81*(meanB-1)^3) -
                   ((8/9)*(meanB+1)^1.5)/(81*(meanB-1)^4.5)) / v -
      sqrt(4*(meanB+1)/(81*(meanB-1)^3))
    betaPar <- alphaPar * (9*(meanB-1)^1.5 / (2*sqrt(meanB+1)) - 1)
    p_val <- pbeta(r_c - zeta - iii,
                   shape1 = alphaPar,
                   shape2 = betaPar,
                   lower.tail = FALSE)
    if (alternative == "two.sided") {
      p_val <- 2 * min(p_val, 1 - p_val)
    }
    p_method_label <- "revised-beta"
  } else if (p.method == "bootstrap") {
    if (!is.null(seed)) set.seed(seed)
    boot_r_c <- numeric(B)
    for (b in 1:B) {
      idx <- sample(n_subj, replace = TRUE)
      boot_data <- data[idx, ]
      R_b <- apply(boot_data, 2, rank, ties.method = "average")
      C_sum_b <- 0
      S3_sum_b <- 0
      for (i in 1:n_subj) {
        row_ranks <- R_b[i, ]
        C_sum_b <- C_sum_b + C_alpha(N, k_items, row_ranks)
        S3_sum_b <- S3_sum_b + S3_alpha(N, k_items)
      }
      boot_r_c[b] <- C_sum_b / S3_sum_b
    }
    if (alternative == "greater") {
      p_val <- mean(boot_r_c >= r_c, na.rm = TRUE)
    } else if (alternative == "two.sided") {
      p_val <- 2 * min(mean(boot_r_c >= r_c, na.rm = TRUE),
                       mean(boot_r_c <= r_c, na.rm = TRUE))
    } else {
      stop("alternative must be 'greater' or 'two.sided'.")
    }
    p_val <- min(p_val, 1)
    p_method_label <- "bootstrap"
  } else {
    stop("p.method must be 'beta' or 'bootstrap'.")
  }

  # --- Confidence intervals (bootstrap) ---

  if (ci) {
    if (!is.null(seed)) set.seed(seed)
    boot_rho <- numeric(B)
    for (b in 1:B) {
      idx <- sample(n_subj, replace = TRUE)
      boot_data <- data[idx, ]
      R_b <- apply(boot_data, 2, rank, ties.method = "average")
      C_sum_b <- 0
      S3_sum_b <- 0
      for (i in 1:n_subj) {
        row_ranks <- R_b[i, ]
        C_sum_b <- C_sum_b + C_alpha(N, k_items, row_ranks)
        S3_sum_b <- S3_sum_b + S3_alpha(N, k_items)
      }
      r_c_b <- C_sum_b / S3_sum_b
      boot_rho[b] <- 2 * sin(pi * (r_c_b - 0.5)) - 1
    }

    alpha <- 1 - conf.level
    if (type == "perc") {
      ci_lwr <- quantile(boot_rho, probs = alpha / 2, na.rm = TRUE)
      ci_upr <- quantile(boot_rho, probs = 1 - alpha / 2, na.rm = TRUE)
    } else if (type == "norm") {
      se <- sd(boot_rho, na.rm = TRUE)
      ci_lwr <- mean(boot_rho, na.rm = TRUE) - qnorm(1 - alpha / 2) * se
      ci_upr <- mean(boot_rho, na.rm = TRUE) + qnorm(1 - alpha / 2) * se
    } else if (type == "bca") {
      z0 <- qnorm(sum(boot_rho < rho, na.rm = TRUE) / B)
      jack_rho <- numeric(n_subj)
      for (i in 1:n_subj) {
        jack_data <- data[-i, , drop = FALSE]
        n_j <- nrow(jack_data)
        N_j <- n_j * k_items
        R_j <- apply(jack_data, 2, rank, ties.method = "average")
        C_sum_j <- 0
        S3_sum_j <- 0
        for (l in 1:n_j) {
          row_ranks <- R_j[l, ]
          C_sum_j <- C_sum_j + C_alpha(N_j, k_items, row_ranks)
          S3_sum_j <- S3_sum_j + S3_alpha(N_j, k_items)
        }
        r_c_j <- C_sum_j / S3_sum_j
        jack_rho[i] <- 2 * sin(pi * (r_c_j - 0.5)) - 1
      }
      mean_jack <- mean(jack_rho, na.rm = TRUE)
      num <- sum((mean_jack - jack_rho)^3, na.rm = TRUE)
      den <- 6 * (sum((mean_jack - jack_rho)^2, na.rm = TRUE))^(3/2)
      a <- num / den
      alpha1 <- pnorm(z0 + (z0 + qnorm(alpha / 2)) / (1 - a * (z0 + qnorm(alpha / 2))))
      alpha2 <- pnorm(z0 + (z0 + qnorm(1 - alpha / 2)) / (1 - a * (z0 + qnorm(1 - alpha / 2))))
      ci_lwr <- quantile(boot_rho, probs = alpha1, na.rm = TRUE)
      ci_upr <- quantile(boot_rho, probs = alpha2, na.rm = TRUE)
    } else {
      stop("type must be 'perc', 'norm', or 'bca'.")
    }
    ci_lwr <- max(min(ci_lwr, 1), -1)
    ci_upr <- max(min(ci_upr, 1), -1)
  } else {
    ci_lwr <- ci_upr <- NA_real_
  }

  # --- Output ---
  out <- data.frame(
    rho = round(rho, 4),
    r_c = round(r_c, 4),
    ci.lower = if (ci) round(ci_lwr, 4) else NA,
    ci.upper = if (ci) round(ci_upr, 4) else NA,
    p.value = round(p_val, 6),
    method = p_method_label,
    alternative = alternative,
    stringsAsFactors = FALSE
  )
  return(out)
}
