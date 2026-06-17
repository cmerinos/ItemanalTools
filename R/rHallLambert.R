#' @title Hall-Lambert Correlation between Ordinal Items and a Dichotomous Criterion
#'
#' @description
#' Computes the \eqn{r_{rd}} correlation coefficient proposed by Hall and Lambert (1962)
#' to measure the association between an ordinal variable (e.g., Likert-type items)
#' and a dichotomous criterion. The coefficient does not assume bivariate normality
#' and is appropriate when the underlying continuous variable has a rectangular
#' distribution or when only ranks are available.
#'
#' @param data.items A data frame or matrix where each column is an item (ordinal variable)
#'        and each row is a respondent. All columns must be numeric.
#' @param criterion A vector of the same length as the number of rows in
#'        \code{data.items}, containing the dichotomous criterion. It can be a
#'        two-level factor, a logical vector, or a numeric vector with values 0/1.
#' @param lower.level Specifies which level of the criterion is considered the
#'        "lower" group. By default, for factors the first level is taken; for
#'        logical vectors, \code{FALSE}; for numeric vectors, the smallest value.
#'        Can be overridden with a specific value (e.g., \code{lower.level = 0}
#'        or \code{lower.level = "Control"}).
#' @param ci Logical. If \code{TRUE}, bootstrap confidence intervals are computed
#'        for the coefficient. Default is \code{FALSE}.
#' @param B Number of bootstrap replicates for the confidence intervals. Default is 1000.
#' @param type Method for the confidence interval: \code{"perc"} (percentile, default),
#'        \code{"norm"} (normal approximation), or \code{"bca"} (bias-corrected and accelerated).
#' @param ties.method Method for handling ties in the item ranks. Passed to \code{\link{rank}}.
#'        Default is \code{"average"}.
#'
#' @return A data frame with one row per item and the following columns:
#' \itemize{
#'   \item \code{Item}: Item name (taken from \code{colnames(data.items)}).
#'   \item \code{r.rd}: Hall-Lambert correlation coefficient (range approximately -1 to 1).
#'   \item \code{H}: Chi-square statistic for testing \eqn{H_0: \rho_{rd}=0}.
#'   \item \code{p.value}: p-value associated with the H statistic (\eqn{\chi^2} with 1 df).
#'   \item \code{lwr.ci}: Lower bound of the confidence interval (if \code{ci=TRUE}).
#'   \item \code{upr.ci}: Upper bound of the confidence interval (if \code{ci=TRUE}).
#' }
#' All numeric values are rounded to three decimal places.
#'
#' @details
#' The coefficient \eqn{r_{rd}} is based on the following formula (equation 6a in Hall & Lambert, 1962):
#' \deqn{r_{rd} = \frac{n(N+1) - 2Y}{n(N-n)}}
#' where \eqn{N} is the total number of observations, \eqn{n} is the size of the
#' lower group on the criterion, and \eqn{Y} is the sum of the ranks (with possible ties)
#' of the subjects in that group.
#'
#' Significance is assessed using the H statistic (equation 10):
#' \deqn{H = \frac{12 n (N-n) (\bar{P} - \bar{Q})^2}{N(N+1)}}
#' which, under the null hypothesis of independence, is approximately distributed as
#' \eqn{\chi^2} with 1 degree of freedom.
#'
#' The coefficient is designed for situations where the ordinal variable (item)
#' can be treated as a series of ranks and the dichotomous variable is an external
#' criterion. It is a robust alternative to the point-biserial and tetrachoric
#' coefficients, especially when normality assumptions are not met.
#'
#' @references
#' Hall, A. C., & Lambert, J. A. (1962). A modification of biserial r for the
#' correlation of a ranked variable with a dichotomous variable. \emph{The British
#' Journal of Statistical Psychology}, 15(2), 129–132.
#' \doi{10.1111/j.2044-8317.1962.tb00237.x}
#'
#' @seealso \code{\link{SpearItems}} for correlation between items and a continuous variable,
#'          \code{\link{RotheryItems}} for concordance among items.
#'
#' @examples
#' if (requireNamespace("boot", quietly = TRUE)) {
#'   # Simulated data: 100 respondents, 3 Likert items (1-5), dichotomous criterion (0/1)
#'   set.seed(123)
#'   items <- data.frame(
#'     Item1 = sample(1:5, 100, replace = TRUE),
#'     Item2 = sample(1:5, 100, replace = TRUE),
#'     Item3 = sample(1:5, 100, replace = TRUE)
#'   )
#'   # Criterion: group (0 = control, 1 = treatment)
#'   group <- sample(c(0,1), 100, replace = TRUE)
#'
#'   # Without confidence intervals (fast)
#'   rHallLambert(items, group)
#'
#'   # With bootstrap confidence intervals (slower)
#'   rHallLambert(items, group, ci = TRUE, B = 500)
#' }
#'
#' @importFrom stats pchisq qnorm quantile
#' @importFrom boot boot boot.ci
#' @export
rHallLambert <- function(data.items,
                         criterion,
                         lower.level = NULL,
                         ci = FALSE,
                         B = 1000,
                         type = "perc",
                         ties.method = "average") {

  # --- 1. Initial validation ---
  if (!is.data.frame(data.items) && !is.matrix(data.items)) {
    stop("'data.items' must be a data frame or a matrix.")
  }
  data.items <- as.data.frame(data.items)

  if (!all(sapply(data.items, is.numeric))) {
    stop("All columns in 'data.items' must be numeric.")
  }

  if (length(criterion) != nrow(data.items)) {
    stop("The criterion must have the same length as the number of rows in 'data.items'.")
  }

  # --- 2. Process criterion to a dichotomous 0/1 vector ---
  # and define the "lower" group
  if (is.factor(criterion)) {
    if (nlevels(criterion) != 2) {
      stop("The criterion must be a factor with exactly two levels.")
    }
    if (is.null(lower.level)) {
      lower.level <- levels(criterion)[1]
    }
    dic <- ifelse(criterion == lower.level, 0, 1)
    if (!any(dic == 0) || !any(dic == 1)) {
      stop("The level specified in 'lower.level' is not present in the criterion.")
    }
  } else if (is.logical(criterion)) {
    if (is.null(lower.level)) lower.level <- FALSE
    dic <- ifelse(criterion == lower.level, 0, 1)
    if (!any(dic == 0) || !any(dic == 1)) {
      stop("The value specified in 'lower.level' is not present in the criterion.")
    }
  } else if (is.numeric(criterion)) {
    if (length(unique(criterion)) != 2) {
      stop("The numeric criterion must have exactly two unique values.")
    }
    if (is.null(lower.level)) {
      lower.level <- min(criterion, na.rm = TRUE)
    }
    dic <- ifelse(criterion == lower.level, 0, 1)
    if (!any(dic == 0) || !any(dic == 1)) {
      stop("The value specified in 'lower.level' is not present in the criterion.")
    }
  } else {
    stop("The criterion must be a two-level factor, a logical vector, or a numeric vector with two values.")
  }

  # --- 3. Remove rows with missing values ---
  complete <- complete.cases(data.items, dic)
  if (!all(complete)) {
    warning(sprintf("Removed %d rows with missing values.", sum(!complete)))
    data.items <- data.items[complete, , drop = FALSE]
    dic <- dic[complete]
  }

  N <- nrow(data.items)
  if (N < 3) {
    stop("At least 3 complete observations are required.")
  }
  n <- sum(dic == 0)
  if (n == 0 || n == N) {
    stop("The lower group cannot be empty or the entire sample.")
  }

  # --- 4. Helper function to compute r.rd, H, and p-value for one item ---
  calc_rrd <- function(x) {
    r <- rank(x, ties.method = ties.method)
    total_rank_sum <- sum(r)
    Y <- sum(r[dic == 0])
    meanQ <- Y / n
    meanP <- (total_rank_sum - Y) / (N - n)

    r_rd <- (n * (N + 1) - 2 * Y) / (n * (N - n))

    H <- 12 * n * (N - n) * (meanP - meanQ)^2 / (N * (N + 1))
    p_valor <- pchisq(H, df = 1, lower.tail = FALSE)

    return(c(r_rd = r_rd, H = H, p.value = p_valor))
  }

  # --- 5. Apply to each item ---
  resultados <- t(sapply(data.items, calc_rrd))
  df_out <- data.frame(
    Item = rownames(resultados),
    r.rd = resultados[, "r_rd"],
    H = resultados[, "H"],
    p.value = resultados[, "p.value"],
    stringsAsFactors = FALSE
  )

  # --- 6. Bootstrap confidence intervals (if requested) ---
  if (ci) {
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("The 'boot' package is required for confidence intervals. Please install it.")
    }

    boot_rrd <- function(data, indices, items, dic_orig) {
      d <- data[indices, , drop = FALSE]
      dic_b <- dic_orig[indices]
      sapply(d, function(col) {
        r <- rank(col, ties.method = ties.method)
        Y <- sum(r[dic_b == 0])
        n_b <- sum(dic_b == 0)
        N_b <- length(dic_b)
        if (n_b == 0 || n_b == N_b) return(NA)
        (n_b * (N_b + 1) - 2 * Y) / (n_b * (N_b - n_b))
      })
    }

    boot_obj <- boot::boot(data = data.items,
                           statistic = boot_rrd,
                           R = B,
                           items = data.items,
                           dic_orig = dic)

    ic_list <- lapply(seq_len(ncol(data.items)), function(i) {
      if (type == "perc") {
        ci_boot <- boot::boot.ci(boot_obj, index = i, type = "perc")
        c(lwr = ci_boot$perc[4], upr = ci_boot$perc[5])
      } else if (type == "norm") {
        ci_boot <- boot::boot.ci(boot_obj, index = i, type = "norm")
        c(lwr = ci_boot$norm[2], upr = ci_boot$norm[3])
      } else if (type == "bca") {
        ci_boot <- boot::boot.ci(boot_obj, index = i, type = "bca")
        c(lwr = ci_boot$bca[4], upr = ci_boot$bca[5])
      } else {
        stop("Unrecognised CI type. Use 'perc', 'norm', or 'bca'.")
      }
    })
    ic_mat <- do.call(rbind, ic_list)
    df_out$lwr.ci <- ic_mat[, "lwr"]
    df_out$upr.ci <- ic_mat[, "upr"]
  }

  # --- 7. Round numeric columns to 3 decimal places ---
  numeric_cols <- sapply(df_out, is.numeric)
  df_out[numeric_cols] <- lapply(df_out[numeric_cols], round, 3)

  rownames(df_out) <- NULL
  return(df_out)
}
