#' Hall–Lambert Correlation for Ordinal Items and a Dichotomous Criterion
#'
#' @param data.items A data frame or matrix with ordinal items in columns and subjects in rows.
#' @param criterion A dichotomous vector (factor, logical, or numeric with two unique values).
#' @param lower.level Specifies which level of the criterion is considered the "lower" group.
#' @param ci Character or logical: `NULL`/`FALSE`, `"chi"`, or `"boot"`.
#' @param B Number of bootstrap replicates (if `ci = "boot"`). Default 1000.
#' @param type Method for bootstrap confidence intervals: `"perc"`, `"norm"`, or `"bca"`.
#' @param ties.method Method for handling ties in ranks. Passed to [rank()]. Default `"average"`.
#'
#' @return A data frame with columns depending on `ci`.
#'
#' @references
#' Hall, A. C., & Lambert, J. A. (1962). A modification of biserial r for the
#' correlation of a ranked variable with a dichotomous variable.
#' *The British Journal of Statistical Psychology*, 15(2), 129–132.
#'
#' @examples
#' set.seed(123)
#' items <- data.frame(
#'   Item1 = sample(1:5, 100, replace = TRUE),
#'   Item2 = sample(1:5, 100, replace = TRUE)
#' )
#' group <- sample(c(0, 1), 100, replace = TRUE)
#' rHallLambert(items, group)
#' rHallLambert(items, group, ci = "boot", B = 500)
#'
#' @importFrom stats pchisq
#' @export
rHallLambert <- function(data.items,
                         criterion,
                         lower.level = NULL,
                         ci = "chi",
                         B = 1000,
                         type = "perc",
                         ties.method = "average") {

  # --- 1. Validations ---
  if (!is.data.frame(data.items) && !is.matrix(data.items)) {
    stop("'data.items' must be a data frame or matrix.")
  }
  data.items <- as.data.frame(data.items)
  if (!all(sapply(data.items, is.numeric))) {
    stop("All columns in 'data.items' must be numeric.")
  }
  if (length(criterion) != nrow(data.items)) {
    stop("'criterion' must have the same length as the number of rows in 'data.items'.")
  }

  # --- 2. Process criterion into a binary 0/1 vector ---
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
  n <- sum(dic == 0)
  if (n == 0 || n == N) stop("The lower group cannot be empty or the whole sample.")

  # --- 3. Internal function to compute r.rd, H, and chi p-value ---
  calc_stats <- function(x) {
    r <- rank(x, ties.method = ties.method)
    Y <- sum(r[dic == 0])
    r_rd <- (n * (N + 1) - 2 * Y) / (n * (N - n))
    meanQ <- Y / n
    meanP <- (sum(r) - Y) / (N - n)
    H <- 12 * n * (N - n) * (meanP - meanQ)^2 / (N * (N + 1))
    p_chi <- pchisq(H, df = 1, lower.tail = FALSE)
    return(c(r_rd = r_rd, H = H, p_chi = p_chi))
  }

  res <- t(sapply(data.items, calc_stats))
  df_out <- data.frame(Item = rownames(res),
                       r.rd = res[, "r_rd"],
                       stringsAsFactors = FALSE)

  # --- 4. Return depending on ci ---
  if (is.null(ci) || isFALSE(ci)) {
    # Round numeric columns to 3 decimals
    num_cols <- sapply(df_out, is.numeric)
    df_out[num_cols] <- lapply(df_out[num_cols], round, 3)
    return(df_out)
  }

  if (ci == "chi") {
    df_out$H <- res[, "H"]
    df_out$p.value <- res[, "p_chi"]
    # Round numeric columns
    num_cols <- sapply(df_out, is.numeric)
    df_out[num_cols] <- lapply(df_out[num_cols], round, 3)
    return(df_out)
  }

  if (ci == "boot") {
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for bootstrap CIs. Please install it.")
    }

    boot_rrd <- function(data, indices, dic_orig) {
      d <- data[indices, , drop = FALSE]
      dic_b <- dic_orig[indices]
      N_b <- length(dic_b)
      n_b <- sum(dic_b == 0)
      if (n_b == 0 || n_b == N_b) return(rep(NA, ncol(d)))
      sapply(d, function(col) {
        r <- rank(col, ties.method = ties.method)
        Y <- sum(r[dic_b == 0])
        (n_b * (N_b + 1) - 2 * Y) / (n_b * (N_b - n_b))
      })
    }

    boot_obj <- boot::boot(data = data.items,
                           statistic = boot_rrd,
                           R = B,
                           dic_orig = dic)

    # Extract confidence intervals
    ic_list <- lapply(seq_len(ncol(data.items)), function(i) {
      ci_boot <- switch(type,
                        perc = boot::boot.ci(boot_obj, index = i, type = "perc"),
                        norm = boot::boot.ci(boot_obj, index = i, type = "norm"),
                        bca  = boot::boot.ci(boot_obj, index = i, type = "bca"),
                        stop("Invalid 'type'. Use 'perc', 'norm', or 'bca'."))
      if (type == "perc") c(lwr = ci_boot$perc[4], upr = ci_boot$perc[5])
      else if (type == "norm") c(lwr = ci_boot$norm[2], upr = ci_boot$norm[3])
      else c(lwr = ci_boot$bca[4], upr = ci_boot$bca[5])
    })
    ic_mat <- do.call(rbind, ic_list)
    df_out$lwr.ci <- ic_mat[, "lwr"]
    df_out$upr.ci <- ic_mat[, "upr"]

    # Two-sided bootstrap p-value
    boot_coefs <- boot_obj$t
    p_boot <- sapply(seq_len(ncol(data.items)), function(i) {
      obs <- df_out$r.rd[i]
      mean(abs(boot_coefs[, i]) >= abs(obs), na.rm = TRUE)
    })
    df_out$p.value <- p_boot

    # Round numeric columns
    num_cols <- sapply(df_out, is.numeric)
    df_out[num_cols] <- lapply(df_out[num_cols], round, 3)
    return(df_out)
  }

  stop("Invalid 'ci' argument. Use FALSE/NULL, 'chi', or 'boot'.")
}
