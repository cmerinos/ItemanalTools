#' Chambers Excess‑Area Correlation for Biserial Data
#'
#' @param data.items A data frame or matrix with items in columns and subjects in rows.
#' @param criterion A numeric vector (continuous criterion) to be dichotomized at `cut.point`.
#' @param cut.point Numeric value used to dichotomize the criterion.
#' @param ci Logical. If `TRUE` (default), bootstrap CIs and a two‑sided p‑value are returned.
#' @param B Number of bootstrap replicates (if `ci = TRUE`). Default 1000.
#' @param type Method for bootstrap confidence intervals: `"perc"`, `"norm"`, or `"bca"`.
#'
#' @return A data frame with columns depending on `ci`.
#'
#' @references
#' Chambers, R. G. (1982). Correlation coefficients from 2×2 tables and from biserial data.
#' *British Journal of Mathematical and Statistical Psychology*, 35(2), 216–227.
#'
#' Alexander, R. A., Alliger, G. M., Carson, K. P., & Barrett, G. V. (1985).
#' The Empirical Performance of Measures of Association in the 2×2 Table.
#' *Educational and Psychological Measurement*, 45(1), 79–87.
#'
#' @examples
#' set.seed(123)
#' items <- data.frame(
#'   Item1 = sample(1:5, 100, replace = TRUE),
#'   Item2 = sample(1:5, 100, replace = TRUE)
#' )
#' total_score <- rowSums(items) + rnorm(100, 0, 2)
#' rChambers(items, total_score, cut.point = 10, ci = FALSE)
#' rChambers(items, total_score, cut.point = 10, ci = TRUE, B = 500)
#'
#' @export
rChambers <- function(data.items,
                      criterion,
                      cut.point,
                      ci = TRUE,
                      B = 1000,
                      type = "perc") {

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
  if (!is.numeric(criterion)) {
    stop("'criterion' must be numeric.")
  }
  if (!is.numeric(cut.point) || length(cut.point) != 1) {
    stop("'cut.point' must be a single numeric value.")
  }

  # Remove rows with NAs in items or criterion
  complete <- complete.cases(data.items, criterion)
  if (!all(complete)) {
    warning(sprintf("Removed %d rows with missing values.", sum(!complete)))
    data.items <- data.items[complete, , drop = FALSE]
    criterion <- criterion[complete]
  }

  if (nrow(data.items) < 3) {
    stop("At least 3 complete observations are required.")
  }
  if (sum(criterion >= cut.point) == 0 || sum(criterion < cut.point) == 0) {
    stop("'cut.point' does not split the criterion into both groups.")
  }

  # --- 2. Internal function to compute r.E ---
  calc_rE <- function(x) {
    ord <- order(x)
    x_ord <- x[ord]
    crit_ord <- criterion[ord]
    N <- length(x_ord)
    py <- mean(crit_ord >= cut.point)
    qy <- 1 - py
    if (py == 0 || py == 1) return(NA)
    n1 <- cumsum(crit_ord >= cut.point)
    AE <- sum(n1 - (1:N) * py) / N^2
    2 * AE / (py * qy)
  }

  r_vals <- sapply(data.items, calc_rE)
  df_out <- data.frame(Item = colnames(data.items),
                       r.E = r_vals,
                       stringsAsFactors = FALSE)

  # --- 3. Return if no CI ---
  if (isFALSE(ci)) {
    # Round numeric columns
    num_cols <- sapply(df_out, is.numeric)
    df_out[num_cols] <- lapply(df_out[num_cols], round, 3)
    return(df_out)
  }

  # --- 4. Bootstrap CIs and p-value ---
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for bootstrap CIs. Please install it.")
  }

  boot_rE <- function(data, indices, crit, cut) {
    d <- data[indices, , drop = FALSE]
    crit_b <- crit[indices]
    sapply(d, function(col) {
      ord <- order(col)
      x_ord <- col[ord]
      crit_ord <- crit_b[ord]
      N <- length(x_ord)
      py <- mean(crit_ord >= cut)
      qy <- 1 - py
      if (py == 0 || py == 1) return(NA)
      n1 <- cumsum(crit_ord >= cut)
      AE <- sum(n1 - (1:N) * py) / N^2
      2 * AE / (py * qy)
    })
  }

  boot_obj <- boot::boot(data = data.items,
                         statistic = boot_rE,
                         R = B,
                         crit = criterion,
                         cut = cut.point)

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
    obs <- df_out$r.E[i]
    mean(abs(boot_coefs[, i]) >= abs(obs), na.rm = TRUE)
  })
  df_out$p.value <- p_boot

  # Round numeric columns
  num_cols <- sapply(df_out, is.numeric)
  df_out[num_cols] <- lapply(df_out[num_cols], round, 3)

  return(df_out)
}
