#' @title Item Thresholds, Monotonicity and Equidistance for Multiple Items
#'
#' @description
#' Applies \code{itemThresholds()} to multiple Likert-type items in a data
#' set and produces per-item summaries, a formatted table per item with
#' thresholds, threshold differences, monotonicity and equidistance
#' diagnostics, and an aggregate view at the scale level.
#'
#' @details
#' For each selected item (column) in \code{data}, this function calls
#' \code{itemThresholds()} and collects:
#' \itemize{
#'   \item Threshold estimates and standard errors.
#'   \item Confidence intervals for each threshold.
#'   \item Differences between consecutive thresholds and their z-tests
#'         (with standard errors that account for the covariance between
#'         thresholds).
#'   \item Monotonicity diagnostics: raw sign, z-test, one-tailed p-value,
#'         and a binary flag for each adjacent pair.
#'   \item Standardized effect sizes (Cohen's d) and a redundancy flag for
#'         each adjacent pair.
#'   \item Equidistance diagnostics: second differences, Z-tests, Wald test,
#'         and the Equidistance Index (Spratto).
#' }
#'
#' In addition, the function returns an aggregate view of the scale
#' (\code{scale_summary}) containing counts of disordered items and redundant
#' pairs, and summary statistics of the Equidistance Index and Cohen's d
#' across all items. This aggregate view is **descriptive**: it summarizes
#' the item-level diagnostics and does not assume a latent dimension
#' underlying the items.
#'
#' ## Formatted table per item
#'
#' The formatted table for each item has the following columns:
#' \describe{
#'   \item{row}{Label for the row (threshold name, contrast,
#'         \code{"Diff-Thr. Item"}, or \code{"Z.local (d*)"}).}
#'   \item{Threshold}{Threshold estimate (if applicable).}
#'   \item{SE.threshold}{Standard error of the threshold.}
#'   \item{Diff.thresholds}{Difference between consecutive thresholds.}
#'   \item{SE.diff}{Standard error of the difference (accounting for
#'         covariance between thresholds).}
#'   \item{Z}{z statistic. For contrast rows this is the z-test of the
#'         difference between consecutive thresholds. For equidistance rows
#'         this is the Z-test of the corresponding second difference.}
#'   \item{p.diff}{Two-tailed p-value for the difference.}
#'   \item{p.disord}{One-tailed p-value for the monotonicity test
#'         (\eqn{H_1: \Delta < 0}).}
#'   \item{d.cohen}{Standardized effect size (Cohen's d) for the difference.
#'         Interpretation depends on \code{link}: see \code{itemThresholds}.}
#'   \item{disordered}{1 if the raw difference is negative, 0 otherwise.}
#'   \item{Wald}{Global Wald statistic for equidistance (in the
#'         \code{"Diff-Thr. Item"} row).}
#'   \item{df}{Degrees of freedom of the Wald test.}
#'   \item{EI.spratto}{Equidistance Index (Spratto) for the item.}
#' }
#'
#' ## Z-tests of second differences (local equidistance tests)
#'
#' For items with four response categories (three thresholds), there is only
#' one second difference, and its Z-test is shown in the
#' \code{"Diff-Thr. Item"} row. For items with more than four categories,
#' the function adds one extra row per second difference
#' (\code{"Z.local (d1)"}, \code{"Z.local (d2)"}, ...). The global Wald
#' test and the Equidistance Index remain reported in the
#' \code{"Diff-Thr. Item"} row.
#'
#' @param data A \code{data.frame} or matrix containing Likert-type items.
#' @param items Optional character vector or numeric indices indicating
#'   which columns of \code{data} to analyse. If \code{NULL} (default), all
#'   columns are used.
#' @param link The link function to use in the cumulative link model.
#'   Passed to \code{itemThresholds()} (e.g., \code{"probit"},
#'   \code{"logit"}). Default is \code{"probit"}.
#' @param conf.level Confidence level for the threshold intervals. Passed to
#'   \code{itemThresholds()}. Default is 0.95.
#' @param nd Number of digits used to round numeric outputs. Default is 3.
#'
#' @return A list with four components:
#' \item{items}{Named list of full \code{itemThresholds()} results, one per
#'   item.}
#' \item{tables}{Named list of formatted tables (data frames), one per item.}
#' \item{summary}{Data frame with one row per item, containing
#'   \code{item}, \code{K}, \code{EI.spratto}, \code{Wald.stat}, \code{df},
#'   \code{p.value}, \code{any_disord}, and \code{n_redundant}.}
#' \item{scale_summary}{One-row data frame with aggregate diagnostics of the
#'   scale: number of items, disordered items, redundant pairs, and summary
#'   statistics of the Equidistance Index and Cohen's d.}
#'
#' @examples
#' ## Example 1: Item with 4 categories
#' set.seed(123)
#' x4 <- ordered(sample(1:4, 300, replace = TRUE,
#'                     prob = c(.2, .3, .3, .2)))
#' data4 <- data.frame(ItemA = x4)
#' res4 <- itemThresholdsMult(data4)
#' res4$summary
#' res4$scale_summary
#' res4$tables$ItemA
#'
#' ## Example 2: Item with 5 categories (multiple Z.local rows)
#' set.seed(456)
#' x5 <- ordered(sample(1:5, 400, replace = TRUE,
#'                     prob = c(.1, .2, .3, .25, .15)))
#' data5 <- data.frame(ItemB = x5)
#' res5 <- itemThresholdsMult(data5)
#' res5$summary
#' res5$scale_summary
#' res5$tables$ItemB
#'
#' @seealso \code{\link{itemThresholds}}
#'
#' @export
itemThresholdsMult <- function(data,
                               items = NULL,
                               link = "probit",
                               conf.level = 0.95,
                               nd = 3) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (is.null(items)) {
    items <- names(data)
  } else if (is.numeric(items)) {
    items <- names(data)[items]
  }

  round_df <- function(df, digits) {
    if (nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  items_results <- list()
  tables_list   <- list()
  summary_rows  <- list()

  # Accumulators for scale-level summary
  EI_vec            <- numeric(0)
  d_vec             <- numeric(0)
  n_disord_items    <- 0L
  n_redundant_pairs <- 0L
  n_ok_items        <- 0L

  for (nm in items) {
    x <- data[[nm]]

    if (all(is.na(x))) {
      warning(sprintf("Item '%s' contains only NA values. Skipping.", nm))
      next
    }

    res_item <- tryCatch(
      itemThresholds(x = x, link = link, conf.level = conf.level, nd = nd),
      error = function(e) {
        warning(sprintf("itemThresholds() failed for item '%s': %s", nm, e$message))
        return(NULL)
      }
    )

    if (is.null(res_item)) next

    items_results[[nm]] <- res_item

    thr <- res_item$thresholds
    dft <- res_item$diff_tests
    eq  <- res_item$equidistance
    mono <- res_item$monotonicity

    x_ord <- if (is.ordered(x)) x else ordered(x)
    K <- length(levels(x_ord))

    # ---------- Formatted table per item ----------

    # 1) Threshold rows
    tab_thr <- data.frame(
      row             = thr$threshold,
      Threshold       = thr$estimate,
      SE.threshold    = thr$se,
      Diff.thresholds = NA_real_,
      SE.diff         = NA_real_,
      Z               = NA_real_,
      p.diff          = NA_real_,
      p.disord        = NA_real_,
      d.cohen         = NA_real_,
      disordered      = NA_integer_,
      Wald            = NA_real_,
      df              = NA_integer_,
      EI.spratto      = NA_real_,
      stringsAsFactors = FALSE
    )

    # 2) Difference rows (with monotonicity and effect size)
    if (nrow(dft) > 0L) {
      tab_diff <- data.frame(
        row             = dft$contrast,
        Threshold       = NA_real_,
        SE.threshold    = NA_real_,
        Diff.thresholds = dft$difference,
        SE.diff         = dft$se.diff,
        Z               = dft$z,
        p.diff          = dft$p.value,
        p.disord        = dft$p.one.sided,
        d.cohen         = dft$d.cohen,
        disordered      = dft$disordered,
        Wald            = NA_real_,
        df              = NA_integer_,
        EI.spratto      = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      tab_diff <- data.frame(
        row = character(0), Threshold = numeric(0), SE.threshold = numeric(0),
        Diff.thresholds = numeric(0), SE.diff = numeric(0), Z = numeric(0),
        p.diff = numeric(0), p.disord = numeric(0), d.cohen = numeric(0),
        disordered = integer(0), Wald = numeric(0), df = integer(0),
        EI.spratto = numeric(0), stringsAsFactors = FALSE
      )
    }

    # 3) Global equidistance row
    second_diffs <- eq$second_diffs

    Z_eq <- if (!is.null(second_diffs) && nrow(second_diffs) == 1L) {
      second_diffs$Z[1]
    } else {
      NA_real_
    }

    if (!is.null(eq$wald) && nrow(eq$wald) > 0L) {
      Wald_stat <- eq$wald$stat[1]
      Wald_df   <- eq$wald$df[1]
      Wald_p    <- eq$wald$p.value[1]
    } else {
      Wald_stat <- NA_real_
      Wald_df   <- NA_integer_
      Wald_p    <- NA_real_
    }

    tab_eq <- data.frame(
      row             = "Diff-Thr. Item",
      Threshold       = NA_real_,
      SE.threshold    = NA_real_,
      Diff.thresholds = NA_real_,
      SE.diff         = NA_real_,
      Z               = Z_eq,
      p.diff          = NA_real_,
      p.disord        = NA_real_,
      d.cohen         = NA_real_,
      disordered      = NA_integer_,
      Wald            = Wald_stat,
      df              = Wald_df,
      EI.spratto      = eq$EI.spratto,
      stringsAsFactors = FALSE
    )

    # 4) Extra Z.local rows for K > 4
    extra_rows <- NULL
    if (!is.null(second_diffs) && nrow(second_diffs) > 1L) {
      for (j in seq_len(nrow(second_diffs))) {
        extra_rows <- rbind(
          extra_rows,
          data.frame(
            row             = paste0("Z.local (d", j, ")"),
            Threshold       = NA_real_,
            SE.threshold    = NA_real_,
            Diff.thresholds = NA_real_,
            SE.diff         = NA_real_,
            Z               = second_diffs$Z[j],
            p.diff          = NA_real_,
            p.disord        = NA_real_,
            d.cohen         = NA_real_,
            disordered      = NA_integer_,
            Wald            = NA_real_,
            df              = NA_integer_,
            EI.spratto      = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    tab_item <- rbind(tab_thr, tab_diff, tab_eq, extra_rows)
    tab_item <- round_df(tab_item, nd)
    tab_item[is.na(tab_item)] <- ""
    tables_list[[nm]] <- tab_item

    # ---------- Per-item summary ----------
    n_redundant_item <- if (nrow(dft) > 0L) sum(dft$redundant, na.rm = TRUE) else 0L

    summary_rows[[nm]] <- data.frame(
      item        = nm,
      K           = K,
      EI.spratto  = eq$EI.spratto,
      Wald.stat   = Wald_stat,
      df          = Wald_df,
      p.value     = Wald_p,
      any_disord  = res_item$any_disord,
      n_redundant = as.integer(n_redundant_item),
      stringsAsFactors = FALSE
    )

    # ---------- Accumulate for scale-level summary ----------
    n_ok_items <- n_ok_items + 1L
    if (!is.na(eq$EI.spratto)) EI_vec <- c(EI_vec, eq$EI.spratto)
    if (nrow(dft) > 0L) {
      d_valid <- dft$d.cohen[!is.na(dft$d.cohen)]
      if (length(d_valid) > 0L) d_vec <- c(d_vec, d_valid)
    }
    if (res_item$any_disord == 1L) n_disord_items <- n_disord_items + 1L
    n_redundant_pairs <- n_redundant_pairs + as.integer(n_redundant_item)
  }

  # ---------- Combine per-item summary ----------
  if (length(summary_rows) > 0L) {
    summary_df <- do.call(rbind, summary_rows)
    summary_df <- round_df(summary_df, nd)
  } else {
    summary_df <- data.frame(
      item = character(0), K = integer(0), EI.spratto = numeric(0),
      Wald.stat = numeric(0), df = integer(0), p.value = numeric(0),
      any_disord = integer(0), n_redundant = integer(0),
      stringsAsFactors = FALSE
    )
  }

  # ---------- Scale-level summary ----------
  safe_mean <- function(v) if (length(v) > 0L) mean(v, na.rm = TRUE) else NA_real_
  safe_med  <- function(v) if (length(v) > 0L) median(v, na.rm = TRUE) else NA_real_
  safe_min  <- function(v) if (length(v) > 0L) min(v, na.rm = TRUE) else NA_real_
  safe_max  <- function(v) if (length(v) > 0L) max(v, na.rm = TRUE) else NA_real_
  safe_sd   <- function(v) if (length(v) > 1L) sd(v, na.rm = TRUE) else NA_real_

  scale_ok <- as.integer(n_disord_items == 0L && n_redundant_pairs == 0L)

  scale_summary <- data.frame(
    n_items           = n_ok_items,
    n_disord_items    = n_disord_items,
    n_redundant_pairs = n_redundant_pairs,
    EI.mean           = safe_mean(EI_vec),
    EI.median         = safe_med(EI_vec),
    EI.min            = safe_min(EI_vec),
    EI.max            = safe_max(EI_vec),
    EI.sd             = safe_sd(EI_vec),
    d.mean            = safe_mean(d_vec),
    d.median          = safe_med(d_vec),
    d.min             = safe_min(d_vec),
    d.max             = safe_max(d_vec),
    scale_ok          = scale_ok,
    stringsAsFactors  = FALSE
  )

  scale_summary <- round_df(scale_summary, nd)
  # Ensure scale_ok remains integer after rounding
  scale_summary$scale_ok <- as.integer(scale_ok)

  return(list(
    items         = items_results,
    tables        = tables_list,
    summary       = summary_df,
    scale_summary = scale_summary
  ))
}
