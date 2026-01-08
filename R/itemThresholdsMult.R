#' @title Item Thresholds and equidistance for Multiple Likert-type Items
#'
#' @description
#' Applies \code{itemThresholds()} to multiple Likert-type items in a data set
#' and produces per-item summaries, including a formated table with
#' thresholds, threshold differences, and equidistance tests (Spratto index).
#'
#' @details
#' For each selected item (column) in \code{data}, this function calls
#' \code{itemThresholds()} and collects:
#' \itemize{
#'   \item Threshold estimates and standard errors.
#'   \item Confidence intervals for each threshold.
#'   \item Differences between consecutive thresholds and their z-tests.
#'   \item Equidistance diagnostics: second differences, Z-tests, Wald test,
#'         and the Equidistance Index (Spratto).
#' }
#'
#' The output includes:
#' \itemize{
#'   \item A list of full \code{itemThresholds()} results per item.
#'   \item A formated table per item with rows for thresholds,
#'         threshold differences, and an overall equidistance row.
#'   \item A summary data frame with one row per item, containing the number of
#'         categories, the Spratto index, and the global Wald test of
#'         equidistance.
#' }
#'
#' The formated table for each item has the following columns:
#' \describe{
#'   \item{row}{Label for the row (threshold name, contrast, "Diff-Thr. Item",
#'         or "Z.local (d1)", "Z.local (d2)", etc.).}
#'   \item{Threshold}{Threshold estimate (if applicable).}
#'   \item{SE.threshold}{Standard error of the threshold.}
#'   \item{Diff.thresholds}{Difference between consecutive thresholds
#'         (if applicable).}
#'   \item{SE.diff}{Standard error of the difference between thresholds.}
#'   \item{Z}{Z-test of equidistance. For items with four categories (i.e.,
#'         three thresholds), this column is used in the "Diff-Thr. Item" row.
#'         For items with more than four categories, local Z-tests are shown
#'         in additional "Z.local (d*)" rows.}
#'   \item{Wald}{Global Wald statistic for equidistance (in the
#'         "Diff-Thr. Item" row).}
#'   \item{df}{Degrees of freedom of the Wald test.}
#'   \item{EI.spratto}{Equidistance Index (Spratto) for the item.}
#' }
#'
#' ## Z-tests of second differences (local equidistance tests)
#'
#' For items with four response categories (i.e., three thresholds), there is
#' only one second difference, so the Z-test of equidistance is shown directly
#' in the "Diff-Thr. Item" row.
#'
#' For items with more than four categories (i.e., more than one second
#' difference), the function adds one extra row per second difference:
#'
#' \itemize{
#'   \item \code{Z.local (d1)}
#'   \item \code{Z.local (d2)}
#'   \item \code{Z.local (d3)}, etc.
#' }
#'
#' These rows display the individual Z-tests of equidistance for each second
#' difference. The global Wald test and the Equidistance Index (Spratto)
#' remain reported in the "Diff-Thr. Item" row.
#'
#' This format mirrors the presentation style found in Spratto (2018) and
#' similar analyses of threshold equidistance.
#'
#' @param data A \code{data.frame} or matrix containing Likert-type items
#'   (preferably ordered factors, but numeric or factor variables are also
#'   accepted and coerced).
#' @param items Optional character vector or numeric indices indicating which
#'   columns of \code{data} to analyse. If \code{NULL} (default), all columns
#'   are used.
#' @param link The link function to use in the cumulative link model. Passed to
#'   \code{itemThresholds()} (e.g., "probit", "logit").
#' @param conf.level Confidence level for the threshold intervals. Passed to
#'   \code{itemThresholds()}. Default is 0.95.
#' @param nd Number of digits used to round numeric outputs. Default is 3.
#'
#' @return A list with three components:
#' \item{items}{Named list of full \code{itemThresholds()} results, one per
#'   item.}
#' \item{tables}{Named list of "figure-style" tables (data frames), one per
#'   item.}
#' \item{summary}{Data frame with one row per item, containing:
#'   \code{item}, \code{K} (number of response categories),
#'   \code{EI.spratto}, \code{Wald.stat}, \code{df}, and \code{p.value}.}
#'
#' @references
#' Christensen, R. H. B. (2019). \emph{ordinal: Regression Models for
#' Ordinal Data}. R package version 2019.12-10.
#' https://CRAN.R-project.org/package=ordinal
#'
#' Spratto, E. M. (2018). \emph{In search of equality: Developing an equal interval Likert response scale
#' (Doctoral dissertations)}. https://commons.lib.jmu.edu/diss201019/172
#'
#' Sideridis, G., Tsaousis, I., & Ghamdi, H. (2022). Equidistant Response Options on Likert-Type
#' Instruments: Testing the Interval Scaling Assumption Using Mplus. \emph{Educational and Psychological Measurement},
#' 83(5), 885-906. https://doi.org/10.1177/00131644221130482
#'
#' @examples
#' ## Example 1: Item with 4 categories
#' set.seed(123)
#' x4 <- ordered(sample(1:4, 300, replace = TRUE,
#'                     prob = c(.2, .3, .3, .2)))
#' data4 <- data.frame(ItemA = x4)
#' res4 <- itemThresholdsMult(data4)
#' res4$summary
#' res4$tables$ItemA
#'
#' ## Example 2: Item with 5 categories (multiple Z.local rows)
#' set.seed(456)
#' x5 <- ordered(sample(1:5, 400, replace = TRUE,
#'                     prob = c(.1, .2, .3, .25, .15)))
#' data5 <- data.frame(ItemB = x5)
#' res5 <- itemThresholdsMult(data5)
#' res5$summary
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
  # basic checks and coercion
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (is.null(items)) {
    items <- names(data)
  } else {
    # allow numeric indices
    if (is.numeric(items)) {
      items <- names(data)[items]
    }
  }

  # helper to round data.frames
  round_df <- function(df, digits) {
    if (nrow(df) == 0L) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, digits = digits)
    df
  }

  items_results <- list()
  tables_list   <- list()
  summary_rows  <- list()

  for (nm in items) {
    x <- data[[nm]]

    # skip if all NA
    if (all(is.na(x))) {
      warning(sprintf("Item '%s' contains only NA values. Skipping.", nm))
      next
    }

    # analyse item
    res_item <- tryCatch(
      itemThresholds(x = x, link = link, conf.level = conf.level, nd = nd),
      error = function(e) {
        warning(sprintf("itemThresholds() failed for item '%s': %s", nm, e$message))
        return(NULL)
      }
    )

    if (is.null(res_item)) next

    items_results[[nm]] <- res_item

    thr  <- res_item$thresholds
    dft  <- res_item$diff_tests
    eq   <- res_item$equidistance

    # number of categories
    x_ord <- if (is.ordered(x)) x else ordered(x)
    K <- length(levels(x_ord))

    ## -------- Figure-style table per item --------

    # 1) Threshold rows
    tab_thr <- data.frame(
      row             = thr$threshold,
      Threshold       = thr$estimate,
      SE.threshold    = thr$se,
      Diff.thresholds = NA_real_,
      SE.diff         = NA_real_,
      Z               = NA_real_,
      Wald            = NA_real_,
      df              = NA_integer_,
      EI.spratto      = NA_real_,
      stringsAsFactors = FALSE
    )

    # 2) Difference rows
    if (nrow(dft) > 0L) {
      tab_diff <- data.frame(
        row             = dft$contrast,
        Threshold       = NA_real_,
        SE.threshold    = NA_real_,
        Diff.thresholds = dft$difference,
        SE.diff         = dft$se.diff,
        Z               = NA_real_,
        Wald            = NA_real_,
        df              = NA_integer_,
        EI.spratto      = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      tab_diff <- data.frame(
        row             = character(0),
        Threshold       = numeric(0),
        SE.threshold    = numeric(0),
        Diff.thresholds = numeric(0),
        SE.diff         = numeric(0),
        Z               = numeric(0),
        Wald            = numeric(0),
        df              = integer(0),
        EI.spratto      = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    # 3) Global equidistance row "Diff-Thr. Item"
    second_diffs <- eq$second_diffs

    # Z in Diff-Thr. Item only if there is a single second difference
    if (!is.null(second_diffs) && nrow(second_diffs) == 1L) {
      Z_eq <- second_diffs$Z[1]
    } else {
      Z_eq <- NA_real_
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
      Wald            = Wald_stat,
      df              = Wald_df,
      EI.spratto      = eq$EI.spratto,
      stringsAsFactors = FALSE
    )

    # 4) Extra Z.local rows for items with multiple second differences (K > 4)
    extra_rows <- NULL
    if (!is.null(second_diffs) && nrow(second_diffs) > 1L) {
      for (j in 1:nrow(second_diffs)) {
        extra_rows <- rbind(
          extra_rows,
          data.frame(
            row             = paste0("Z.local (d", j, ")"),
            Threshold       = NA_real_,
            SE.threshold    = NA_real_,
            Diff.thresholds = NA_real_,
            SE.diff         = NA_real_,
            Z               = second_diffs$Z[j],
            Wald            = NA_real_,
            df              = NA_integer_,
            EI.spratto      = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    # Combine table
    tab_item <- rbind(tab_thr, tab_diff, tab_eq, extra_rows)
    tab_item <- round_df(tab_item, nd)

    # Replace NA with blanks for cleaner presentation
    tab_item[is.na(tab_item)] <- ""

    tables_list[[nm]] <- tab_item

    ## -------- Summary per item --------
    summary_rows[[nm]] <- data.frame(
      item       = nm,
      K          = K,
      EI.spratto = eq$EI.spratto,
      Wald.stat  = Wald_stat,
      df         = Wald_df,
      p.value    = Wald_p,
      stringsAsFactors = FALSE
    )
  }

  # Combine summary
  if (length(summary_rows) > 0L) {
    summary_df <- do.call(rbind, summary_rows)
    summary_df <- round_df(summary_df, nd)
  } else {
    summary_df <- data.frame(
      item       = character(0),
      K          = integer(0),
      EI.spratto = numeric(0),
      Wald.stat  = numeric(0),
      df         = integer(0),
      p.value    = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  return(list(
    items   = items_results,
    tables  = tables_list,
    summary = summary_df
  ))
}
