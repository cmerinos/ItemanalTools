#' @title Item Thresholds for Multiple Likert-type Items
#'
#' @description
#' Applies \code{itemThresholds()} to multiple Likert-type items in a data set
#' and produces per-item summaries, including a "figure-style" table with
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
#'   \item A "figure-style" table per item with rows for thresholds,
#'         threshold differences, and an overall equidistance row.
#'   \item A summary data frame with one row per item, containing the number of
#'         categories, the Spratto index, and the global Wald test of
#'         equidistance.
#' }
#'
#' The "figure-style" table for each item has the following columns:
#' \describe{
#'   \item{row}{Label for the row (threshold name, contrast, or "Diff-Thr. Item").}
#'   \item{Threshold}{Threshold estimate (if applicable).}
#'   \item{SE.threshold}{Standard error of the threshold.}
#'   \item{Diff.thresholds}{Difference between consecutive thresholds
#'         (if applicable).}
#'   \item{SE.diff}{Standard error of the difference between thresholds.}
#'   \item{Z}{Z-test of equidistance (only filled when there is a single
#'            second difference; otherwise NA).}
#'   \item{Wald}{Global Wald statistic for equidistance.}
#'   \item{df}{Degrees of freedom of the Wald test.}
#'   \item{EI.spratto}{Equidistance Index (Spratto) for the item.}
#' }
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
#' @examples
#' \dontrun{
#' data(Science, package = "mirt")
#' res_mult <- itemThresholdsMult(data = bfi[, c("Comfort", "Work", "Future", "Benefit")])
#' res_mult$summary
#' res_mult$tables$Work
#' }
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

    # saltar si la columna es toda NA
    if (all(is.na(x))) {
      warning(sprintf("Item '%s' contains only NA values. Skipping.", nm))
      next
    }

    # intentar análisis del ítem
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

    # número de categorías (niveles) del ítem
    x_ord <- if (is.ordered(x)) x else ordered(x)
    K <- length(levels(x_ord))

    ## -------- Tabla "tipo figura" por ítem --------

    # 1) Filas de thresholds
    tab_thr <- data.frame(
      row           = thr$threshold,
      Threshold     = thr$estimate,
      SE.threshold  = thr$se,
      Diff.thresholds = NA_real_,
      SE.diff       = NA_real_,
      Z             = NA_real_,
      Wald          = NA_real_,
      df            = NA_integer_,
      EI.spratto    = NA_real_,
      stringsAsFactors = FALSE
    )

    # 2) Filas de diferencias entre thresholds
    if (nrow(dft) > 0L) {
      tab_diff <- data.frame(
        row            = dft$contrast,
        Threshold      = NA_real_,
        SE.threshold   = NA_real_,
        Diff.thresholds = dft$difference,
        SE.diff        = dft$se.diff,
        Z              = NA_real_,
        Wald           = NA_real_,
        df             = NA_integer_,
        EI.spratto     = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      tab_diff <- data.frame(
        row            = character(0),
        Threshold      = numeric(0),
        SE.threshold   = numeric(0),
        Diff.thresholds = numeric(0),
        SE.diff        = numeric(0),
        Z              = numeric(0),
        Wald           = numeric(0),
        df             = integer(0),
        EI.spratto     = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    # 3) Fila "Diff-Thr. Item" (equidistancia global)
    #    Z: solo si hay una única segunda diferencia; si hay varias, Z = NA
    second_diffs <- eq$second_diffs
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
      row            = "Diff-Thr. Item",
      Threshold      = NA_real_,
      SE.threshold   = NA_real_,
      Diff.thresholds = NA_real_,
      SE.diff        = NA_real_,
      Z              = Z_eq,
      Wald           = Wald_stat,
      df             = Wald_df,
      EI.spratto     = eq$EI.spratto,
      stringsAsFactors = FALSE
    )

    # Combinar tabla del ítem
    tab_item <- rbind(tab_thr, tab_diff, tab_eq)
    tab_item <- round_df(tab_item, nd)

    # Replace NA with blank for clean presentation
    tab_item[is.na(tab_item)] <- ""

    tables_list[[nm]] <- tab_item

    ## -------- Resumen por ítem --------
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

  # Combinar resumen
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
