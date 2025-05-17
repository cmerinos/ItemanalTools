#' @title summaryItemAssoc: Summary statistics for item association coefficients
#'
#' @description
#' Computes the mean and median of item-level coefficients (and optionally their confidence intervals).
#'
#' @param df Data frame returned by an ItemAn function (e.g., RGWil, SpearItems, epsilonitems).
#' @param coef.col Name of the column with the association coefficient (e.g., "rg", "tau", "eps2").
#' @param lwr.col Optional. Name of the column with lower CI bound.
#' @param upr.col Optional. Name of the column with upper CI bound.
#'
#' @return A data frame with summary statistics (mean, median) for the coefficient and CIs (if given).
#'
#' @examples
#' result <- data.frame(
#'   Item = c("Item1", "Item2", "Item3"),
#'   rg = c(.21, .35, .28),
#'   lwr.ci = c(.12, .22, .15),
#'   lwr.upp = c(.31, .45, .38)
#' )
#' summaryItemAssoc(df = result, coef.col = "rg", lwr.col = "lwr.ci", upr.col = "lwr.upp")
#'
#' @export
summaryItemAssoc <- function(df, coef.col, lwr.col = NULL, upr.col = NULL) {
  out <- data.frame(
    coef.mean = mean(df[[coef.col]], na.rm = TRUE),
    coef.median = median(df[[coef.col]], na.rm = TRUE)
  )

  if (!is.null(lwr.col) && !is.null(upr.col)) {
    out$lwr.mean <- mean(df[[lwr.col]], na.rm = TRUE)
    out$lwr.median <- median(df[[lwr.col]], na.rm = TRUE)
    out$upr.mean <- mean(df[[upr.col]], na.rm = TRUE)
    out$upr.median <- median(df[[upr.col]], na.rm = TRUE)
  }

  return(round(out, 3))
}
