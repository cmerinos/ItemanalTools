#' @title summaryItemAssoc: Summary Statistics for Item-Level Association Coefficients
#'
#' @description
#' Computes descriptive statistics (mean, median) for a vector of item-level association coefficients
#' and, optionally, their confidence intervals.
#'
#' This function is designed to be used with the output of functions like \code{RGitems},
#' \code{SpearItems}, or \code{epsilonItems}, which return a data frame with item-wise statistics.
#'
#' @param df A data frame returned by an item association function.
#' @param coef.col Character. Name of the column containing the association coefficient (e.g., \code{"rg"}, \code{"rho"}, \code{"eps2"}).
#' @param lwr.col Optional. Name of the column with lower confidence interval limits.
#' @param upr.col Optional. Name of the column with upper confidence interval limits.
#'
#' @return A one-row data frame with the following columns:
#' \describe{
#'   \item{\code{coef.mean}}{Mean of the association coefficients.}
#'   \item{\code{coef.median}}{Median of the association coefficients.}
#'   \item{\code{lwr.mean}, \code{lwr.median}}{Mean and median of the lower bounds (if provided).}
#'   \item{\code{upr.mean}, \code{upr.median}}{Mean and median of the upper bounds (if provided).}
#' }
#'
#' @examples
#' df <- data.frame(
#'   Item = paste0("Item", 1:5),
#'   rg = c(.22, .18, .35, .12, .29),
#'   lwr.ci = c(.10, .05, .25, .00, .15),
#'   lwr.upp = c(.33, .29, .45, .20, .40)
#' )
#' summaryItemAssoc(df = df, coef.col = "rg", lwr.col = "lwr.ci", upr.col = "lwr.upp")
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
