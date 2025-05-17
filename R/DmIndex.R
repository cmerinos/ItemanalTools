#' @title DmIndex: Item-Level Deviation Index from a Theoretical Association Value
#'
#' @description
#' Computes the Dm index for each item, quantifying the deviation of observed associations 
#' from a reference value. Also flags whether the reference value falls within each item's confidence interval.
#' The Dm index is calculated as: \cr
#' \code{Dm = 1 - (r_observed - ref.value)^2} \cr
#' A Dm closer to 1 indicates that the observed coefficient is close to the expected value.
#'
#' @param df A data frame of item-level results, typically returned by functions such as \code{RGWil}, \code{SpearItems}, or \code{epsilonitems}.
#' @param coef.col Character. Name of the column containing the association coefficient (e.g., "rg", "tau", "eps2").
#' @param lwr.col Character. Name of the column with the lower bound of the confidence interval.
#' @param upr.col Character. Name of the column with the upper bound of the confidence interval.
#' @param ref.value Numeric. The reference correlation value to test against. Default is 0.30.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{Item}: The item name.
#'   \item \code{Dm}: Deviation index relative to the reference value.
#'   \item \code{solap}: "Yes" if the confidence interval includes the reference value; otherwise "No".
#' }
#'
#' @examples
#' result <- data.frame(
#'   Item = paste0("Item", 1:5),
#'   rg = c(.22, .18, .35, .12, .29),
#'   lwr.ci = c(.10, .05, .25, .00, .15),
#'   lwr.upp = c(.33, .29, .45, .20, .40)
#' )
#' DmIndex(df = result, coef.col = "rg", lwr.col = "lwr.ci", upr.col = "lwr.upp", ref.value = 0.30)
#'
#' @export
DmIndex <- function(df, coef.col, lwr.col, upr.col, ref.value = 0.30) {
  coef.vec <- df[[coef.col]]
  lwr.vec  <- df[[lwr.col]]
  upr.vec  <- df[[upr.col]]
  
  Dm <- 1 - (coef.vec - ref.value)^2
  solap <- ifelse(ref.value >= lwr.vec & ref.value <= upr.vec, "Yes", "No")
  
  out <- data.frame(
    Item = df$Item,
    Dm = round(Dm, 3),
    solap = solap
  )
  
  return(out)
}
