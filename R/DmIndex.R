#' @title DmIndex: Item-Level Deviation Index from a Theoretical Association Value
#'
#' @description
#' Computes the Dm index for each item, which quantifies how much an observed coefficient
#' deviates from a reference value. It also flags whether the reference value falls within
#' each item's confidence interval. The Dm index is calculated as:
#' \code{Dm = 1 - (r_observed - ref.value)^2}.
#'
#' A Dm value closer to 1 indicates that the observed coefficient is close to the expected value.
#'
#'
#' @param df A data frame of item-level results, typically returned by functions such as
#' \code{RGitems()}, \code{SpearItems()}, or \code{epsilonitems()}.
#' @param coef.col Character. Name of the column containing the association coefficient
#' (e.g., "rg", "tau", "eps2").
#' @param lwr.col Character. Name of the column containing the lower bound of the confidence interval.
#' @param upr.col Character. Name of the column containing the upper bound of the confidence interval.
#' @param ref.value Numeric. Reference value to compare against. Default is \code{0.30}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{\code{Item}}{The item name.}
#'   \item{\code{Dm}}{Deviation index relative to the reference value.}
#'   \item{\code{Overlap}}{"Yes" if the confidence interval includes the reference value; otherwise "No".}
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
#' @references
#' Taras, V., & Kline, T. (2010). Scale validation via quantifying item validity using
#' the Dm index. Psychological reports, 107(2), 535–546. https://doi.org/10.2466/03.PR0.107.5.535-546
#'
#' Dominguez Lara, S. (2013). Proposal for the calculation of the DM index for the validity of the items.
#' Interdisciplinaria Revista de Psicología y Ciencias Afines, 30(2), 297–303.
#' https://doi.org/10.16888/interd.2013.30.2.8
#'
#' Dominguez-Lara, S. (2018). Testing a general hypothesis in empirical research using Dm-index.
#' Interacciones, 4(1), 42-47. https://doi.org/10.24016/2018.v4n1.77



#' @export
DmIndex <- function(df, coef.col, lwr.col, upr.col, ref.value = 0.30) {
  coef.vec <- df[[coef.col]]
  lwr.vec  <- df[[lwr.col]]
  upr.vec  <- df[[upr.col]]

  Dm <- 1 - (coef.vec - ref.value)^2
  Overlap <- ifelse(ref.value >= lwr.vec & ref.value <= upr.vec, "Yes", "No")

  out <- data.frame(
    Item = df$Item,
    Dm = round(Dm, 3),
    Overlap = Overlap
  )

  return(out)
}
