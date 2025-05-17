#' @title FWitems: Friedman Test + Kendall's W for Items
#'
#' @description
#' Computes the Friedman test and Kendall’s W to assess differences and concordance among multiple ordinal items.
#' This function assumes a repeated-measures structure (i.e., each row is a subject, columns are items).
#'
#' @param data.items A data frame or matrix of ordinal items (each column is an item, each row is a subject).
#' @param ci Logical. Whether to compute confidence intervals for Kendall's W. Default is TRUE.
#' @param correct Logical. Apply correction for ties in Kendall's W. Default is TRUE.
#' @param type Character. Type of confidence interval ("norm", "basic", "perc", or "bca"). Default is "norm".
#' @param B Number of bootstrap replicates. Default is 1000.
#' @param conf.level Confidence level for the confidence intervals. Default is 0.95.
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{Friedman}: Output of the \code{friedman.test}.
#'   \item \code{KendallW}: A data frame with W and its confidence interval.
#' }
#'
#' @examples
#' 
#' ### Example 1 ---------------------
#' set.seed(123)
#' items <- data.frame(
#'   Item1 = sample(1:5, 50, TRUE),
#'   Item2 = sample(1:5, 50, TRUE),
#'   Item3 = sample(1:5, 50, TRUE)
#' )
#' FWitems(items)
#'
#' ### Example 2 ---------------------
#'
#' # Ejemplo con tendencia de respuesta (mayor concordancia esperada)
#' set.seed(456)
#' responses <- data.frame(
#'   Item1 = sample(c(1, 2, 3), 30, replace = TRUE, prob = c(0.1, 0.2, 0.7)),
#'   Item2 = sample(c(1, 2, 3), 30, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
#'   Item3 = sample(c(1, 2, 3), 30, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
#'   Item4 = sample(c(1, 2, 3), 30, replace = TRUE, prob = c(0.2, 0.4, 0.4)),
#'   Item5 = sample(c(1, 2, 3), 30, replace = TRUE, prob = c(0.6, 0.3, 0.1)))
#'   
#' FWitems(responses, ci = TRUE, B = 500)
#' 
#' @export
FWitems <- function(data.items, ci = TRUE, correct = TRUE, type = "norm", B = 1000, conf.level = 0.95) {
  # Validar clase
  if (!is.data.frame(data.items) && !is.matrix(data.items)) {
    stop("`data.items` must be a data frame or matrix.")
  }
  
  # Asegurar formato matriz
  data.mat <- as.matrix(data.items)
  
  # Friedman test
  friedman <- friedman.test(data.mat)
  
  # Kendall's W
  kw <- rcompanion::kendallW(data.mat,
                             correct = correct,
                             ci = ci,
                             conf = conf.level,
                             type = type,
                             R = B,
                             histogram = FALSE)
  
  # Preparar salida W
  if (ci) {
    W.out <- data.frame(W = round(kw$W, 3),
                        lwr.ci = round(kw$lower.ci, 3),
                        upp.ci = round(kw$upper.ci, 3))
  } else {
    W.out <- data.frame(W = round(kw$W, 3))
  }
  
  return(list(Friedman = friedman, KendallW = W.out))
}
