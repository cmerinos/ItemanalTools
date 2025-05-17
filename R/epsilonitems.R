#' @title epsilonitems: Epsilon Squared (ε²) for Multiple-Group Item Comparison
#'
#' @description 
#' Computes the non-parametric effect size epsilon squared (ε²) for each ordinal item across multiple groups.
#' Based on Kruskal-Wallis test. Bootstrap confidence intervals can be computed using different methods.
#'
#' @param data.items A data frame of ordinal item responses (e.g., Likert-type).
#' @param group A grouping variable (factor, character, or numeric) with 3 or more levels.
#' @param ci Logical. Whether to compute confidence intervals. Default is TRUE.
#' @param B Integer. Number of bootstrap samples. Default is 1000.
#' @param type Character. Type of confidence interval: "perc" (default), "norm", "basic", or "bca".
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{Item}: Item name.
#'   \item \code{eps2}: Epsilon squared effect size.
#'   \item \code{lwr.ci}, \code{lwr.upp}: Lower and upper bounds of the CI (if \code{ci = TRUE}).
#' }
#'
#' @examples
#' set.seed(123)
#' items <- data.frame(Item1 = sample(1:5, 90, TRUE),
#'                     Item2 = sample(1:5, 90, TRUE))
#' grupo <- sample(c("A", "B", "C"), 90, replace = TRUE)
#' epsilonitems(items, grupo, ci = TRUE, B = 500, type = "perc")
#'
#' @export
epsilonitems <- function(data.items, group, ci = TRUE, B = 1000, type = "perc") {
  require(rcompanion)
  
  output <- lapply(names(data.items), function(x) {
    res <- epsilonSquared(x = data.items[[x]], 
                          g = group,
                          ci = ci,
                          type = type, 
                          R = B)
    
    if (ci) {
      data.frame(Item = x,
                 eps2 = res[1],
                 lwr.ci = res[2],
                 lwr.upp = res[3])
    } else {
      data.frame(Item = x,
                 eps2 = as.numeric(res))
    }
  })
  
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  return(result)
}
