#' @title Epsilon Squared (ε²) for Multiple-Group Item Comparison
#'
#' @description
#' Computes the nonparametric effect size epsilon squared (ε²) for each ordinal item
#' across multiple groups, based on the Kruskal-Wallis test. Bootstrap confidence intervals
#' can be computed using different methods.
#'
#' Epsilon squared is an interpretable metric similar to adjusted R² in linear models.
#' They are often treated as equivalent in magnitude and interpretation (Frieman et al., 2017; Vogt, 2005).
#'
#' @param data.items A data frame of ordinal item responses (e.g., Likert-type).
#' @param group A grouping variable (factor, character, or numeric) with 3 or more levels.
#' @param ci Logical. Whether to compute confidence intervals. Default is \code{TRUE}.
#' @param B Integer. Number of bootstrap samples. Default is \code{1000}.
#' @param type Character. Type of confidence interval: one of \code{"perc"} (default),
#' \code{"norm"}, \code{"basic"}, or \code{"bca"}.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{Item}: Item name.
#'   \item \code{eps2}: Epsilon squared effect size.
#'   \item \code{lwr.ci}, \code{lwr.upp}: Lower and upper bounds of the confidence interval (if \code{ci = TRUE}).
#' }
#'
#'#' @details
#' This function relies on \code{rcompanion::epsilonSquared} to compute
#' epsilon-squared, a nonparametric measure of effect size suitable for item-level
#' comparisons.
#'
#' @examples
#' set.seed(123)
#' items <- data.frame(
#'   Item1 = sample(1:5, 90, TRUE),
#'   Item2 = sample(1:5, 90, TRUE)
#' )
#' group <- sample(c("A", "B", "C"), 90, replace = TRUE)
#' epsilonItems(items, group, ci = TRUE, B = 500, type = "perc")
#'
#' @references
#' Frieman, J., Saucier, D. A., & Miller, S. (2017). Principles & Methods of Statistical Analysis.
#'
#' Carroll, R. M., & Nordholm, L. A. (1975). Sampling characteristics of Kelley’s ε² and Hays’ ω².
#' Educational and Psychological Measurement, 35, 541-554.
#'
#' Kelley, T. L. (1935). An unbiased correlation ratio measure. Proceedings of the National Academy of Sciences, 21(9), 554-559.
#'
#' Vogt, W. P. (2005). Dictionary of Statistics & Methodology: A Nontechnical Guide for the Social Sciences. SAGE Publications.
#'
#' @export
epsilonItems <- function(data.items, group, ci = TRUE, B = 1000, type = "perc") {
  if (!requireNamespace("rcompanion", quietly = TRUE)) {
    stop("Package 'rcompanion' is required.")
  }

  output <- lapply(names(data.items), function(x) {
    res <- rcompanion::epsilonSquared(x = data.items[[x]],
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
