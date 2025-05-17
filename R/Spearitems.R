#' @title SpearItems: Correlation Between Items and External Criterion (Nonparametric)
#'
#' @description
#' Computes Spearman's rho or Kendall's tau correlation between each ordinal item and
#' a continuous external criterion. Designed for item-level analysis with ordinal data.
#' Confidence intervals are computed via bootstrap.
#'
#' @param data.items A data frame of ordinal item responses (e.g., Likert-type items).
#' @param criteria A numeric vector (external continuous variable) with the same length as the number of rows in `data.items`.
#' @param ci Logical. Whether to compute confidence intervals. Default is TRUE.
#' @param B Integer. Number of bootstrap replicates. Default is 1000.
#' @param method Character. Correlation method: "spearman" (default) or "kendall".
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{Item}: Name of the item.
#'   \item \code{rho}: Correlation coefficient (Spearman's rho or Kendall's tau).
#'   \item \code{lwr.ci}: Lower bound of the confidence interval (if \code{ci = TRUE}).
#'   \item \code{lwr.upp}: Upper bound of the confidence interval (if \code{ci = TRUE}).
#' }
#'
#' @details
#' For ordinal data with many ties or small samples, Kendall's tau is recommended.
#' If an item shows high tie frequency (less than 50% unique values), a message will suggest using method = "kendall".
#'
#' @examples
#' # Simulated example
#' set.seed(123)
#' items <- data.frame(
#'   Item1 = sample(1:5, 100, replace = TRUE),
#'   Item2 = sample(1:5, 100, replace = TRUE)
#' )
#' edad <- rnorm(100, mean = 30, sd = 10)
#'
#' # Using Spearman's rho (default)
#' SpearItems(items, edad, ci = TRUE, B = 500)
#'
#' # Using Kendall's tau (more conservative, robust to ties)
#' SpearItems(items, edad, ci = TRUE, B = 500, method = "kendall")
#'
#' @export
SpearItems <- function(data.items, criteria, ci = TRUE, B = 1000, method = "spearman") {
  requireNamespace(rcompanion)

  output <- lapply(names(data.items), function(x) {
    vec <- data.items[[x]]

    if (method == "spearman" && length(unique(vec)) / length(vec) < 0.5) {
      message(paste0("⚠️  Item '", x, "' has a high number of ties. Consider using method = 'kendall'."))
    }

    res <- spearmanRho(x = vec, y = criteria,
                       method = method,
                       ci = ci, type = "perc", R = B)

    if (ci) {
      data.frame(Item = x,
                 rho = res[1],
                 lwr.ci = res[2],
                 lwr.upp = res[3])
    } else {
      data.frame(Item = x,
                 rho = as.numeric(res))
    }
  })

  result <- do.call(rbind, output)
  rownames(result) <- NULL
  return(result)
}
