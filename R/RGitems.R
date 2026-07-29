#' @title RGitems: Wilcoxon Rank-Biserial Effect Size for Multiple Items
#'
#' @description
#' Computes the rank-biserial correlation coefficient (\code{r_g}) for multiple items
#' using the Wilcoxon-Mann-Whitney test between two independent groups.
#' Optionally returns confidence intervals computed via bootstrap.
#'
#' @param data.items A data frame or matrix where each column is an item and each row is a subject.
#' @param group A grouping variable (factor or character vector) indicating group membership
#'              for each subject. Must contain exactly two unique levels.
#' @param ci Logical. Should confidence intervals be computed? Default is \code{TRUE}.
#' @param type Character. Type of bootstrap confidence interval. One of \code{"norm"},
#' \code{"basic"}, \code{"perc"}, or \code{"bca"}. Default is \code{"norm"}.
#' @param B Integer. Number of bootstrap replicates for interval estimation. Default is \code{1000}.
#' @param conf.level Confidence level for the interval estimation. Default is \code{0.95}.
#'
#'@details
#' This function internally calls \code{rcompanion::wilcoxonRG} to compute the
#' rank-biserial effect size (RG). The wrapper standardizes the output for
#' item-level psychometric analysis.
#'
#' @return A data frame with one row per item. If \code{ci = TRUE}, columns are:
#' \describe{
#'   \item{\code{Item}}{Item name.}
#'   \item{\code{rg}}{Wilcoxon rank-biserial correlation.}
#'   \item{\code{lwr.ci}}{Lower confidence limit.}
#'   \item{\code{upr.ci}}{Upper confidence limit.}
#' }
#' If \code{ci = FALSE}, only \code{Item} and \code{rg} are returned.
#'
#'@references
#'
#'Willson, V. L. (1976). Critical Values of the Rank-Biserial Correlation Coefficient. \emph{Educational and Psychological Measurement}, 36(2), 297-300. \doi{10.1177/001316447603600207}
#'
#'Cureton, E. E. (1968). Rank-Biserial Correlation when Ties Are Present. \emph{Educational and Psychological Measurement}, 28(1), 77-79. \doi{10.1177/001316446802800107}
#'
#'Glass, G. V. (1966). Note on Rank Biserial Correlation. \emph{Educational and Psychological Measurement}, 26(3), 623-631. \doi{10.1177/001316446602600307}
#'
#'Cureton, E. E. (1956). Rank-biserial correlation. \emph{Psychometrika}, 21, 287-290. \doi{10.1007/BF02289138}
#'
#' @examples
#' \dontrun{
#' data(iris)
#' sub <- subset(iris, Species %in% c("setosa", "versicolor"))
#' RGitems(data.items = sub[, 1:4], group = sub$Species)
#' }
#'
#' @importFrom rcompanion wilcoxonRG
#'
#' @export
RGitems <- function(data.items, group,
                    ci = TRUE,
                    type = "norm",
                    B = 1000,
                    conf.level = 0.95) {

  if (!requireNamespace("rcompanion", quietly = TRUE)) {
    stop("Package 'rcompanion' is required.")
  }

  if (!is.data.frame(data.items) && !is.matrix(data.items)) {
    stop("'data.items' must be a data frame or matrix.")
  }

  data.items <- as.data.frame(data.items)

  if (length(group) != nrow(data.items)) {
    stop("Length of 'group' must match the number of rows in 'data.items'.")
  }

  if (length(unique(group)) != 2) {
    stop("The 'group' variable must have exactly two levels.")
  }

  output <- lapply(names(data.items), function(x) {
    res <- tryCatch({
      rcompanion::wilcoxonRG(x = data.items[[x]],
                             g = group,
                             ci = ci,
                             type = type,
                             R = B,
                             conf = conf.level)
    }, error = function(e) return(NULL))

    if (is.null(res)) {
      if (ci) {
        return(data.frame(Item = x, rg = NA, lwr.ci = NA, upr.ci = NA))
      } else {
        return(data.frame(Item = x, rg = NA))
      }
    }

    if (ci) {
      rg <- res[1, "rg"]
      lower <- if ("lower.ci" %in% colnames(res)) res[1, "lower.ci"] else NA
      upper <- if ("upper.ci" %in% colnames(res)) res[1, "upper.ci"] else NA
      return(data.frame(Item = x,
                        rg = round(rg, 3),
                        lwr.ci = round(lower, 3),
                        upr.ci = round(upper, 3)))
    } else {
      # Cuando ci = FALSE, res es un vector nombrado
      rg <- if ("rg" %in% names(res)) res["rg"] else NA
      return(data.frame(Item = x, rg = round(rg, 3)))
    }

  })

  result <- do.call(rbind, output)
  rownames(result) <- NULL
  return(result)
}
