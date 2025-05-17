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
#' @return A data frame with one row per item. If \code{ci = TRUE}, columns are:
#' \describe{
#'   \item{\code{Item}}{Item name.}
#'   \item{\code{rg}}{Wilcoxon rank-biserial correlation.}
#'   \item{\code{lwr.ci}}{Lower confidence limit.}
#'   \item{\code{lwr.upp}}{Upper confidence limit.}
#' }
#' If \code{ci = FALSE}, only \code{Item} and \code{rg} are returned.
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
        return(data.frame(Item = x, rg = NA, lwr.ci = NA, lwr.upp = NA))
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
                        lwr.upp = round(upper, 3)))
    } else {
      rg <- if ("rg" %in% colnames(res)) res[1, "rg"] else NA
      return(data.frame(Item = x, rg = round(rg, 3)))
    }
  })

  result <- do.call(rbind, output)
  rownames(result) <- NULL
  return(result)
}
