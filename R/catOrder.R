#' @title catOrder: Category Order Diagnostics for Polytomous Items
#'
#' @description
#' Evaluates whether response categories of Likert-type items follow a monotonic trend
#' along a total score or latent proxy. For each item, it computes the mean total score
#' by response category and applies selected nonparametric trend tests (Cuzick, 1985;
#' Jonckheere, 1954).
#'
#' @param data A data frame containing polytomous item responses (e.g., values from 0 to 4 or 1 to 5).
#' Each column represents an item.
#' @param score.total A numeric vector containing the total or proxy score per respondent.
#' Typically the row sum of the item data.
#' @param tests Character vector indicating which trend tests to perform. Options are:
#' \code{"JT"} for Jonckheere-Terpstra test, and \code{"Cu"} for Cuzick's test.
#' Default is \code{c("JT", "Cu")}.
#'
#' @return A list with two main components:
#' \item{summary}{A data frame with one row per item, reporting item name, monotonicity flag,
#' number of violations, Z-statistics, p-values, and effect sizes for each selected test.}
#' \item{means.wide}{A wide-format data frame where each row is an item and each column
#' contains the mean total score for a given response category (e.g., \code{cat0.mean}, \code{cat1.mean}, ...).}
#'
#' An additional internal component, \code{"means.long"}, is stored as an attribute of the returned list.
#' This is not printed by default but can be accessed with \code{attr(result, "means.long")}, and is used
#' by functions like \code{\link{catPlot}}.
#'
#' @examples
#' \dontrun{
#' # Simulated 5-point Likert items
#' set.seed(123)
#' data <- data.frame(
#'   item1 = sample(0:4, 100, TRUE),
#'   item2 = sample(0:4, 100, TRUE),
#'   item3 = sample(0:4, 100, TRUE)
#' )
#' score <- rowSums(data)
#'
#' res <- catOrder(data, score.total = score, tests = "JT")
#'
#' res$summary       # summary table
#' res$means.wide    # category means in wide format
#' head(attr(res, "means.long"))  # internal long-format means
#' }
#'
#' @references
#' Cuzick J. (1985). A Wilcoxon-type test for trend. Statistics in medicine, 4(1),
#' 87–90. https://doi.org/10.1002/sim.4780040112
#'
#' Jonckheere, A. R. (1954) A distribution-free k-sample test against ordered alternatives.
#' Biometrica 41, 133–145. http://dx.doi.org/10.2307/2333011
#'
#' @importFrom PMCMRplus jonckheereTest cuzickTest
#'
#' @export
catOrder <- function(data, score.total, tests = c("JT", "Cu")) {
  if (!requireNamespace("PMCMRplus", quietly = TRUE)) stop("Package 'PMCMRplus' is required.")

  test.columns <- list(
    JT = c("JT.Z", "JT.p", "JT.ES"),
    Cu = c("Cu.Z", "Cu.p", "Cu.ES")
  )

  summary.list <- list()
  means.list <- list()
  means.wide <- list()

  for (item in names(data)) {
    responses <- data[[item]]
    valid <- !is.na(responses) & !is.na(score.total)
    x <- score.total[valid]
    g <- as.factor(responses[valid])
    N <- length(x)

    cat.means <- tapply(x, g, mean, na.rm = TRUE)
    is.mono <- all(diff(cat.means) >= 0)
    violations <- sum(diff(cat.means) < 0)

    cat.wide <- as.list(round(cat.means, 3))
    names(cat.wide) <- paste0("cat", names(cat.means), ".mean")
    cat.wide$Item <- item

    means.list[[item]] <- data.frame(
      Item = item,
      Category = as.numeric(names(cat.means)),
      Mean.Score = as.numeric(cat.means)
    )

    jt.z <- jt.p <- jt.es <- NA_real_
    cu.z <- cu.p <- cu.es <- NA_real_

    if ("JT" %in% tests) {
      jt <- suppressWarnings(PMCMRplus::jonckheereTest(x, g))
      jt.z <- as.numeric(jt$statistic)
      jt.p <- jt$p.value
      jt.es <- jt.z / sqrt(N)
    }

    if ("Cu" %in% tests) {
      cu <- suppressWarnings(PMCMRplus::cuzickTest(x, g))
      cu.z <- as.numeric(cu$statistic)
      cu.p <- cu$p.value
      cu.es <- cu.z / sqrt(N)
    }

    row <- data.frame(
      Item = item,
      Monotonic = ifelse(is.mono, "Yes", "No"),
      Violations = violations,
      JT.Z = jt.z,
      JT.p = jt.p,
      JT.ES = jt.es,
      Cu.Z = cu.z,
      Cu.p = cu.p,
      Cu.ES = cu.es,
      stringsAsFactors = FALSE
    )

    summary.list[[item]] <- row
    means.wide[[item]] <- cat.wide
  }

  summary.df <- do.call(rbind, summary.list)
  columns.to.keep <- c("Item", "Monotonic", "Violations", unlist(test.columns[tests]))
  summary.df <- summary.df[, columns.to.keep, drop = FALSE]

  means.wide.df <- do.call(rbind, lapply(means.wide, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  means.long.df <- do.call(rbind, means.list)

  out <- list(
    summary = summary.df,
    means.wide = means.wide.df
  )
  attr(out, "means.long") <- means.long.df

  return(out)
}
