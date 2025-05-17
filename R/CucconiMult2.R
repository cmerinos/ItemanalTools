#' @title CucconiMult2: Multivariable Cucconi Test for Multiple Groups
#'
#' @description
#' Performs a nonparametric test for assessing simultaneous location and scale differences
#' across multiple groups for several variables (e.g., test items).
#' It applies a multigroup Cucconi test to each variable using permutation.
#'
#' @param data A data.frame or matrix where each column represents a variable (e.g., item).
#' @param group A factor or grouping vector (same length as number of rows in \code{data}).
#' @param B Integer. Number of permutations to compute p-values. Default is 1000.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{Global}}{One row per variable, with the Cucconi test statistic and its p-value.}
#'   \item{\code{Partial}}{A named list with one data frame per group. Each contains
#'   partial statistics by variable.}
#' }
#'
#' @details
#' This function generalizes the original Cucconi test (Cucconi, 1968; Marozzi, 2014)
#' to multivariable designs. Each variable is tested independently.
#'
#' @references
#' Cucconi, O. (1968). Un nuovo test non parametrico per il confronto tra due gruppi campionari.
#' \emph{Giornale degli Economisti e Annali di Economia}, 17(1), 225–248.
#'
#' Marozzi, M. (2014). The multisample Cucconi test.
#' \emph{Statistical Methods and Applications}, 23, 209–227. \doi{10.1007/s10260-014-0255-x}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- data.frame(
#'   item1 = c(rnorm(10), rnorm(10, 1), rnorm(10)),
#'   item2 = c(rnorm(10), rnorm(10, 0), rnorm(10, sd = 2))
#' )
#' g <- factor(rep(c("A", "B", "C"), each = 10))
#' CucconiMult2(dat, g, B = 500)
#' }
#'
#' @export
CucconiMult2 <- function(data, group, B = 1000) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (length(group) != nrow(data)) stop("Length of 'group' must match number of rows in 'data'.")
  if (!is.factor(group)) group <- as.factor(group)

  group_levels <- levels(group)
  global_results <- list()
  partial_by_group <- setNames(vector("list", length(group_levels)), group_levels)
  for (g in group_levels) {
    partial_by_group[[g]] <- data.frame()
  }

  for (var in colnames(data)) {
    x <- data[[var]]
    pooled <- x[order(group)]
    sample_sizes <- table(group)
    sizes_vec <- as.vector(sample_sizes)

    computeStat <- function(pooled.sample, sample.sizes) {
      sample.sizes <- c(0, sample.sizes)
      n <- sum(sample.sizes)
      K <- length(sample.sizes)

      ranks <- rank(pooled.sample)
      contr.ranks <- n + 1 - ranks

      group_ranks <- vector("list", K - 1)
      group_contr <- vector("list", K - 1)
      cum.sizes <- cumsum(sample.sizes)

      for (k in 1:(K - 1)) {
        group_ranks[[k]] <- ranks[(cum.sizes[k] + 1):(cum.sizes[k + 1])]
        group_contr[[k]] <- contr.ranks[(cum.sizes[k] + 1):(cum.sizes[k + 1])]
      }

      sample.sizes <- sample.sizes[-1]

      means <- sample.sizes * (n + 1) * (2 * n + 1) / 6
      sds <- sqrt(sample.sizes * (n - sample.sizes) * (n + 1) * (2 * n + 1) * (8 * n + 11) / 180)
      covariance <- -(30 * n + 14 * n^2 + 19) / ((8 * n + 11) * (2 * n + 1))

      u <- sapply(1:(K - 1), function(k) (sum(group_ranks[[k]]^2) - means[k]) / sds[k])
      v <- sapply(1:(K - 1), function(k) (sum(group_contr[[k]]^2) - means[k]) / sds[k])

      partial <- (u^2 + v^2 - 2 * u * v * covariance) / (2 * (1 - covariance^2))
      c_stat <- mean(partial)
      return(list(c.stat = c_stat, partial = partial))
    }

    perm_stats <- replicate(B, {
      shuffled <- sample(pooled)
      computeStat(shuffled, sizes_vec)$c.stat
    })

    observed <- computeStat(pooled, sizes_vec)
    pval <- mean(perm_stats >= observed$c.stat)

    global_results[[var]] <- data.frame(
      Variable = var,
      TestStatistic = round(observed$c.stat, 4),
      p.value = round(pval, 4)
    )

    for (i in seq_along(group_levels)) {
      g <- group_levels[i]
      partial_by_group[[g]] <- rbind(
        partial_by_group[[g]],
        data.frame(Variable = var, PartialStatistic = round(observed$partial[i], 4))
      )
    }
  }

  Global <- do.call(rbind, global_results)
  Partial <- partial_by_group

  return(list(Global = Global, Partial = Partial))
}
