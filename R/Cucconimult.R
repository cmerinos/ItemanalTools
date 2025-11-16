#' @title Multigroup Cucconi Test for Location and Scale Differences
#'
#' @description
#' Performs a nonparametric test for evaluating whether multiple independent groups
#' differ in both location (central tendency) and scale (dispersion). This function
#' generalizes the original Cucconi test (1968) to multiple groups using permutation.
#'
#' @param x A numeric vector of values (pooled sample).
#' @param group A factor or grouping vector of the same length as \code{x}, indicating group membership.
#' @param B Integer. Number of permutations used to compute empirical p-values. Default is 1000.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{\code{Partial}}{One row per group showing each group's partial Cucconi statistic.}
#'   \item{\code{Global}}{One row with the overall Cucconi statistic and permutation-based p-value.}
#' }
#'
#' @details
#' The test statistic is an average of squared standardized rank contrasts per group.
#' The null distribution is approximated via permutation of the pooled data.
#' Based on the method proposed by Marozzi (2014).
#'
#' @references
#' Cucconi, O. (1968). Un nuovo test non parametrico per il confronto tra due gruppi campionari.
#' \emph{Giornale degli Economisti e Annali di Economia}, 17(1), 225-–248. https://www.jstor.org/stable/23241361
#'
#' Marozzi, M. (2014). The multisample Cucconi test.
#' \emph{Statistical Methods and Applications}, 23, 209--227. https://doi.org/10.1007/s10260-014-0255-x
#'
#' Marozzi, M. (2009). Some notes on the location–scale Cucconi test. \emph{Journal of Nonparametric Statistics},
#' 21(5), 629--647. https://doi.org/10.1080/10485250902952435
#'
#' Marozzi, M. (2013). Nonparametric Simultaneous Tests for Location and Scale Testing: A Comparison of Several Methods.
#' \emph{Communications in Statistics -- Simulation and Computation}, 42(6), 1298–1317. https://doi.org/10.1080/03610918.2012.665546
#'
#' @examples
#' set.seed(123)
#' x <- c(rnorm(10), rnorm(12, mean = 1), rnorm(8, sd = 2))
#' group <- factor(rep(c("A", "B", "C"), times = c(10, 12, 8)))
#' CucconiMult(x, group, B = 999)
#'
#' @export
CucconiMult <- function(x, group, B = 1000) {
  if (length(x) != length(group)) stop("'x' and 'group' must be the same length.")
  if (!is.numeric(x)) stop("'x' must be numeric.")

  group <- as.factor(group)
  group.levels <- levels(group)
  sample.sizes <- table(group)
  K <- length(sample.sizes)

  if (K < 2) stop("At least two groups are required.")

  pooled <- x[order(group)]
  sizes.vec <- as.vector(sample.sizes)

  computeStat <- function(pooled.sample, sample.sizes) {
    sample.sizes <- c(0, sample.sizes)
    n <- sum(sample.sizes)
    K <- length(sample.sizes)

    ranks <- rank(pooled.sample)
    contr.ranks <- n + 1 - ranks

    group.ranks <- vector("list", K - 1)
    group.contr.ranks <- vector("list", K - 1)
    cum.sizes <- cumsum(sample.sizes)

    for (k in 1:(K - 1)) {
      group.ranks[[k]] <- ranks[(cum.sizes[k] + 1):(cum.sizes[k + 1])]
      group.contr.ranks[[k]] <- contr.ranks[(cum.sizes[k] + 1):(cum.sizes[k + 1])]
    }

    sample.sizes <- sample.sizes[-1]

    means <- sample.sizes * (n + 1) * (2 * n + 1) / 6
    sds <- sqrt(sample.sizes * (n - sample.sizes) * (n + 1) * (2 * n + 1) * (8 * n + 11) / 180)
    covariance <- -(30 * n + 14 * n^2 + 19) / ((8 * n + 11) * (2 * n + 1))

    u <- sapply(1:(K - 1), function(k) {
      (sum(group.ranks[[k]]^2) - means[k]) / sds[k]
    })

    v <- sapply(1:(K - 1), function(k) {
      (sum(group.contr.ranks[[k]]^2) - means[k]) / sds[k]
    })

    partial.c <- (u^2 + v^2 - 2 * u * v * covariance) / (2 * (1 - covariance^2))
    c.stat <- mean(partial.c)

    return(list(c.stat = c.stat, partial = partial.c))
  }

  # Permutations
  perm.stats <- replicate(B, {
    shuffled <- sample(pooled)
    computeStat(shuffled, sizes.vec)$c.stat
  })

  observed <- computeStat(pooled, sizes.vec)
  pval <- mean(perm.stats >= observed$c.stat)

  partial.df <- data.frame(
    Group = group.levels,
    PartialStatistic = round(observed$partial, 4)
  )

  global.df <- data.frame(
    TestStatistic = round(observed$c.stat, 4),
    p.value = round(pval, 4)
  )

  return(list(Partial = partial.df, Global = global.df))
}
