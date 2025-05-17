#' @title CucconiMult: Multigroup Cucconi Test for Location and Scale Differences
#'
#' @description
#' Performs a nonparametric test for assessing whether multiple independent groups
#' differ in both location (median) and scale (dispersion). This function generalizes
#' the original Cucconi test (1968) to multiple samples, using a permutation-based approach.
#'
#' @param x A numeric vector of observations (pooled sample).
#' @param group A factor or grouping vector of the same length as `x`, indicating group membership.
#' @param B Integer. Number of permutations to compute the empirical p-value. Default is 1000.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{\code{Partial}}{A data frame with one row per group showing the partial Cucconi statistic.}
#'   \item{\code{Global}}{A data frame with the overall test statistic and permutation-based p-value.}
#' }
#'
#' @details
#' This function implements a multisample version of the Cucconi test for location and scale
#' differences, using the permutation-based method proposed by Marozzi (2014).
#' The test statistic is an average of standardized squared rank-based contrasts for each group,
#' and the null distribution is approximated via permutation of the pooled sample.
#' This R implementation is a modified version based on the approach described in:
#' Marozzi (2014), \emph{The multisample Cucconi test}, see reference below.
#'
#' @references
#' Cucconi, O. (1968). Un nuovo test non parametrico per il confronto tra due gruppi campionari.
#' \emph{Giornale degli Economisti e Annali di Economia}, 17(1), 225–248.
#'
#' Marozzi, M. (2009). Some notes on the location-scale Cucconi test.
#' \emph{Journal of Nonparametric Statistics}, 21(5), 629–647.
#'
#' Marozzi, M. (2012). A modified Cucconi Test for Location and Scale Change Alternatives.
#' \emph{Revista Colombiana de Estadística}, 35(3), 371–384.
#'
#' Marozzi, M. (2013). Nonparametric simultaneous tests for location and scale testing: A comparison of several methods.
#' \emph{Communications in Statistics - Simulation and Computation}, 42(6), 1298–1317.
#'
#' Marozzi, M. (2014). The multisample Cucconi test.
#' \emph{Statistical Methods and Applications}, 23, 209–227. \doi{10.1007/s10260-014-0255-x}
#'
#' @examples
#' set.seed(123)
#' x <- c(rnorm(10), rnorm(12, mean = 1), rnorm(8, sd = 2))
#' group <- factor(rep(c("A", "B", "C"), times = c(10, 12, 8)))
#' CucconiMult(x, group, B = 999)
#'
#' @export
CucconiMult <- function(x, group, B = 1000) {
  if (length(x) != length(group)) stop("x and group must have the same length")
  if (!is.numeric(x)) stop("x must be numeric")
  
  group <- as.factor(group)
  sample.sizes <- table(group)
  group.levels <- levels(group)
  K <- length(sample.sizes)
  
  if (K < 2) stop("At least two groups are required")
  
  pooled.sample <- x[order(group)]
  sizes.vec <- as.vector(sample.sizes)
  
  MultiSampleCucconiStat <- function(pooled.sample, sample.sizes) {
    sample.sizes <- c(0, sample.sizes)
    n <- sum(sample.sizes)
    K <- length(sample.sizes)
    
    ranks <- rank(pooled.sample)
    contrary.ranks <- n + 1 - ranks
    
    groups.ranks <- vector("list", K - 1)
    groups.contrary.ranks <- vector("list", K - 1)
    cum.sample.sizes <- cumsum(sample.sizes)
    
    for (k in 1:(K - 1)) {
      groups.ranks[[k]] <- ranks[(cum.sample.sizes[k] + 1):(cum.sample.sizes[k + 1])]
      groups.contrary.ranks[[k]] <- contrary.ranks[(cum.sample.sizes[k] + 1):(cum.sample.sizes[k + 1])]
    }
    
    sample.sizes <- sample.sizes[-1]
    
    means <- sample.sizes * (n + 1) * (2 * n + 1) / 6
    st.dev <- sqrt(sample.sizes * (n - sample.sizes) * (n + 1) * (2 * n + 1) * (8 * n + 11) / 180)
    covariance <- -(30 * n + 14 * n^2 + 19) / ((8 * n + 11) * (2 * n + 1))
    
    u.stat <- sapply(1:(K - 1), function(k) {
      (sum(groups.ranks[[k]]^2) - means[k]) / st.dev[k]
    })
    
    v.stat <- sapply(1:(K - 1), function(k) {
      (sum(groups.contrary.ranks[[k]]^2) - means[k]) / st.dev[k]
    })
    
    partial.c.stat <- (u.stat^2 + v.stat^2 - 2 * u.stat * v.stat * covariance) /
      (2 * (1 - covariance^2))
    
    c.stat <- mean(partial.c.stat)
    
    return(list(c.stat = c.stat, partial = partial.c.stat))
  }
  
  permutation.stats <- replicate(B, {
    shuffled <- sample(pooled.sample)
    MultiSampleCucconiStat(shuffled, sizes.vec)$c.stat
  })
  
  observed <- MultiSampleCucconiStat(pooled.sample, sizes.vec)
  pval <- mean(permutation.stats >= observed$c.stat)
  
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
