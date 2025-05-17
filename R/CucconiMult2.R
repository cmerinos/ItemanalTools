#' @title CucconiMult2: Multigroup Cucconi Test for Multiple Variables
#'
#' @description
#' Applies a multisample Cucconi test across multiple variables (e.g., test items),
#' evaluating simultaneous location and scale differences among groups for each variable.
#'
#' @param data A data.frame or matrix. Each column represents a variable or item.
#' @param group A factor or grouping vector (same length as number of rows in `data`).
#' @param B Integer. Number of permutations (default = 1000).
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{\code{Global}}{One row per variable with test statistic and p-value.}
#'   \item{\code{Partial}}{One row per variable × group with partial statistic.}
#' }
#'
#' @export
CucconiMult2 <- function(data, group, B = 1000) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (length(group) != nrow(data)) stop("Length of group must match number of rows in data")
  if (!is.factor(group)) group <- as.factor(group)
  
  group_levels <- levels(group)
  results.global <- list()
  
  # Preparar lista vacía para cada grupo
  results.partial.bygroup <- setNames(vector("list", length(group_levels)), group_levels)
  for (g in group_levels) {
    results.partial.bygroup[[g]] <- data.frame()
  }
  
  for (var in colnames(data)) {
    x <- data[[var]]
    pooled.sample <- x[order(group)]
    sample.sizes <- table(group)
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
    
    results.global[[var]] <- data.frame(
      Variable = var,
      TestStatistic = round(observed$c.stat, 4),
      p.value = round(pval, 4)
    )
    
    # Registrar parciales por grupo
    for (i in seq_along(group_levels)) {
      g <- group_levels[i]
      results.partial.bygroup[[g]] <- rbind(
        results.partial.bygroup[[g]],
        data.frame(Variable = var, PartialStatistic = round(observed$partial[i], 4))
      )
    }
  }
  
  Global <- do.call(rbind, results.global)
  Partial <- results.partial.bygroup
  
  return(list(Global = Global, Partial = Partial))
}
