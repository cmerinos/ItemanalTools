#' @title MurakamiPairs: Pairwise Cucconi Test for Location and Scale Differences
#'
#' @description
#' Performs all-pairs nonparametric comparisons of groups using a Cucconi test,
#' which is sensitive to both location and scale differences.
#' The method follows the pairwise ranking approach recommended by Murakami (2016),
#' and p-values are computed via permutation for each group comparison.
#'
#' @param x A numeric vector containing the observed values.
#' @param group A factor or vector indicating group membership, with at least two levels.
#' @param B Integer. Number of permutations to compute the p-values (default = 1000).
#' @param adjust.method Character string indicating the method of p-value adjustment.
#' Options are those accepted by \code{\link{p.adjust}}: \code{"holm"}, \code{"bonferroni"},
#' \code{"BH"}, \code{"none"}, etc. Default is \code{"holm"}.
#'
#' @return A data frame with one row per pairwise group comparison, and the following columns:
#' \describe{
#'   \item{\code{Group1, Group2}}{The two groups being compared.}
#'   \item{\code{CucconiStatistic}}{The Cucconi statistic \( C_{st} \) computed for that pair.}
#'   \item{\code{p.value}}{Empirical p-value obtained via permutation.}
#'   \item{\code{p.adjusted}}{P-value adjusted for multiple comparisons using the specified method.}
#' }
#'
#' @details
#' This function implements the pairwise version of the Cucconi test based on the method
#' described in Murakami (2016). The test statistic is designed to simultaneously detect
#' differences in both location (central tendency) and scale (dispersion) between two groups.
#'
#' Each comparison is carried out using only the two groups involved (pairwise ranking),
#' and permutation resampling is used to obtain a p-value.
#'
#' The p-value returned in \code{p.adjusted} can be interpreted in the usual way:
#' a value below 0.05 suggests a significant difference between the two groups,
#' controlling for the familywise error rate or false discovery rate, depending on the
#' adjustment method used.
#'
#' The largest value of the Cucconi statistic among all pairs (i.e., \code{max(C_{st})})
#' is sometimes referred to as \code{SDC} in the literature, but should still be interpreted
#' as the result of a specific comparison, not a global test.
#'
#' @references
#' Murakami, H. (2016). All-pairs multiple comparisons based on the Cucconi test.
#' \emph{AStA Advances in Statistical Analysis}, 100, 355–368.
#' \doi{10.1007/s10182-016-0268-x}
#'
#' @examples
#' # Reproducing Example 1 from Murakami (2016)
#' g1 <- c(5.7, 10.2, 13.9, 10.3, 1.3, 12.0, 14.0,
#'         15.1, 8.8, 12.7, 15.0, 17.0, 15.2)
#' g2 <- c(8.3, 12.3, 6.1, 10.1, 6.3, 12.1, 13.0,
#'         13.4, 11.9, 9.9, 5.0, 6.0, 4.0)
#' g3 <- c(9.5, 8.1, 7.0, 7.8, 9.3, 12.2, 6.7,
#'         10.6, 6.6, 7.0, 8.0, 9.0, 3.0)
#'
#' x <- c(g1, g2, g3)
#' group <- factor(rep(1:3, each = 13))
#'
#' set.seed(123)
#' MurakamiPairs(x = x, group = group, B = 1000, adjust.method = "holm")
#'
#' @export
MurakamiPairs <- function(x, group, B = 1000, adjust.method = "holm") {
  if (!is.numeric(x)) stop("x must be numeric.")
  if (length(x) != length(group)) stop("x and group must be the same length.")
  
  group <- as.factor(group)
  groups <- levels(group)
  k <- length(groups)
  if (k < 2) stop("At least two groups are required.")
  
  results <- data.frame()
  
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      g1 <- groups[i]
      g2 <- groups[j]
      
      x1 <- x[group == g1]
      x2 <- x[group == g2]
      
      x_pair <- c(x1, x2)
      n1 <- length(x1)
      n2 <- length(x2)
      N <- n1 + n2
      
      ranks <- rank(x_pair)
      r1 <- ranks[1:n1]
      r1_contrary <- N + 1 - r1
      
      # Constantes
      mean_r2 <- n1 * (N + 1) * (2 * N + 1) / 6
      sd_r2 <- sqrt(n1 * n2 * (N + 1) * (2 * N + 1) * (8 * N + 11) / 180)
      rho <- 2 * (N^2 - 4) / ((2 * N + 1) * (8 * N + 11)) - 1
      
      # Estadísticos
      U <- (sum(r1^2) - mean_r2) / sd_r2
      V <- (sum(r1_contrary^2) - mean_r2) / sd_r2
      
      Cst_obs <- (U^2 + V^2 - 2 * rho * U * V) / (2 * (1 - rho^2))
      
      # Permutaciones
      C_perm <- numeric(B)
      for (b in 1:B) {
        idx_perm <- sample(N)
        x1_perm <- x_pair[idx_perm[1:n1]]
        x2_perm <- x_pair[idx_perm[(n1 + 1):N]]
        x_perm <- c(x1_perm, x2_perm)
        
        ranks_perm <- rank(x_perm)
        r1_perm <- ranks_perm[1:n1]
        r1c_perm <- N + 1 - r1_perm
        
        U_perm <- (sum(r1_perm^2) - mean_r2) / sd_r2
        V_perm <- (sum(r1c_perm^2) - mean_r2) / sd_r2
        
        C_perm[b] <- (U_perm^2 + V_perm^2 - 2 * rho * U_perm * V_perm) / (2 * (1 - rho^2))
      }
      
      p_val <- mean(C_perm >= Cst_obs)
      
      results <- rbind(results, data.frame(
        Group1 = g1,
        Group2 = g2,
        CucconiStatistic = round(Cst_obs, 4),
        p.value = round(p_val, 4)
      ))
    }
  }
  
  # Ajuste de p-valores
  results$p.adjusted <- round(p.adjust(results$p.value, method = adjust.method), 4)
  
  return(results)
}
