#' Guilford's average rank-order correlation for item sets
#'
#' Computes Guilford's average rank-order correlation coefficient for a set of items,
#' treating items as "judges" that rank respondents. The function returns the average
#' Spearman rank correlation between items, together with a Friedman-type test statistic
#' and, optionally, bootstrap confidence intervals.
#'
#' @param data A numeric matrix or \code{data.frame} with respondents in rows and items
#'   in columns. Typically, these are ordinal or Likert-type item scores. Only rows
#'   with complete data are used; rows containing \code{NA} are removed with a warning.
#' @param ties.correct Logical. If \code{TRUE} (default), the Friedman test is computed
#'   using \code{\link[stats]{friedman.test}} on the transposed data, which applies a
#'   standard tie correction to the test statistic. If \code{FALSE}, a simple
#'   Friedman-type statistic is computed without tie correction, based on ranked data.
#' @param ci Logical. If \code{TRUE}, bootstrap confidence intervals for the average
#'   correlation are computed. Default is \code{FALSE}.
#' @param conf.level Confidence level for the confidence interval when \code{ci = TRUE}.
#'   Default is \code{0.95}.
#' @param B Integer. Number of bootstrap resamples used when \code{ci = TRUE}.
#'   Default is \code{1000}.
#' @param ci.method Character string indicating the type of bootstrap confidence interval
#'   to compute when \code{ci = TRUE}. One of \code{"perc"}, \code{"basic"}, or
#'   \code{"bca"}. The default is \code{"perc"}. See \code{\link[boot]{boot.ci}}.
#' @param nd Integer. Number of decimal places used to round the numeric results.
#'   Default is \code{3}.
#'
#' @details
#' Guilford (1954) proposed an average rank-order correlation coefficient to quantify
#' the agreement among several judges who rank a set of stimuli. In the context of
#' item analysis, items can be viewed as judges that rank respondents. If item scores
#' are converted to ranks within each item, the pairwise Spearman correlations between
#' items reflect how similarly the items order the respondents.
#'
#' Let \eqn{k} denote the number of items (columns) and \eqn{N} the number of respondents
#' (rows). For each pair of items, a Spearman rank correlation is computed, and
#' \eqn{p_{\mathrm{avg}}} is defined as the average of all distinct pairwise Spearman
#' correlations among the \eqn{k} items. This average rank correlation is closely
#' related to Guilford's original coefficient, to Kendall's coefficient of concordance
#' \eqn{W}, and to the Friedman test statistic. Under suitable conditions, the
#' Friedman-type statistic can be approximated by a chi-square distribution with
#' \eqn{N - 1} degrees of freedom, which allows testing the null hypothesis of
#' zero average rank-order correlation among items (see Willson, 1980).
#'
#' Conceptually, \eqn{p_{\mathrm{avg}}} can be interpreted as a global measure of
#' how consistently a set of items orders respondents, in the same spirit in which
#' Kendall's \eqn{W} summarizes concordance among multiple judges. For ordinal or
#' Likert-type scales, the average Spearman rank correlation provides a robust,
#' nonparametric summary of item similarity. When \code{ci = TRUE}, the function
#' uses nonparametric bootstrap resampling of respondents to estimate confidence
#' intervals for \eqn{p_{\mathrm{avg}}}, which can be useful for inference and
#' comparison across scales or item sets.
#'
#' The function:
#' \itemize{
#'   \item Converts item scores to ranks internally via Spearman correlations.
#'   \item Computes the average pairwise Spearman correlation between items
#'         (\eqn{p_{\mathrm{avg}}}).
#'   \item Computes a Friedman-type statistic \eqn{S} and approximates its
#'         \emph{p}-value using a chi-square distribution with \eqn{N - 1} degrees
#'         of freedom.
#'   \item Optionally computes bootstrap confidence intervals for \eqn{p_{\mathrm{avg}}}
#'         using \pkg{boot}.
#' }
#'
#' This implementation is inspired by the relationship between Guilford's coefficient,
#' Friedman’s test, and Kendall’s \eqn{W} described in Guilford (1954), Hays (1981),
#' Howell (1997), and Willson (1980).
#'
#' @return
#' A \code{data.frame} with one row named \code{"pavg"} and the following columns:
#' \itemize{
#'   \item \code{estimate}: Average Spearman rank-order correlation between items.
#'   \item \code{S}: Friedman-type test statistic associated with the rank-order
#'         agreement among items.
#'   \item \code{df}: Degrees of freedom used for the chi-square approximation
#'         (equal to \eqn{N - 1}).
#'   \item \code{p.value}: Approximate \emph{p}-value for testing the null hypothesis
#'         of zero average rank-order correlation among items.
#'   \item \code{lwr.ci}, \code{upr.ci}: Lower and upper bounds of the bootstrap
#'         confidence interval for \code{estimate}, returned only when \code{ci = TRUE}.
#' }
#'
#' @references
#' Guilford, J. P. (1954). \emph{Psychometric methods} (2nd ed.). McGraw-Hill.
#'
#' Hays, W. L. (1981). \emph{Statistics} (3rd ed.). Holt, Rinehart and Winston.
#'
#' Howell, D. C. (1997). \emph{Statistical methods for psychology} (4th ed.).
#' Duxbury Press.
#'
#' Willson, V. L. (1980). The distribution of Guilford's average rank-order
#' correlation coefficient. \emph{Educational and Psychological Measurement, 40}(2), 427-431. \doi{10.1177/001316448004000209}
#'
#' @seealso
#' \code{\link[stats]{friedman.test}}, \code{\link[stats]{cor}}, \code{\link[boot]{boot}},
#' \code{\link[boot]{boot.ci}}.
#'
#' @examples
#' \dontrun{
#' ## Example with Likert-type item responses from the mirt package
#' if (requireNamespace("mirt", quietly = TRUE)) {
#'   data(Science, package = "mirt")
#'
#'   ## Average rank-order correlation and Friedman-type test
#'   Guilforditems(Science)
#'
#'   ## With bootstrap confidence interval
#'   Guilforditems(Science, ci = TRUE, B = 500, ci.method = "perc")
#' }
#'
#' ## Example with simulated ordinal data
#' set.seed(123)
#' dat <- matrix(sample(1:5, 200 * 5, replace = TRUE), ncol = 5)
#' colnames(dat) <- paste0("item", 1:5)
#'
#' Guilforditems(dat)
#' Guilforditems(dat, ci = TRUE, B = 200)
#' }
#'
#' @export
Guilforditems <- function(data,
                          ties.correct = TRUE,
                          ci = FALSE,
                          conf.level = 0.95,
                          B = 1000,
                          ci.method = c("perc", "basic", "bca"),
                          nd = 3) {

  ## --- Preparación de datos ---
  if (is.data.frame(data)) {
    X <- as.matrix(data)
  } else if (is.matrix(data)) {
    X <- data
  } else {
    stop("'data' must be a matrix or data.frame.")
  }

  if (!is.numeric(X)) {
    stop("'data' must contain numeric values (e.g., Likert scores).")
  }

  # Eliminar filas con NA (solo filas, no columnas)
  cc <- stats::complete.cases(X)
  removed.rows <- sum(!cc)

  if (removed.rows > 0) {
    X <- X[cc, , drop = FALSE]
    warning(sprintf("Removed %d rows containing missing values.", removed.rows))
  }

  N <- nrow(X)
  k <- ncol(X)

  if (N < 3L) stop("Need at least 3 subjects (rows) after removing missing data.")
  if (k < 2L) stop("Need at least 2 items (columns).")

  ## --- Coeficiente Guilford (pavg): promedio de correlaciones de Spearman ---
  R_s <- suppressWarnings(stats::cor(X, method = "spearman"))
  if (any(is.na(R_s))) {
    stop("Cannot compute Spearman correlations (check for constant columns/items).")
  }
  pavg <- mean(R_s[upper.tri(R_s)])

  ## --- Estadístico S tipo Friedman ---
  if (ties.correct) {
    ft <- stats::friedman.test(t(X))
    S  <- unname(ft$statistic)
    df <- unname(ft$parameter)  # debería ser N - 1
  } else {
    Rmat <- apply(X, 2, rank, ties.method = "average")
    if (is.vector(Rmat)) Rmat <- matrix(Rmat, ncol = 1L)
    Ri_dot <- rowSums(Rmat)
    S <- (12 / (k * N * (N + 1))) * sum(Ri_dot^2) - 3 * k * (N + 1)
    df <- N - 1L
  }

  # p-valor (aproximación chi-cuadrado)
  p.value <- stats::pchisq(S, df = df, lower.tail = FALSE)

  ## --- Intervalos de confianza por bootstrap (opcional) ---
  ci.method <- match.arg(ci.method)

  if (ci) {

    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for confidence intervals. Please install 'boot'.")
    }

    stat_fun <- function(dat, idx) {
      dd  <- dat[idx, , drop = FALSE]
      R_s <- suppressWarnings(stats::cor(dd, method = "spearman"))
      mean(R_s[upper.tri(R_s)])
    }

    boot_obj <- boot::boot(X, statistic = stat_fun, R = B)
    ci_obj   <- boot::boot.ci(boot_obj, conf = conf.level, type = ci.method)

    if (ci.method == "perc") {
      lwr <- ci_obj$percent[4]
      upr <- ci_obj$percent[5]
    } else if (ci.method == "basic") {
      lwr <- ci_obj$basic[4]
      upr <- ci_obj$basic[5]
    } else if (ci.method == "bca") {
      lwr <- ci_obj$bca[4]
      upr <- ci_obj$bca[5]
    }

    coefficients <- data.frame(
      estimate = pavg,
      S        = S,
      df       = df,
      p.value  = p.value,
      lwr.ci   = lwr,
      upr.ci   = upr
    )

  } else {

    coefficients <- data.frame(
      estimate = pavg,
      S        = S,
      df       = df,
      p.value  = p.value
    )
  }

  ## --- Redondeo ---
  coefficients <- as.data.frame(
    lapply(coefficients, function(x) if (is.numeric(x)) round(x, nd) else x)
  )
  rownames(coefficients) <- "pavg"

  return(coefficients)
}
