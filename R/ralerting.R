#'  @title ralerting-Construct Validity with permutation test and bootstrap CIs
#'
#' @description
#' Computes the alerting index (ralerting-CV) for the Quantifying Construct Validity procedure (QCV) as the correlation between
#' Fisher-z transformed observed correlations and mean-centered predicted
#' correlations. Provides a one-sided randomization (permutation) p-value and
#' bootstrap confidence intervals for ralerting. Includes a minimal-variance
#' jitter mechanism to handle constant predicted profiles.
#'
#' @param actr Numeric vector. Observed validity correlations (one per criterion).
#' @param predr Numeric vector. Predicted validity correlations (same length as \code{actr}).
#' @param iter Integer. Number of permutations for the randomization p-value
#'   (default \code{1000}). Use \code{0} to skip permutation testing.
#' @param B Integer. Number of bootstrap replicates for the CI (default \code{1000}).
#'   Use \code{0} to skip the CI.
#' @param conf.level Numeric in (0,1). Confidence level for the CI (default \code{0.95}).
#' @param type Character. Bootstrap CI type: \code{"perc"} (default), \code{"basic"}, or \code{"norm"}.
#' @param jitter Numeric scalar. Amplitude of uniform noise added to \code{predr} if it is constant.
#'   Interpreted as \eqn{\pm}jitter. Set to \code{0} to disable (default \code{0.001}).
#' @param seed Integer or \code{NULL}. Seed for reproducibility (optional).
#' @param digits Integer. Number of decimals to round the reported statistics (default \code{3}).
#'
#' @details
#' This function implements the QCV alerting index (ralerting-CV) following
#' the algorithm described by Westen & Rosenthal (2003) and clarified by
#' Furr & Heuckeroth (2019), and mirrors the core logic of \code{qcv::ralertingp()}:
#' (a) Fisher z-transform the observed validity correlations, and
#' (b) mean-center the predicted correlations to obtain contrast weights
#'     (\eqn{\lambda = r_{pred} - \bar{r}_{pred}}), then
#' (c) compute \eqn{r_{alerting} = cor(Zr, \lambda)}.
#'
#' The present implementation extends \code{qcv::ralertingp()} in two ways:
#' \enumerate{
#'   \item It offers a minimal-variance jitter mechanism for constant predicted
#'         profiles (default \eqn{\pm 0.001}) to enable computation while
#'         transparently reporting the perturbation used.
#'   \item It provides nonparametric bootstrap confidence intervals (percentile,
#'         basic, or normal) for \code{ralerting}, based on paired resampling of
#'         \code{(actr_i, predr_i)}.
#' }
#' The one-sided randomization p-value is computed by permuting \eqn{\lambda}
#' (as in \code{qcv::ralertingp()}), estimating \eqn{p = mean(r_{null} \ge r_{obs})}.
#'
#' \strong{Interpretation.} \code{ralerting} quantifies shape similarity between the
#' hypothesized profile and the observed validity profile (Furr & Heuckeroth, 2019).
#' Large values indicate that criteria predicted to correlate more strongly with the
#' focal test indeed do so, and vice versa.
#'
#' \strong{Caveat on constant predictions.} When \code{predr} is constant (e.g., using a
#' single reference value, as in some Dm-index applications), ralerting-CV is
#' undefined without variance in \eqn{\lambda}. The jitter option is a pragmatic
#' device to permit computation; use with interpretive caution.
#'
#' @seealso \code{qcv::ralertingp()} for the original randomization test of ralerting-CV,
#'   and \code{qcv} package vignettes/examples for QCV workflows.
#'
#' @references
#' Furr, R. M., & Heuckeroth, S. (2019). The “quantifying construct validity” procedure: Its role,
#' value, interpretations, and computation. \emph{Assessment}, 26(4), 555–566. https://doi.org/10.1177/1073191118820638
#'
#'Westen, D., & Rosenthal, R. (2003). Quantifying construct validity: Two simple measures.
#'\emph{Journal of Personality and Social Psychology}, 84(3), 608–618. https://doi.org/10.1037/0022-3514.84.3.608
#'
#' See also the \code{qcv} R package for an implementation of QCV indices.
#'
#' @examples
#' # Example 1 -----------
#'
#' # Constant predictions (e.g., ref value from DmIndex): jitter enables computation
#' ralerting(actr = c(.233, .449, .257, .092, -.002, .174),
#'            predr = rep(.50, 6),
#'            iter = 1000, B = 1000, seed = 123)
#'
#' # Example 2 -----------
#' # Strict mode (no jitter) -> NA if predictions are constant
#' ralerting(actr = c(.233, .449, .257, .092, -.002, .174),
#'            predr = rep(.50, 6),
#'            jitter = 0, iter = 2000, B = 0, seed = 123)
#'
#' # Example 3 -----------
#' # Non-constant predictions: no jitter applied
#' ralerting(actr = c(.10, .35, .08, .31, .17, .33),
#'            predr = c(.19, .201, .21, .201, .20, .201),
#'            iter = 2000, B = 2000, seed = 321)
#'
#' @export
ralerting <- function(actr, predr, iter = 1000, B = 1000, conf.level = 0.95,
                       type = c("perc","basic","norm"), jitter = 0.001,
                       seed = NULL, digits = 3) {
  type <- match.arg(type)
  stopifnot(is.numeric(actr), is.numeric(predr), length(actr) == length(predr))
  n <- length(actr)

  # paired NA removal
  ok <- is.finite(actr) & is.finite(predr)
  if (!all(ok)) {
    actr <- actr[ok]; predr <- predr[ok]; n <- length(actr)
    if (n < 3) stop("Insufficient number of criteria after removing NAs.")
  }

  # Fisher z, clipping to avoid +/-Inf
  eps <- 1e-7
  actr <- pmin(pmax(actr, -1 + eps), 1 - eps)
  zactr <- atanh(actr)

  # centered predictions
  lambda <- predr - mean(predr)
  const.pred <- (stats::sd(lambda) == 0)

  jitter.used <- 0
  note <- ""

  if (const.pred) {
    if (jitter <= 0) {
      warning("ralerting-CV not computable: 'predr' is constant and jitter = 0.")
      return(data.frame(
        ralerting = NA_real_, p.perm = NA_real_,
        ci.lwr = NA_real_, ci.upr = NA_real_,
        method.ci = type, conf.level = conf.level,
        n.criteria = n, iter = iter, B = B,
        jitter.used = 0, note = "predr constant; no jitter; ralerting undefined."
      ))
    } else {
      if (!is.null(seed)) set.seed(seed)
      predr <- predr + stats::runif(n, min = -jitter, max = jitter)
      lambda <- predr - mean(predr)
      jitter.used <- jitter
      note <- sprintf("predr was constant; added jitter +/-%.4f to allow computation.", jitter)
    }
  }

  if (stats::sd(lambda) == 0) {
    warning("ralerting-CV not computable: lambda has zero variance after jitter.")
    return(data.frame(
      ralerting = NA_real_, p.perm = NA_real_,
      ci.lwr = NA_real_, ci.upr = NA_real_,
      method.ci = type, conf.level = conf.level,
      n.criteria = n, iter = iter, B = B,
      jitter.used = jitter.used, note = "lambda no variance after jitter."
    ))
  }

  # observed ralerting
  r.obs <- as.numeric(stats::cor(zactr, lambda))

  # permutation p-value (one-sided)
  p.perm <- NA_real_
  if (iter > 0) {
    if (!is.null(seed)) set.seed(seed + 1L)
    rnull <- replicate(iter, {
      lam.shuf <- sample(lambda, replace = FALSE)
      as.numeric(stats::cor(zactr, lam.shuf))
    })
    p.perm <- mean(rnull >= r.obs)
  }

  # bootstrap CI
  ci.lwr <- ci.upr <- NA_real_
  if (B > 0) {
    if (!is.null(seed)) set.seed(seed + 2L)
    r.boot <- numeric(B)
    idx_all <- seq_len(n)
    for (b in seq_len(B)) {
      idx  <- sample(idx_all, size = n, replace = TRUE)
      z_b  <- zactr[idx]
      pr_b <- predr[idx]  # use (possibly jittered) predictions, paired
      lam_b <- pr_b - mean(pr_b)
      if (stats::sd(lam_b) == 0 && jitter > 0) {
        pr_b <- pr_b + stats::runif(n, min = -jitter, max = jitter)
        lam_b <- pr_b - mean(pr_b)
      }
      r.boot[b] <- if (stats::sd(lam_b) == 0) NA_real_ else as.numeric(stats::cor(z_b, lam_b))
    }
    r.boot <- r.boot[is.finite(r.boot)]
    if (length(r.boot) < max(10, 0.25*B)) {
      warning("Too many NA bootstrap replicates; CI not computed.")
    } else {
      alpha <- 1 - conf.level
      if (type == "perc") {
        q <- stats::quantile(r.boot, probs = c(alpha/2, 1 - alpha/2), names = FALSE, na.rm = TRUE)
        ci.lwr <- q[1]; ci.upr <- q[2]
      } else if (type == "basic") {
        q <- stats::quantile(r.boot, probs = c(alpha/2, 1 - alpha/2), names = FALSE, na.rm = TRUE)
        ci.lwr <- 2*r.obs - q[2]; ci.upr <- 2*r.obs - q[1]
      } else if (type == "norm") {
        m <- mean(r.boot, na.rm = TRUE); s <- stats::sd(r.boot, na.rm = TRUE)
        z <- stats::qnorm(c(alpha/2, 1 - alpha/2))
        ci.lwr <- r.obs + z[1]*s; ci.upr <- r.obs + z[2]*s
      }
    }
  }

  data.frame(
    ralerting = round(r.obs, digits),
    p.perm     = ifelse(is.na(p.perm), NA, round(p.perm, digits)),
    ci.lwr     = ifelse(is.na(ci.lwr), NA, round(ci.lwr, digits)),
    ci.upr     = ifelse(is.na(ci.upr), NA, round(ci.upr, digits)),
    method.ci  = if (B > 0) type else "",
    conf.level = if (B > 0) conf.level else NA_real_,
    n.criteria = n,
    iter       = iter,
    B          = B,
    jitter.used = jitter.used,
    note       = note,
    row.names  = NULL
  )
}
