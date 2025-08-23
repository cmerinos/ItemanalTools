#' @title Item Thresholds for Likert-type Items
#'
#' @description 
#' Estimates thresholds (cutpoints) for an ordinal Likert-type item 
#' using a cumulative link model (CLM). 
#' Provides estimates, confidence intervals, differences between 
#' consecutive thresholds, and z-tests for those differences.
#'
#' @details 
#' This function fits an ordinal regression model without predictors 
#' (\code{ordinal::clm}) and extracts threshold estimates. 
#' The differences between consecutive thresholds are directly 
#' interpretable as standardized effect sizes on the latent variable scale, 
#' especially when the \code{link = "probit"} option is used 
#' (default). Larger differences indicate more distinct response 
#' categories, while small differences suggest redundancy.
#' 
#' The z-tests for threshold differences are computed as
#' \deqn{ z = \frac{ \tau_{i+1} - \tau_i }{ \sqrt{SE(\tau_{i+1})^2 + SE(\tau_i)^2} } }
#' where \eqn{\tau_i} is the i-th threshold and \eqn{SE} its standard error.
#' 
#' The function is designed as a lightweight, exploratory alternative 
#' to IRT or SEM approaches. For more exact threshold estimation, 
#' especially in large-scale applications, IRT models (e.g., GRM, PCM) 
#' or CFA with categorical items are recommended. 
#'
#' Thresholds are estimated using the \pkg{ordinal} package, which 
#' implements cumulative link models via maximum likelihood.
#'
#' @param x An ordered factor representing a Likert-type item.
#' @param link The link function to use in the CLM. Options include 
#' "logit", "probit" (default), "cloglog", etc.
#' @param conf.level Confidence level for the threshold intervals. 
#' Default is 0.95.
#'
#' @return A list with four data frames:
#' \item{thresholds}{Threshold estimates and their standard errors.}
#' \item{ci}{Confidence intervals for each threshold.}
#' \item{differences}{Differences between consecutive thresholds 
#' (interpretable as standardized effect sizes on the latent scale).}
#' \item{diff_tests}{z-tests for the differences between consecutive 
#' thresholds, including p-values.}
#'
#' @examples
#' set.seed(123)
#' y <- ordered(sample(1:5, 200, replace = TRUE,
#'                    prob = c(.1, .2, .3, .25, .15)))
#' res <- itemThresholds(y)
#' res$thresholds
#' res$ci
#' res$differences
#' res$diff_tests
#'
#' @seealso \code{\link[ordinal]{clm}}
#'
#' @references 
#' Christensen, R. H. B. (2019). *ordinal: Regression Models for 
#' Ordinal Data*. R package version 2019.12-10. 
#' https://CRAN.R-project.org/package=ordinal
#'
#' @importFrom stats coef pnorm qnorm
#' @export
itemThresholds <- function(x, link = "probit", conf.level = 0.95) {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("Package 'ordinal' is required. Please install it.")
  }
  if (!is.ordered(x)) stop("The argument 'x' must be an ordered factor.")
  
  # Fit CLM without predictors
  fit <- ordinal::clm(x ~ 1, link = link)
  
  # Extract thresholds and SE
  coefs <- coef(summary(fit))
  thr   <- coefs[grepl("\\|", rownames(coefs)), , drop = FALSE]
  est   <- thr[, "Estimate"]
  se    <- thr[, "Std. Error"]
  thr_names <- rownames(thr)
  
  thresholds <- data.frame(
    threshold = thr_names,
    estimate  = est,
    se        = se,
    stringsAsFactors = FALSE
  )
  
  # Confidence intervals
  zval <- qnorm(1 - (1 - conf.level) / 2)
  ci <- data.frame(
    threshold = thr_names,
    lwr.ci = est - zval * se,
    upr.ci = est + zval * se,
    stringsAsFactors = FALSE
  )
  
  # Differences and z-tests
  if (length(est) > 1) {
    diffs <- diff(est)
    diffs_df <- data.frame(
      contrast   = paste0(thr_names[-1], " - ", thr_names[-length(thr_names)]),
      difference = diffs,
      stringsAsFactors = FALSE
    )
    
    se_diffs <- sqrt(se[-1]^2 + se[-length(se)]^2)
    z_vals   <- diffs / se_diffs
    p_vals   <- 2 * (1 - pnorm(abs(z_vals)))
    
    diff_tests <- data.frame(
      contrast   = paste0(thr_names[-1], " - ", thr_names[-length(thr_names)]),
      difference = diffs,
      se.diff    = se_diffs,
      z          = z_vals,
      p.value    = p_vals,
      stringsAsFactors = FALSE
    )
  } else {
    diffs_df   <- data.frame()
    diff_tests <- data.frame()
  }
  
  return(list(
    thresholds  = thresholds,
    ci          = ci,
    differences = diffs_df,
    diff_tests  = diff_tests
  ))
}
