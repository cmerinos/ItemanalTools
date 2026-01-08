#' @title Actual Equivalent Number of Options for one item
#'
#' @description
#' Calculates the effective number of response categories (or actual equivalent number of
#' options) used in a polytomous item, based on Shannon entropy (Samejima, 1980;
#' Sato & Morimoto, 1976). The result reflects the dispersion or actual utilization
#' of the available response options.
#'
#' Response categories with zero frequency are treated as *incidental zeros*
#' (Linacre, 1999; 2002). To avoid issues with undefined logarithms, such frequencies
#' are replaced with a small value defined by the argument \code{pseudo_zero}.
#' Within a content validity framework (Merino-Soto et al., 2022), this indicator provides
#' objective information to assess the functionality of the response options.
#'
#' @param vector A numeric, factor, or character vector of item responses.
#' @param k The total number of theoretical response options (e.g., \code{k = 6} for 0–5).
#'          This value is required, even if some categories were unused.
#' @param min.cat Possible minimum value of theoretical response options. Default is \code{1}.
#' @param pseudo_zero A small positive value used to replace zero frequencies,
#'        avoiding undefined logarithms. Default is \code{0.05}.
#'
#' @return A list with one element:
#' \item{AENO}{A numeric value indicating the effective number of options used.}
#'
#' @examples
#' responses <- c(0, 1, 1, 2, 3, 3, 5, NA)
#' AENOuni(responses, k = 6, min.cat = 0)
#'
#' @references
#' Sato, T., & Morimoto, U. (1976). Sentaku-shi keishiki tesuto kaitou bunpu no bunseki
#' (Analyzing endorsement distribution of selected-response items).
#' In *Proceedings of the 4th Annual Meeting of the Behaviometric Society of Japan*, Tokyo, Japan.
#'
#' Samejima, F. (1980). \emph{Research on the Multiple-Choice Test Item in Japan: Toward the Validation of Mathematical Models}.
#' Scientific Monograph, ONRT-1M3; Department of the Navy Office of Naval Research: Tokyo, Japan.
#'
#' Merino-Soto, C., Juárez-García, A., Salinas-Escudero, G., & Toledano-Toledano, F. (2022).
#' Item-Level Psychometric Analysis of the Psychosocial Processes at Work Scale (PROPSIT) in Workers.
#' \emph{International Journal of Environmental Research and Public Health}, 19(13), 7972. https://doi.org/10.3390/ijerph19137972
#'
#' Linacre, J. M. (1999). Investigating rating scale category utility. \emph{Journal of Outcome Measurement}, 3(2), 103–122.
#'
#' Linacre, J. M. (2002). Optimizing rating scale category effectiveness. \emph{Journal of Applied Measurement}, 3(1), 85–106.
#'
#' @export
AENOuni <- function(vector, k, min.cat = 1, pseudo_zero = 0.05) {

  if (missing(k)) stop("You must specify 'k' (total number of response options).")
  if (!min.cat %in% c(0, 1)) stop("'min.cat' must be 0 or 1.")

  # observed proportions (omit NA automatically)
  prop_responses <- table(vector) / sum(!is.na(vector))

  # define theoretical categories
  full_categories <- as.character(min.cat:(min.cat + k - 1))

  # IMPORTANT: prevent vector extension by matching
  prop_full <- setNames(rep(0, k), full_categories)
  common <- intersect(names(prop_responses), full_categories)
  prop_full[common] <- prop_responses[common]

  # if there are observed categories not in the theoretical set, stop
  extra <- setdiff(names(prop_responses), full_categories)
  if (length(extra) > 0) {
    stop(
      "Observed categories not covered by (min.cat:(min.cat+k-1)): ",
      paste(extra, collapse = ", "),
      ". Please set 'min.cat' correctly or recode the item."
    )
  }

  # incidental zeros
  prop_full[prop_full == 0] <- pseudo_zero
  prop_full <- prop_full / sum(prop_full)

  entropy <- -sum(prop_full * log2(prop_full))
  AENO <- 2^entropy

  return(list(AENO = AENO))
}
