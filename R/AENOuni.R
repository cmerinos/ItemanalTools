#' @title AENOuni: Actual Equivalent Number of Options for one item
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
#' @param pseudo_zero A small positive value used to replace zero frequencies,
#'        avoiding undefined logarithms. Default is \code{0.05}.
#'
#' @return A list with one element:
#' \item{AENO}{A numeric value indicating the effective number of options used.}
#'
#' @examples
#' responses <- c(0, 1, 1, 2, 3, 3, 5, NA)
#' AENOuni(responses, k = 6)
#'
#' @references
#' Sato, T., & Morimoto, U. (1976). Sentaku-shi keishiki tesuto kaitou bunpu no bunseki
#' (Analyzing endorsement distribution of selected-response items).
#' In *Proceedings of the 4th Annual Meeting of the Behaviometric Society of Japan*, Tokyo, Japan.
#'
#' Samejima, F. (1980). *Research on the Multiple-Choice Test Item in Japan: Toward the Validation of Mathematical Models*.
#' Scientific Monograph, ONRT-1M3; Department of the Navy Office of Naval Research: Tokyo, Japan.
#'
#' Merino-Soto, C., Juárez-García, A., Salinas-Escudero, G., & Toledano-Toledano, F. (2022).
#' Item-Level Psychometric Analysis of the Psychosocial Processes at Work Scale (PROPSIT) in Workers.
#' *International Journal of Environmental Research and Public Health, 19*(13), 7972. https://doi.org/10.3390/ijerph19137972
#'
#' Linacre, J. M. (1999). Investigating rating scale category utility. *Journal of Outcome Measurement, 3*(2), 103–122.
#'
#' Linacre, J. M. (2002). Optimizing rating scale category effectiveness. *Journal of Applied Measurement, 3*(1), 85–106.
#'
#' @export
AENOuni <- function(vector, k, pseudo_zero = 0.05) {
  if (missing(k)) stop("You must specify the total number of response options with the 'k' argument.")

  prop_responses <- table(vector) / sum(!is.na(vector))

  full_categories <- as.character(0:(k - 1))
  prop_full <- rep(0, k)
  names(prop_full) <- full_categories
  prop_full[names(prop_responses)] <- prop_responses

  prop_full[prop_full == 0] <- pseudo_zero
  prop_full <- prop_full / sum(prop_full)

  entropy_terms <- prop_full * log2(prop_full)
  entropy <- -sum(entropy_terms)
  AENO <- 2^entropy

  return(list(AENO = AENO))
}
