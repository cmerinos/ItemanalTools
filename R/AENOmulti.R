#' @title AENOmulti: Effective Number of Nominal Options for Multiple Items
#'
#' @description
#' Calculates the AENO index (Effective Number of Nominal Options; Samejima, 1980) for each item
#' in a dataset. Useful for exploring how well response options are being utilized
#' across polytomous items. AENO is a summary practical index for content validity
#' focused in item analysis (Merino-Soto et al., 2022)
#'
#' Internally calls \code{\link{AENOuni}}. Response categories with zero frequency
#' are treated as *incidental zeros* (Linacre, 1999; 2002) and replaced with a small
#' value (default \code{0.05}) to allow entropy-based computations.
#'
#' @param data A \code{data.frame} of item responses, where each column is an item.
#' @param k Integer. Total number of response options possible (e.g., \code{k = 6} if options range from 0 to 5).
#'
#' @return A \code{data.frame} with the following columns for each item:
#' \item{item}{Item name (column name from input).}
#' \item{AENO}{Effective number of response options used (rounded to 3 decimals).}
#' \item{cat.eff}{Number of observed categories actually used.}
#' \item{n}{Number of valid (non-missing) responses.}
#'
#' @examples
#' data <- data.frame(
#'   item1 = c(0, 1, 2, 3, 4, 5, NA),
#'   item2 = c(1, 2, 2, 2, 3, 4, 5)
#' )
#' AENOmulti(data, k = 6)
#'
#' @seealso \code{\link{AENOuni}}
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
AENOmulti <- function(data, k) {

  if (!is.data.frame(data)) {
    stop("Input must be a data.frame with items as columns.")
  }
  if (missing(k)) {
    stop("You must specify the total number of response options with the argument 'k'.")
  }

  results <- lapply(names(data), function(item_name) {
    vector <- data[[item_name]]
    vector_no_na <- vector[!is.na(vector)]
    cat_eff <- length(unique(vector_no_na))
    n_valid <- length(vector_no_na)

    AENO <- round(AENOuni(vector, k = k)$AENO, 3)

    data.frame(
      item = item_name,
      AENO = AENO,
      cat.eff = cat_eff,
      n = n_valid,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results)
  return(result_df)
}
