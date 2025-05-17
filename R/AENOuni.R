#' AENOuni: Numero efectivo de opciones nominales para un ítem
#'
#' Calcula el numero efectivo de categorías utilizadas en un ítem politómico,
#' basado en la entropía de Shannon. El valor resultante refleja la dispersión
#' o utilización efectiva de las opciones de respuesta.
#'
#' Las categorías con frecuencia cero son tratadas como *ceros incidentales* (Linacre, 1999; 2002).
#' Para evitar problemas con logaritmos indefinidos, su frecuencia se reemplaza por un valor pequeño
#' definido en el argumento \code{pseudo_zero}, permitiendo continuar con el cálculo de la entropía.
#'
#' @param vector Un vector numérico, de factor o carácter con las respuestas de un ítem.
#' @param k Número total de opciones de respuesta posibles. Este valor es obligatorio
#'        y debe incluir todas las categorías teóricas, incluso si algunas no fueron utilizadas.
#' @param pseudo_zero Valor pequeño que se asigna a categorías con frecuencia cero,
#'        para evitar problemas con logaritmos indefinidos. Por defecto es \code{0.05}.
#'
#' @return Una lista con un solo elemento:
#' \item{AENO}{Valor numérico que indica el número efectivo de opciones utilizadas.}
#'
#' @examples
#' respuestas <- c(0, 1, 1, 2, 3, 3, 5, NA)
#' AENOuni(respuestas, k = 6)
#'
#' @references
#' Sato, T., & Morimoto, U. (1976). Sentaku-shi keishiki tesuto kaitou bunpu no bunseki [Analyzing endorsement distribution of selected-response items].
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
  if (missing(k)) stop("Debe especificar el número total de opciones posibles con el argumento 'k'.")

  prop_respuestas <- table(vector) / sum(!is.na(vector))

  categorias_completas <- as.character(0:(k - 1))
  prop_respuestas_completo <- rep(0, k)
  names(prop_respuestas_completo) <- categorias_completas
  prop_respuestas_completo[names(prop_respuestas)] <- prop_respuestas

  prop_respuestas_completo[prop_respuestas_completo == 0] <- pseudo_zero
  prop_respuestas_completo <- prop_respuestas_completo / sum(prop_respuestas_completo)

  OPLOG <- prop_respuestas_completo * log2(prop_respuestas_completo)
  SUMLOG <- -sum(OPLOG)
  AENO <- 2^SUMLOG

  return(list(AENO = AENO))
}
