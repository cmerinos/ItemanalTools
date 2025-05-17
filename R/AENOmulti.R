#' AENOmulti: Número efectivo de opciones para múltiples ítems
#'
#' Calcula el índice AENO (Número Efectivo de Opciones Nominales) para cada ítem 
#' en un conjunto de datos. Es útil para explorar la calidad de uso de las opciones 
#' de respuesta en escalas politómicas.
#'
#' Internamente llama a \code{\link{AENOuni}}. Las categorías con frecuencia cero 
#' son tratadas como *ceros incidentales* (Linacre, 1999; 2002) y su frecuencia es 
#' reemplazada por un valor pequeño (por defecto \code{0.05}) para evitar problemas 
#' en el cálculo del logaritmo durante la estimación de entropía.
#'
#' @param data Un \code{data.frame} con las respuestas a varios ítems (cada columna es un ítem).
#' @param k Número total de opciones de respuesta posibles (por ejemplo, \code{k = 6} si las opciones van de 0 a 5).
#'
#' @return Un \code{data.frame} con los siguientes campos para cada ítem:
#' \item{item}{Nombre del ítem.}
#' \item{AENO}{Número efectivo de opciones utilizadas (redondeado a 3 decimales).}
#' \item{cat.eff}{Número de categorías efectivamente utilizadas (observadas).}
#' \item{n}{Número de respuestas válidas (sin NA).}
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
AENOmulti <- function(data, k) {
  
  if (!is.data.frame(data)) {
    stop("La entrada debe ser un data.frame con ítems en columnas.")
  }
  if (missing(k)) {
    stop("Debe especificar el número total de opciones posibles con el argumento 'k'.")
  }
  
  resultados <- lapply(names(data), function(nombre_item) {
    vector <- data[[nombre_item]]
    vector_no_na <- vector[!is.na(vector)]
    cat_eff <- length(unique(vector_no_na))
    n_valido <- length(vector_no_na)
    
    # Calcular AENO redondeado a 3 decimales usando el valor explícito de k
    AENO <- round(AENOuni(vector, k = k)$AENO, 3)
    
    data.frame(
      item = nombre_item,
      AENO = AENO,
      cat.eff = cat_eff,
      n = n_valido,
      stringsAsFactors = FALSE
    )
  })
  
  df_resultados <- do.call(rbind, resultados)
  return(df_resultados)
}
