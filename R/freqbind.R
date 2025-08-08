#' @title Frequency and Proportion Table for Items
#' 
#' @description
#' Computes frequency tables for each item (column) in a data frame, with optional proportion tables.
#' Allows setting specific response levels and includes missing values if desired.
#'
#' @param data A data frame where each column is an item (e.g., from a questionnaire).
#' @param levels A vector of expected response categories (e.g., \code{1:5}). If \code{NULL}, all unique values are used.
#' @param show.na Logical. If \code{TRUE}, includes missing values in the frequency table.
#' @param na.label Character. Label to use for missing values (default is \code{"Missing"}).
#' @param proportion Logical. If \code{TRUE}, also computes proportions per column.
#' @param digits Integer. Number of decimals to round the proportions (default is \code{3}).
#'
#' @return
#' If \code{proportion = FALSE}, returns a data frame of frequencies.
#' If \code{proportion = TRUE}, returns a list with:
#' \itemize{
#'   \item \code{$freq}: Frequency table.
#'   \item \code{$prop}: Proportion table.
#' }
#'
#' @examples
#' data <- data.frame(item1 = c(1, 2, 2, 3, 3, NA),
#'                    item2 = c(1, 1, 2, 2, 3, 3))
#' freq.bind(data, levels = 1:3, proportion = TRUE)
#'
#' @export
freqbind <- function(data, levels = NULL, show.na = TRUE, na.label = "Missing",
                      proportion = FALSE, digits = 3) {
  if (!is.data.frame(data)) stop("El argumento debe ser un data.frame.")
  
  # Detectar niveles si no se especifican
  if (is.null(levels)) {
    levels <- sort(unique(unlist(data)))
    if (show.na && any(is.na(levels))) levels <- levels[!is.na(levels)]
  }
  
  # Crear lista de tablas por columna
  freq.list <- lapply(data, function(x) {
    x <- factor(x, levels = levels, exclude = NULL)
    table(x, useNA = if (show.na) "ifany" else "no")
  })
  
  # Unir en una tabla
  freq.table <- do.call(cbind, freq.list)
  freq.df <- as.data.frame(freq.table)
  colnames(freq.df) <- names(data)
  
  # Renombrar filas (NA a "Missing")
  rn <- rownames(freq.df)
  rn[is.na(rn)] <- na.label
  rownames(freq.df) <- rn
  
  # Si no se solicita proporciones
  if (!proportion) {
    return(freq.df)
  }
  
  # Calcular proporciones
  prop.mat <- prop.table(as.matrix(freq.df), margin = 2)
  prop.df <- round(as.data.frame(prop.mat), digits = digits)
  rownames(prop.df) <- rownames(freq.df)
  colnames(prop.df) <- names(data)
  
  # Devolver ambos resultados
  return(list(freq = freq.df, prop = prop.df))
}
