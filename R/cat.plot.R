#' @title Plot of Category Means by Item
#'
#' @description 
#' Creates a panel plot showing the progression of mean total scores across response 
#' categories for each item. This helps visually assess the monotonicity of response 
#' categories, in line with rating scale diagnostics.
#'
#' The function expects as input the object returned by \code{\link{cat.order}}, 
#' which stores a long-format data frame of category means as an internal attribute.
#'
#' @param data An object returned by \code{\link{cat.order}}. It must contain the 
#' internal attribute \code{"means.long"}, which holds the long-format table of category means.
#' @param items Optional character vector with the names of specific items to plot. 
#' If \code{NULL} (default), all items are included.
#' @param smooth Logical. If \code{TRUE}, adds a dashed LOESS smoothing line to each plot. 
#' Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object with one facet per item. Each facet shows:
#' \itemize{
#'   \item The observed means of the total score for each response category (points and lines),
#'   \item An optional LOESS line to highlight trend consistency.
#' }
#'
#' @details 
#' This function uses the internal \code{"means.long"} attribute stored by \code{cat.order()} 
#' and does not require the user to extract or prepare it manually.
#'
#' @examples
#' \dontrun{
#' # Fit diagnostics
#' res <- cat.order(data = your_data, score.total = rowSums(your_data), tests = "JT")
#' 
#' # Default plot (all items)
#' cat.plot(res)
#'
#' # Plot with LOESS smoothing
#' cat.plot(res, smooth = TRUE)
#'
#' # Subset of items
#' cat.plot(res, items = c("item1", "item3"))
#' }
#'
#' @import ggplot2
#' 
#' @export
cat.plot <- function(data, items = NULL, smooth = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  
  # Recuperar 'means.long' desde atributo
  df <- attr(data, "means.long")
  
  if (is.null(df)) {
    stop("El objeto debe provenir de 'cat.order()' y contener el atributo 'means.long'.")
  }
  
  # Filtrar ítems si se especifica
  if (!is.null(items)) {
    df <- df[df$Item %in% items, ]
  }
  
  # Crear gráfico base
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Category, y = Mean.Score, group = Item)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ Item, scales = "free_y") +
    ggplot2::labs(
      x = "Response Category",
      y = "Mean of Total Score",
      title = "Progression of Category Means by Item"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Línea suavizada (opcional)
  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, linetype = "dashed")
  }
  
  return(p)
}
