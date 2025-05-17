#' @title catPlot: Plot of Category Means by Item
#'
#' @description
#' Creates a faceted plot showing the progression of mean total scores across
#' response categories for each item. This helps visually assess the monotonicity
#' of response categories in rating scales.
#'
#' The function expects the object returned by \code{\link{catOrder}}, which stores
#' a long-format data frame of category means as an internal attribute.
#'
#' @param data An object returned by \code{\link{catOrder}}. It must contain the
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
#' This function uses the internal \code{"means.long"} attribute stored by \code{catOrder}
#' and does not require the user to extract or prepare it manually.
#'
#' @examples
#' \dontrun{
#' # Fit diagnostics
#' res <- catOrder(data = your_data, score.total = rowSums(your_data), tests = "JT")
#'
#' # Default plot (all items)
#' catPlot(res)
#'
#' # Plot with LOESS smoothing
#' catPlot(res, smooth = TRUE)
#'
#' # Subset of items
#' catPlot(res, items = c("item1", "item3"))
#' }
#'
#' @import ggplot2
#' @export
catPlot <- function(data, items = NULL, smooth = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  df <- attr(data, "means.long")
  if (is.null(df)) stop("The object must be returned by 'catOrder()' and contain the 'means.long' attribute.")

  if (!is.null(items)) df <- df[df$Item %in% items, ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Category, y = Mean.Score, group = Item)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ Item, scales = "free_y") +
    ggplot2::labs(x = "Response Category", y = "Mean of Total Score", title = "Progression of Category Means by Item") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, linetype = "dashed")
  }

  return(p)
}

# Declare global variables for ggplot2
utils::globalVariables(c("Category", "Item", "Mean.Score"))
