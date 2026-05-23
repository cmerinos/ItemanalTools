#' @title itemPlot: Plot of item coefficients with confidence intervals
#'
#' @description
#' Creates a graphical display of item-level coefficients with their corresponding
#' confidence intervals. The function supports optional grouping, item ordering,
#' reference thresholds, and customizable value labels.
#'
#' This function is designed for outputs from item analysis procedures where each
#' row corresponds to an item and includes a coefficient estimate and its lower
#' and upper confidence interval bounds.
#'
#' @param data A data.frame containing item-level results.
#'
#' @param item.col Character. Name of the column containing item labels.
#'
#' @param value.col Character. Name of the column containing the coefficient estimates.
#'
#' @param lwr.ci.col Character. Name of the column containing the lower bound of the confidence interval.
#'
#' @param upr.ci.col Character. Name of the column containing the upper bound of the confidence interval.
#'
#' @param group Optional. Either:
#' \itemize{
#'   \item A character string indicating a column in `data`, or
#'   \item A vector of group labels with length equal to the number of rows in `data`.
#' }
#' When provided, coefficients and intervals are plotted separately by group.
#'
#' @param threshold Optional numeric value. Adds a reference line (vertical or horizontal depending on orientation).
#'
#' @param sort.items Character. Sorting method for items:
#' \describe{
#'   \item{"none"}{Keeps original order.}
#'   \item{"asc"}{Sorts items in ascending order of the coefficient.}
#'   \item{"desc"}{Sorts items in descending order.}
#' }
#'
#' @param orientation Character. Plot orientation:
#' \describe{
#'   \item{"h"}{Horizontal (recommended).}
#'   \item{"v"}{Vertical.}
#' }
#'
#' @param theme.style Character. Plot theme:
#' \describe{
#'   \item{"academic"}{Clean academic-style theme (default).}
#'   \item{"classic"}{Classic ggplot2 theme.}
#'   \item{"minimal"}{Minimal ggplot2 theme.}
#' }
#'
#' @param show.values Logical. If TRUE, displays value labels on the plot.
#'
#' @param value.size Numeric. Font size of value labels.
#'
#' @param value.digits Integer. Number of decimal places in labels.
#'
#' @param value.format Character. Format of value labels:
#' \describe{
#'   \item{"estimate"}{Shows only the coefficient.}
#'   \item{"ci"}{Shows only the confidence interval.}
#'   \item{"both"}{Shows coefficient and confidence interval.}
#' }
#'
#' @param value.hjust Numeric. Horizontal justification of value labels.
#'
#' @param value.vjust Numeric. Vertical justification of value labels.
#'
#' @param xlim Optional numeric vector of length 2. Sets the limits of the x-axis (for horizontal plots).
#'
#' @param point.size Numeric. Size of the point representing the coefficient.
#'
#' @param errorbar.width Numeric. Width (thickness) of the confidence interval bars.
#'
#' @param dodge.width Numeric. Controls separation between groups.
#'
#' @param point.color Character. Color of points (when no grouping is used).
#'
#' @param bar.color Character. Color of confidence interval bars (when no grouping is used).
#'
#' @param line.color Character. Color of the reference threshold line.
#'
#' @param xlab Optional character. Label for x-axis.
#'
#' @param ylab Optional character. Label for y-axis.
#'
#' @param title Optional character. Plot title.
#'
#' @details
#' When `group` is provided, value labels (if enabled) are displayed for each
#' group-specific estimate. Labels are positioned over the corresponding confidence
#' interval lines to improve readability.
#'
#' The horizontal orientation is recommended, especially when `show.values = TRUE`,
#' as it improves label readability and avoids overlap.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' # Example without grouping
#' res <- data.frame(
#'   Item = paste0("Item", 1:5),
#'   coef = c(.62, .55, .48, .41, .36),
#'   lwr.ci = c(.54, .47, .39, .32, .27),
#'   upr.ci = c(.70, .63, .57, .50, .45)
#' )
#'
#' itemPlot(res,
#'          show.values = TRUE,
#'          value.format = "both",
#'          threshold = 0.40)
#'
#' # Example with grouping
#' res_g <- data.frame(
#'   Item = rep(paste0("Item", 1:5), 2),
#'   group = rep(c("A", "B"), each = 5),
#'   coef = c(.62, .55, .48, .41, .36,
#'            .58, .50, .45, .39, .33),
#'   lwr.ci = c(.54, .47, .39, .32, .27,
#'              .49, .42, .36, .30, .24),
#'   upr.ci = c(.70, .63, .57, .50, .45,
#'              .67, .58, .54, .48, .42)
#' )
#'
#' itemPlot(res_g,
#'          group = "group",
#'          show.values = TRUE,
#'          value.format = "both",
#'          threshold = 0.40)
#'
#' @export
ItemPlot <- function(data,
                        item.col = "Item",
                        value.col = "coef",
                        lwr.ci.col = "lwr.ci",
                        upr.ci.col = "upr.ci",
                        group = NULL,
                        threshold = NULL,
                        sort.items = c("none", "asc", "desc"),
                        orientation = c("h", "v"),
                        theme.style = c("academic", "classic", "minimal"),
                        show.values = FALSE,
                        value.size = 3,
                        value.digits = 2,
                        value.format = c("estimate", "ci", "both"),
                        value.hjust = 0.5,
                        value.vjust = -0.6,
                        xlim = NULL,
                        point.size = 2.5,
                        errorbar.width = 0.20,
                        dodge.width = 0.50,
                        point.color = "black",
                        bar.color = "gray40",
                        line.color = "red",
                        xlab = NULL,
                        ylab = NULL,
                        title = NULL) {
  
  # =========================
  # Checks
  # =========================
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  sort.items  <- match.arg(sort.items)
  orientation <- match.arg(orientation)
  theme.style <- match.arg(theme.style)
  value.format <- match.arg(value.format)
  
  needed.cols <- c(item.col, value.col, lwr.ci.col, upr.ci.col)
  if (!all(needed.cols %in% names(data))) {
    stop("Missing required columns.")
  }
  
  df <- data
  
  # =========================
  # Group handling
  # =========================
  group.name <- NULL
  
  if (!is.null(group)) {
    if (is.character(group) && group %in% names(df)) {
      group.name <- group
    } else {
      if (length(group) != nrow(df)) stop("Invalid group length.")
      df$.group <- as.factor(group)
      group.name <- ".group"
    }
    df[[group.name]] <- as.factor(df[[group.name]])
  }
  
  # =========================
  # Ordering
  # =========================
  df[[item.col]] <- as.character(df[[item.col]])
  
  if (sort.items != "none") {
    if (is.null(group.name)) {
      ord <- df[order(df[[value.col]]), item.col]
    } else {
      tmp <- aggregate(df[[value.col]], list(df[[item.col]]), mean)
      ord <- tmp[order(tmp$x), 1]
    }
    if (sort.items == "desc") ord <- rev(ord)
    df[[item.col]] <- factor(df[[item.col]], levels = ord)
  } else {
    df[[item.col]] <- factor(df[[item.col]], levels = unique(df[[item.col]]))
  }
  
  # =========================
  # Labels
  # =========================
  if (show.values) {
    fmt <- function(x) format(round(x, value.digits), nsmall = value.digits)
    
    est <- fmt(df[[value.col]])
    lwr <- fmt(df[[lwr.ci.col]])
    upr <- fmt(df[[upr.ci.col]])
    
    df$.label <- switch(
      value.format,
      estimate = est,
      ci = paste0("[", lwr, ", ", upr, "]"),
      both = paste0(est, " [", lwr, ", ", upr, "]")
    )
  }
  
  # =========================
  # Themes
  # =========================
  theme_academic <- function() {
    ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }
  
  theme_selected <- switch(
    theme.style,
    academic = theme_academic(),
    classic = ggplot2::theme_classic(),
    minimal = ggplot2::theme_minimal()
  )
  
  dodge <- ggplot2::position_dodge(width = dodge.width)
  
  # =========================
  # Plot
  # =========================
  if (orientation == "h") {
    
    if (is.null(group.name)) {
      
      p <- ggplot2::ggplot(df,
                           ggplot2::aes(x = .data[[value.col]], y = .data[[item.col]])) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = .data[[lwr.ci.col]], xmax = .data[[upr.ci.col]]),
          height = errorbar.width,
          color = bar.color
        ) +
        ggplot2::geom_point(size = point.size, color = point.color)
      
      if (show.values) {
        p <- p + ggplot2::geom_text(
          ggplot2::aes(label = .data[[".label"]]),
          hjust = value.hjust,
          vjust = value.vjust,
          size = value.size
        )
      }
      
    } else {
      
      p <- ggplot2::ggplot(df,
                           ggplot2::aes(x = .data[[value.col]],
                                        y = .data[[item.col]],
                                        color = .data[[group.name]])) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = .data[[lwr.ci.col]],
                       xmax = .data[[upr.ci.col]]),
          height = errorbar.width,
          position = dodge
        ) +
        ggplot2::geom_point(size = point.size, position = dodge)
      
      if (show.values) {
        p <- p + ggplot2::geom_text(
          ggplot2::aes(label = .data[[".label"]]),
          position = dodge,
          hjust = value.hjust,
          vjust = value.vjust,
          size = value.size,
          show.legend = FALSE
        )
      }
    }
    
    if (!is.null(threshold)) {
      p <- p + ggplot2::geom_vline(xintercept = threshold,
                                   linetype = "dashed",
                                   color = line.color)
    }
    
    if (!is.null(xlim)) {
      p <- p + ggplot2::coord_cartesian(xlim = xlim)
    }
    
    p <- p + ggplot2::labs(x = "Coefficient", y = NULL, title = title) +
      theme_selected
    
  } else {
    
    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = .data[[item.col]], y = .data[[value.col]])) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lwr.ci.col]],
                     ymax = .data[[upr.ci.col]])
      )
    
    if (show.values) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[[".label"]]),
        hjust = value.hjust,
        vjust = value.vjust,
        size = value.size
      )
    }
    
    if (!is.null(threshold)) {
      p <- p + ggplot2::geom_hline(yintercept = threshold,
                                   linetype = "dashed",
                                   color = line.color)
    }
    
    p <- p + theme_selected
  }
  
  return(p)
}