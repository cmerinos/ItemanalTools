#' @title Category Order Diagnostics for Polytomous Items
#'
#' @description 
#' Evaluates whether the response categories of Likert-type items show a monotonic 
#' progression along a total score or latent trait proxy. For each item, it computes 
#' the mean of the total score by response category and applies selected nonparametric 
#' trend tests.
#'
#' @param data A data frame containing polytomous items (e.g., with values 0–4 or 1–5). 
#' Each column represents an item.
#' @param score.total A numeric vector with the total or composite score for each respondent. 
#' Typically the row sum of the item data.
#' @param tests Character vector specifying which trend tests to compute. Options are:
#' \code{"JT"} for the Jonckheere-Terpstra test and \code{"Cu"} for Cuzick's test. 
#' Default is \code{c("JT", "Cu")}.
#'
#' @return A list with two visible components:
#' \item{summary}{A data frame with one row per item, including: item name, monotonicity flag,
#' number of violations, Z-statistics, p-values, and effect sizes for each selected test.}
#' \item{means.wide}{A data frame in wide format. Each row corresponds to an item and each column
#' contains the mean score for a given response category (e.g., \code{cat0.mean}, \code{cat1.mean}, ...).}
#'
#' Internally, the returned object stores a third component \code{"means.long"} as an attribute. 
#' This is not printed automatically but can be accessed manually using:
#' \code{attr(result, "means.long")}. This component is used internally by functions such as \code{\link{cat.plot}}.
#'
#' @examples
#' \dontrun{
#' # Simulated 5-point Likert items
#' set.seed(123)
#' data <- data.frame(
#'   item1 = sample(0:4, 100, TRUE),
#'   item2 = sample(0:4, 100, TRUE),
#'   item3 = sample(0:4, 100, TRUE)
#' )
#' score <- rowSums(data)
#' 
#' res <- cat.order(data, score.total = score, tests = "JT")
#' 
#' # View results
#' res$summary       # Summary table
#' res$means.wide    # Category means in wide format
#' 
#' # Access internal long-format table (for plotting)
#' head(attr(res, "means.long"))
#' }
#'
#' @importFrom PMCMRplus jonckheereTest cuzickTest
#' 
#' @export
cat.order <- function(data, score.total, tests = c("JT", "Cu")) {
  if (!requireNamespace("PMCMRplus", quietly = TRUE)) stop("Package 'PMCMRplus' is required.")
  
  # Diccionario de columnas para cada test
  cols.by.test <- list(
    JT = c("JT.Z", "JT.p", "JT.ES"),
    Cu = c("Cu.Z", "Cu.p", "Cu.ES")
  )
  
  summary.list <- list()
  means.list <- list()
  means.wide <- list()
  
  for (item in names(data)) {
    responses <- data[[item]]
    valid <- !is.na(responses) & !is.na(score.total)
    x <- score.total[valid]
    g <- as.factor(responses[valid])
    N <- length(x)
    
    # Medias por categoría
    cat.means <- tapply(x, g, mean, na.rm = TRUE)
    is.mono <- all(diff(cat.means) >= 0)
    violations <- sum(diff(cat.means) < 0)
    
    # Guardar medias
    cat.wide <- as.list(round(cat.means, 3))
    names(cat.wide) <- paste0("cat", names(cat.means), ".mean")
    cat.wide$Item <- item
    
    means.list[[item]] <- data.frame(
      Item = item,
      Category = as.numeric(names(cat.means)),
      Mean.Score = as.numeric(cat.means)
    )
    
    # Inicializar estadísticos
    jt.z <- jt.p <- jt.es <- NA_real_
    cu.z <- cu.p <- cu.es <- NA_real_
    
    # Test JT
    if ("JT" %in% tests) {
      jt <- suppressWarnings(PMCMRplus::jonckheereTest(x, g))
      jt.z <- as.numeric(jt$statistic)
      jt.p <- jt$p.value
      jt.es <- jt.z / sqrt(N)
    }
    
    # Test Cu
    if ("Cu" %in% tests) {
      cu <- suppressWarnings(PMCMRplus::cuzickTest(x, g))
      cu.z <- as.numeric(cu$statistic)
      cu.p <- cu$p.value
      cu.es <- cu.z / sqrt(N)
    }
    
    row <- data.frame(
      Item = item,
      Monotonic = ifelse(is.mono, "Yes", "No"),
      Violations = violations,
      JT.Z = jt.z,
      JT.p = jt.p,
      JT.ES = jt.es,
      Cu.Z = cu.z,
      Cu.p = cu.p,
      Cu.ES = cu.es,
      stringsAsFactors = FALSE
    )
    
    summary.list[[item]] <- row
    means.wide[[item]] <- cat.wide
  }
  
  # Combinar resultados
  summary.df <- do.call(rbind, summary.list)
  columns.to.keep <- c("Item", "Monotonic", "Violations", unlist(cols.by.test[tests]))
  summary.df <- summary.df[, columns.to.keep, drop = FALSE]
  
  means.wide.df <- do.call(rbind, lapply(means.wide, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  means.long.df <- do.call(rbind, means.list)
  
  # Armar lista y ocultar means.long como atributo
  out <- list(
    summary = summary.df,
    means.wide = means.wide.df
  )
  attr(out, "means.long") <- means.long.df
  
  return(out)
}
