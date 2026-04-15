#' Information loss after recoding item response categories
#'
#' @description
#' Computes information-loss indices for comparing an original item-response dataset
#' with a recoded version of the same dataset. The function evaluates information loss
#' at three complementary levels: distribution, correlation structure, and latent structure.
#'
#' @param data A data.frame or matrix containing the original item responses.
#' Rows represent subjects and columns represent items.
#' @param data.recoded A data.frame or matrix containing the recoded item responses.
#' It must have the same number of rows and columns as `data`, and preserve the same
#' column order. Column names do not need to match those in `data`.
#' @param cor.type Type of correlation matrix used for correlation-structure and
#' latent-structure indices. Options are `"pearson"` or `"poly"`.
#' @param nfact Optional integer indicating the number of leading eigenvalues to be
#' retained in `Latent.structure`. If `NULL`, all eigenvalues are used.
#' @param log.base Base of the logarithm used to compute Shannon entropy.
#' Default is `exp(1)`.
#' @param use Character string indicating how missing values are handled for Pearson
#' correlations. Options are `"pairwise"` or `"complete"`.
#' @param correct Continuity correction used when `cor.type = "poly"`.
#' Passed to `psych::polychoric()`. Default is `0.5`.
#' @param smooth Logical value indicating whether smoothing should be applied when
#' `cor.type = "poly"` and the estimated polychoric matrix is not positive definite.
#' Passed to `psych::polychoric()`. Default is `TRUE`.
#' @param digits Number of decimal places used to round the results.
#'
#' @details
#' Information loss is assessed at three levels:
#'
#' \strong{Distribution.}
#' Shannon entropy is computed for each item before and after recoding. Absolute and
#' relative entropy loss are reported. Smaller values of `Relative.loss.H` indicate
#' lower distributional information loss. Implementation inspired from Hundepool et al. (2024)
#' and Templ et al. (2015)
#'
#' \strong{Correlation structure.}
#' The original and recoded inter-item correlation matrices are compared using the
#' root mean square error (`RMSE.R`), the mean absolute difference (`MAD.R`), the
#' correlation between the vectorized lower triangles of both matrices (`r.vec`),
#' and the `IL.correl` index. Smaller values of `RMSE.R`, `MAD.R`, and `IL.correl`,
#' together with larger values of `r.vec`, indicate better preservation of the
#' correlation structure.Implementation inspired from Hundepool et al. (2024)
#' and Templ et al. (2015)
#'
#' \strong{Latent structure.}
#' Eigenvalues from the original and recoded correlation matrices are compared to
#' quantify latent-structure distortion. Only the first `nfact` eigenvalues are
#' reported. The `Summary` row aggregates the retained eigenvalues. Smaller values
#' of `RL.lambda` indicate lower latent-structure information loss.
#'
#' If `cor.type = "poly"`, the function computes polychoric correlations using
#' `psych::polychoric()` with the specified continuity correction (`correct`) and
#' smoothing option (`smooth`).
#'
#'@references
#'Mlodak, A. (2019). Using the Complex Measure in an Assessment of the Information Loss Due to the Microdata Disclosure Control. Przegląd Statystyczny, 2019, 66(1), 7-26. DOI: 10.5604/01.3001.0013.8285
#'
#'Mlodak, A. (2020). Information loss resulting from statistical disclosure control of output data, Wiadomosci Statystyczne. The Polish Statistician, 2020, 65(9), 7-27, DOI: 10.5604/01.3001.0014.4121
#'
#'Templ, M., Kowarik, A., & Meindl, B. (2015). Statistical Disclosure Control for Micro-Data Using the R Package sdcMicro. Journal of Statistical Software, 67(4), 1-36. doi:10.18637/jss.v067.i04.
#'
#'Hundepool, A., Domingo-Ferrer, J., Franconi, L., Giessing, S., Lenz, R., Naylor, J., Nordholt, E. S., Seri, G., De Wolf, P., Tent, R., Młodak, A., Gussenbauer, J., & Wilak, K. (2024). Handbook on Statistical Disclosure Control.https://julienjamme.github.io/handbook_sdc_from_doc_to_md/03-microdata.html#sec-informationloss-microdata
#'
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{Distribution}: data.frame with item-level entropy indices and a final
#'   row with total entropy values.
#'   \item \code{Correlation.structure}: data.frame with global indices comparing the
#'   original and recoded correlation matrices.
#'   \item \code{Latent.structure}: data.frame with eigenvalue-based indices for the
#'   retained components and a final summary row.
#' }
#'
#' @examples
#' data.orig <- data.frame(
#'   Item1 = c(1, 2, 3, 4, 4),
#'   Item2 = c(1, 1, 2, 3, 4),
#'   Item3 = c(2, 2, 3, 3, 4)
#' )
#'
#' data.rec <- data.frame(
#'   Item1 = c(1, 1, 2, 2, 2),
#'   Item2 = c(1, 1, 1, 2, 2),
#'   Item3 = c(1, 1, 2, 2, 2)
#' )
#'
#' infoLoss(data = data.orig, data.recoded = data.rec, nfact = 2)
#'
#' @export
infoLoss <- function(data,
                     data.recoded,
                     cor.type = c("pearson", "poly"),
                     nfact = NULL,
                     log.base = exp(1),
                     use = c("pairwise", "complete"),
                     correct = 0.5,
                     smooth = TRUE,
                     digits = 3) {
  
  cor.type <- match.arg(cor.type)
  use <- match.arg(use)
  
  if (is.matrix(data)) data <- as.data.frame(data)
  if (is.matrix(data.recoded)) data.recoded <- as.data.frame(data.recoded)
  
  if (!is.data.frame(data) || !is.data.frame(data.recoded)) {
    stop("`data` and `data.recoded` must be data.frames or matrices.")
  }
  
  if (nrow(data) != nrow(data.recoded)) {
    stop("`data` and `data.recoded` must have the same number of rows.")
  }
  
  if (ncol(data) != ncol(data.recoded)) {
    stop("`data` and `data.recoded` must have the same number of columns.")
  }
  
  p <- ncol(data)
  
  if (p < 2) stop("At least two items are required.")
  
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("Item", seq_len(p))
  }
  
  if (is.null(colnames(data.recoded))) {
    colnames(data.recoded) <- paste0("Item.recoded", seq_len(p))
  }
  
  non.numeric.data <- !vapply(data, is.numeric, logical(1))
  non.numeric.rec  <- !vapply(data.recoded, is.numeric, logical(1))
  if (any(non.numeric.data) || any(non.numeric.rec)) {
    stop("All columns in `data` and `data.recoded` must be numeric.")
  }
  
  if (!is.null(nfact)) {
    if (length(nfact) != 1L || !is.numeric(nfact) || is.na(nfact)) {
      stop("`nfact` must be a single positive integer or NULL.")
    }
    nfact <- as.integer(nfact)
    if (nfact < 1L || nfact > p) {
      stop("`nfact` must be between 1 and the number of items.")
    }
  } else {
    nfact <- p
  }
  
  if (!is.numeric(log.base) || length(log.base) != 1L ||
      is.na(log.base) || log.base <= 0 || log.base == 1) {
    stop("`log.base` must be a single positive number different from 1.")
  }
  
  if (!is.numeric(correct) || length(correct) != 1L ||
      is.na(correct) || correct < 0) {
    stop("`correct` must be a single non-negative numeric value.")
  }
  
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth)) {
    stop("`smooth` must be TRUE or FALSE.")
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative number.")
  }
  digits <- as.integer(digits)
  
  .entropy.item <- function(x, log.base = exp(1)) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_real_)
    probs <- prop.table(table(x))
    probs <- as.numeric(probs)
    -sum(probs * (log(probs) / log(log.base)))
  }
  
  .lower.tri.vec <- function(R) {
    R[lower.tri(R, diag = FALSE)]
  }
  
  .compute.cor <- function(dat, cor.type, use, correct, smooth) {
    if (cor.type == "pearson") {
      use.cor <- if (use == "pairwise") "pairwise.complete.obs" else "complete.obs"
      return(stats::cor(dat, use = use.cor, method = "pearson"))
    }
    
    if (!requireNamespace("psych", quietly = TRUE)) {
      stop("Package 'psych' is required when `cor.type = \"poly\"`.")
    }
    
    out <- tryCatch(
      psych::polychoric(dat, correct = correct, smooth = smooth),
      error = function(e) {
        stop(
          "Polychoric correlation matrix could not be estimated. ",
          "Check for sparse categories, zero-variance items, or try different ",
          "`correct`/`smooth` settings."
        )
      }
    )
    
    out$rho
  }
  
  .IL.correl <- function(R1, R2) {
    inv1 <- tryCatch(solve(R1), error = function(e) NULL)
    inv2 <- tryCatch(solve(R2), error = function(e) NULL)
    
    if (is.null(inv1) || is.null(inv2)) return(NA_real_)
    
    d1 <- diag(inv1)
    d2 <- diag(inv2)
    
    if (anyNA(d1) || anyNA(d2)) return(NA_real_)
    if (sum(d1^2) == 0 || sum(d2^2) == 0) return(NA_real_)
    
    d1 <- d1 / sqrt(sum(d1^2))
    d2 <- d2 / sqrt(sum(d2^2))
    
    (1 / sqrt(2)) * sqrt(sum((d1 - d2)^2))
  }
  
  .safe.round <- function(x, digits) {
    if (is.numeric(x)) round(x, digits) else x
  }
  
  ## Distribution
  H.original <- vapply(data, .entropy.item, numeric(1), log.base = log.base)
  H.recoded  <- vapply(data.recoded, .entropy.item, numeric(1), log.base = log.base)
  
  Delta.H <- H.original - H.recoded
  Relative.loss.H <- ifelse(H.original == 0, NA_real_, Delta.H / H.original)
  
  Distribution <- data.frame(
    Item = colnames(data),
    H.original = H.original,
    H.recoded = H.recoded,
    Delta.H = Delta.H,
    Relative.loss.H = Relative.loss.H,
    stringsAsFactors = FALSE
  )
  
  Distribution <- rbind(
    Distribution,
    data.frame(
      Item = "Total",
      H.original = sum(H.original, na.rm = TRUE),
      H.recoded = sum(H.recoded, na.rm = TRUE),
      Delta.H = sum(Delta.H, na.rm = TRUE),
      Relative.loss.H = ifelse(sum(H.original, na.rm = TRUE) == 0,
                               NA_real_,
                               sum(Delta.H, na.rm = TRUE) / sum(H.original, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  )
  
  ## Correlation structure
  R.original <- .compute.cor(data, cor.type = cor.type, use = use,
                             correct = correct, smooth = smooth)
  R.recoded  <- .compute.cor(data.recoded, cor.type = cor.type, use = use,
                             correct = correct, smooth = smooth)
  
  vec.original <- .lower.tri.vec(R.original)
  vec.recoded  <- .lower.tri.vec(R.recoded)
  
  ok <- stats::complete.cases(vec.original, vec.recoded)
  
  if (!any(ok)) {
    RMSE.R <- NA_real_
    MAD.R <- NA_real_
    r.vec <- NA_real_
  } else {
    diffs <- vec.original[ok] - vec.recoded[ok]
    RMSE.R <- sqrt(mean(diffs^2))
    MAD.R <- mean(abs(diffs))
    r.vec <- if (sum(ok) < 2) NA_real_ else stats::cor(vec.original[ok], vec.recoded[ok])
  }
  
  IL.correl <- .IL.correl(R.original, R.recoded)
  
  Correlation.structure <- data.frame(
    cor.type = cor.type,
    RMSE.R = RMSE.R,
    MAD.R = MAD.R,
    r.vec = r.vec,
    IL.correl = IL.correl,
    stringsAsFactors = FALSE
  )
  
  ## Latent structure
  eig.original <- sort(eigen(R.original, symmetric = TRUE, only.values = TRUE)$values,
                       decreasing = TRUE)
  eig.recoded  <- sort(eigen(R.recoded, symmetric = TRUE, only.values = TRUE)$values,
                       decreasing = TRUE)
  
  idx <- seq_len(nfact)
  
  eig.original <- eig.original[idx]
  eig.recoded  <- eig.recoded[idx]
  
  Delta.lambda <- eig.original - eig.recoded
  Abs.delta.lambda <- abs(Delta.lambda)
  
  Prop.original <- eig.original / p
  Prop.recoded  <- eig.recoded / p
  
  RL.lambda <- ifelse(eig.original == 0, NA_real_, Delta.lambda / eig.original)
  
  Latent.structure <- data.frame(
    Component = paste0("Comp", idx),
    Eigen.original = eig.original,
    Eigen.recoded = eig.recoded,
    Delta.lambda = Delta.lambda,
    Abs.delta.lambda = Abs.delta.lambda,
    Prop.original = Prop.original,
    Prop.recoded = Prop.recoded,
    RL.lambda = RL.lambda,
    stringsAsFactors = FALSE
  )
  
  Latent.structure <- rbind(
    Latent.structure,
    data.frame(
      Component = "Summary",
      Eigen.original = sum(eig.original, na.rm = TRUE),
      Eigen.recoded = sum(eig.recoded, na.rm = TRUE),
      Delta.lambda = sum(Delta.lambda, na.rm = TRUE),
      Abs.delta.lambda = sum(Abs.delta.lambda, na.rm = TRUE),
      Prop.original = sum(Prop.original, na.rm = TRUE),
      Prop.recoded = sum(Prop.recoded, na.rm = TRUE),
      RL.lambda = ifelse(sum(eig.original, na.rm = TRUE) == 0,
                         NA_real_,
                         (sum(eig.original, na.rm = TRUE) - sum(eig.recoded, na.rm = TRUE)) /
                           sum(eig.original, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  )
  
  Distribution[] <- lapply(Distribution, .safe.round, digits = digits)
  Correlation.structure[] <- lapply(Correlation.structure, .safe.round, digits = digits)
  Latent.structure[] <- lapply(Latent.structure, .safe.round, digits = digits)
  
  list(
    Distribution = Distribution,
    Correlation.structure = Correlation.structure,
    Latent.structure = Latent.structure
  )
}