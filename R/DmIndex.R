#' @title DmIndex: Deviation Index from a Theoretical Association
#'
#' @description
#' Computes the **Deviation Index (Dm)** for each item, defined as:
#' \deqn{Dm = \sqrt{1 - (r - r_{ref})^2}}
#' or, if `metric="fisher"`,
#' \deqn{Dm = \sqrt{1 - (z(r) - z(r_{ref}))^2}},
#' where \eqn{z()} is the Fisher transformation.
#'
#' Confidence intervals are estimated by bootstrap in the Fisher z
#' domain, then transformed back to the Dm scale.
#'
#' @param df A data frame containing item-level results.
#' @param coef.col Character. Name of the column with observed correlations (e.g. `"rg"`, `"tau"`, `"r"`).
#' @param ref.value Numeric. Reference value(s). May be a single scalar (applied to all items)
#'   or a vector of the same length as `nrow(df)`. Default = 0.50.
#' @param metric Character. Scale for computation:
#'   * `"r"`: raw correlation scale
#'   * `"fisher"`: Fisher z scale (default).
#' @param ci.method Character. Confidence interval method: `"none"` or `"zboot"`.
#'   * `"none"`: only point estimate of Dm is reported.
#'   * `"zboot"`: **parametric bootstrap** confidence interval, assuming the
#'     Fisher z transformation of the correlation is normally distributed with
#'     variance 1/(n.eff - 3). Requires `n.eff`.
#' @param conf.level Confidence level for bootstrap CI. Default = 0.95.
#' @param B Integer. Number of bootstrap replicates for `"zboot"`. Default = 2000.
#' @param n.eff Integer. Effective sample size used to approximate SE in Fisher z.
#'   Required for bootstrap (`ci.method="zboot"`).
#' @param round.digits Integer. Number of decimal places to report. Default = 3.
#' @param clamp.zero Logical. If TRUE and `metric="r"`, truncate negative Dm at 0.
#'   Default = TRUE.
#' @param item.col Character. Name of the column with item labels. Default `"Item"`.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item `Item`: item name.
#'   \item `Dm`: Deviation Index point estimate.
#'   \item `Dm.lwr`, `Dm.upr`: lower and upper bootstrap CI bounds (if `ci.method="zboot"`).
#'   \item `coef.mean`: observed coefficient.
#'   \item `abs.diff.mean`: absolute difference between observed and reference value
#'     (in `r` or Fisher z scale, depending on `metric`).
#' }
#'
#' @details
#' The Deviation Index (Taras & Kline, 2010) quantifies the proximity of
#' an observed correlation to a theoretical reference. Values closer to 1 indicate
#' better agreement with the reference.
#'
#' The `"zboot"` option performs a **parametric bootstrap** in the Fisher z domain:
#' for each replicate, a Fisher z value is sampled from
#' \eqn{N(\hat{z}, 1/(n.eff-3))}, transformed back, and used to compute Dm.
#' The empirical quantiles across replicates form the confidence interval.
#' This is not a nonparametric bootstrap on raw data, but an approximation
#' suitable when only correlation coefficients are available.
#'
#' @references
#' Dominguez Lara, S. (2013). Propuesta para el cálculo del índice DM para la validez de los
#' ítemes. \emph{Interdisciplinaria Revista de Psicología y Ciencias Afines}, 30(2), 297–303.
#' https://doi.org/10.16888/interd.2013.30.2.8
#'
#' Taras, V., & Kline, T. (2010). Scale validation via quantifying item validity using the
#' Dm index. \emph{Psychological Reports}, 107(2), 535–546. https://doi.org/10.2466/03.PR0.107.5.535-546
#'
#'@importFrom stats rnorm
#'
#' @examples
#' df <- data.frame(
#'   Item = paste0("Item", 1:5),
#'   tau = c(.22, .18, .35, .12, .29)
#' )
#' # Only point estimate, for a dataframe ("df", with a target colunm "tau")
#' DmIndex(df, coef.col = "tau", ref.value = 0.30, ci.method = "none")
#'
#' # With parametric bootstrap CI (requires n.eff)
#' DmIndex(df, coef.col = "tau", ref.value = 0.30,
#'         ci.method = "zboot", n.eff = 250, B = 500)
#'
#' @export
DmIndex <- function(df,
                    coef.col,
                    ref.value = 0.50,
                    metric = c("fisher", "r"),
                    ci.method = c("none", "zboot"),
                    conf.level = 0.95,
                    B = 2000,
                    n.eff = NULL,
                    round.digits = 3,
                    clamp.zero = TRUE,
                    item.col = "Item") {

  metric <- match.arg(metric)
  ci.method <- match.arg(ci.method)

  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!coef.col %in% names(df)) stop("`coef.col` not found in `df`.")

  n <- nrow(df)
  # Handle item labels
  Items <- if (item.col %in% names(df)) df[[item.col]] else paste0("Item", seq_len(n))

  # Reference value vector
  if (length(ref.value) == 1L) {
    ref.vec <- rep(ref.value, n)
  } else if (length(ref.value) == n) {
    ref.vec <- ref.value
  } else {
    stop("`ref.value` must be scalar or of length equal to nrow(df).")
  }

  r <- df[[coef.col]]

  # Fisher transform
  zfun <- function(x) atanh(pmin(pmax(x, -0.999999), 0.999999))

  # Compute Dm (item-level)
  Dm_fun <- function(ri, refi, use_z, clamp0) {
    if (use_z) {
      zr <- zfun(ri); zref <- zfun(refi)
      sqrt(1 - (zr - zref)^2)
    } else {
      out <- sqrt(1 - (ri - refi)^2)
      if (clamp0) out <- pmax(0, out)
      out
    }
  }

  Dm_obs <- mapply(function(ri, refi) Dm_fun(ri, refi, metric == "fisher", clamp.zero), r, ref.vec)

  # Bootstrap CI (item-level)
  Dm_lwr <- Dm_upr <- rep(NA_real_, n)
  if (ci.method == "zboot") {
    if (is.null(n.eff)) stop("`n.eff` (effective sample size) is required for bootstrap.")
    q_hi <- (1 + conf.level) / 2
    zq   <- qnorm(q_hi)

    for (i in seq_len(n)) {
      ri <- r[i]; refi <- ref.vec[i]
      if (is.na(ri) || is.na(refi)) next
      zhat <- zfun(ri)
      se   <- 1 / sqrt(n.eff - 3)
      if (!is.finite(se) || se <= 0) next
      zsim <- rnorm(B, mean = zhat, sd = se)
      if (metric == "fisher") {
        Dsim <- sqrt(1 - (zsim - zfun(refi))^2)
      } else {
        rsim <- tanh(zsim)
        Dsim <- Dm_fun(rsim, refi, use_z = FALSE, clamp0 = clamp.zero)
      }
      qs <- quantile(Dsim, c(1 - q_hi, q_hi), na.rm = TRUE, names = FALSE)
      Dm_lwr[i] <- qs[1]; Dm_upr[i] <- qs[2]
    }
  }

  out <- data.frame(
    Item = Items,
    Dm = round(Dm_obs, round.digits),
    Dm.lwr = if (all(is.na(Dm_lwr))) NA_real_ else round(Dm_lwr, round.digits),
    Dm.upr = if (all(is.na(Dm_upr))) NA_real_ else round(Dm_upr, round.digits),
    coef.mean = round(r, round.digits),
    abs.diff.mean = if (metric == "fisher") {
      round(abs(zfun(r) - zfun(ref.vec)), round.digits)
    } else {
      round(abs(r - ref.vec), round.digits)
    },
    stringsAsFactors = FALSE
  )

  ## ---- NUEVO: Dm.total (siempre) ----
  X <- if (metric == "fisher") zfun(r) else r
  T <- if (metric == "fisher") zfun(ref.vec) else ref.vec
  rmse_total <- sqrt(mean((X - T)^2, na.rm = TRUE))
  Dm_total   <- 1 - rmse_total
  if (clamp.zero) Dm_total <- max(0, Dm_total)

  # --- NUEVO: IC para Dm.total cuando zboot
  Dm_total_lwr <- Dm_total_upr <- NA_real_
  if (ci.method == "zboot") {
    q_hi <- (1 + conf.level) / 2
    # simulación conjunta de toda la muestra (k=1 por ítem)
    zhat_vec <- zfun(r)
    se <- 1 / sqrt(n.eff - 3)
    if (is.finite(se) && se > 0) {
      Dm_total_star <- numeric(B)
      for (b in seq_len(B)) {
        zsim_i <- rnorm(n, mean = zhat_vec, sd = se)
        if (metric == "fisher") {
          Xsim <- zsim_i
          Tsim <- zfun(ref.vec)
        } else {
          Xsim <- tanh(zsim_i)
          Tsim <- ref.vec
        }
        rmse_b <- sqrt(mean((Xsim - Tsim)^2, na.rm = TRUE))
        Dm_b   <- 1 - rmse_b
        if (clamp.zero) Dm_b <- max(0, Dm_b)
        Dm_total_star[b] <- Dm_b
      }
      qsT <- quantile(Dm_total_star, c(1 - q_hi, q_hi), na.rm = TRUE, names = FALSE)
      Dm_total_lwr <- qsT[1]; Dm_total_upr <- qsT[2]
    }
  }

  total_row <- data.frame(
    Item = "Dm.total",
    Dm   = round(Dm_total, round.digits),
    Dm.lwr = if (is.na(Dm_total_lwr)) NA_real_ else round(Dm_total_lwr, round.digits),
    Dm.upr = if (is.na(Dm_total_upr)) NA_real_ else round(Dm_total_upr, round.digits),
    coef.mean = round(mean(r, na.rm = TRUE), round.digits),
    abs.diff.mean = if (metric == "fisher") {
      round(mean(abs(zfun(r) - zfun(ref.vec)), na.rm = TRUE), round.digits)
    } else {
      round(mean(abs(r - ref.vec), na.rm = TRUE), round.digits)
    },
    stringsAsFactors = FALSE
  )
  out <- rbind(out, total_row)
  ## -------------------------------

  attr(out, "ref.value")   <- ref.value
  attr(out, "coef.col")    <- coef.col
  attr(out, "metric")      <- metric
  attr(out, "ci.method")   <- ci.method
  attr(out, "conf.level")  <- conf.level
  attr(out, "B")           <- B
  attr(out, "n.eff")       <- n.eff
  attr(out, "Dm.total")    <- unname(total_row$Dm)
  attr(out, "Dm.total.lwr")<- if (!is.na(total_row$Dm.lwr)) unname(total_row$Dm.lwr) else NA_real_
  attr(out, "Dm.total.upr")<- if (!is.na(total_row$Dm.upr)) unname(total_row$Dm.upr) else NA_real_

  out
}
