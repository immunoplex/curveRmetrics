# Parameter normalization layer
# All parameters converted to natural-scale (concentration domain) with
# unified family names before inflection point calculation.
# =============================================================================

# -----------------------------------------------------------------------------
# .canonical_family
# -----------------------------------------------------------------------------

#' Resolve a raw model name to a canonical family string
#'
#' Maps both Frequentist and Bayesian family
#' codes to a single controlled vocabulary used
#' throughout this package.
#'
#' @param raw_name Character scalar. The raw `model_name` or `curve_family`
#'   value as stored in the database.
#'
#' @return Character scalar. One of `"logistic4"`, `"logistic5"`,
#'   `"loglogistic4"`, `"loglogistic5"`, or `"gompertz4"`.
#'
#' @keywords internal
#' @noRd
.canonical_family <- function(raw_name) {
  lookup <- c(
    # NLS names  -> canonical
    "logistic4"   = "logistic4",
    "logistic5"   = "logistic5",
    "loglogistic4" = "loglogistic4",
    "loglogistic5" = "loglogistic5",
    "gompertz4"   = "gompertz4",
    # Bayes family codes  -> canonical
    "4pl"         = "logistic4",
    "5pl"         = "logistic5",
    "gompertz"    = "gompertz4"
    # Note: Bayes does not have separate loglogistic families;
    # "4pl" and "5pl" cover both symmetric and asymmetric logistic forms.
    # If your Stan model later adds loglogistic variants, extend here.
  )
  canonical <- lookup[raw_name]
  if (is.na(canonical)) {
    stop(sprintf("[.normalise_params] Unknown family/model name: '%s'", raw_name),
         call. = FALSE)
  }
  unname(canonical)
}

# -----------------------------------------------------------------------------
# .normalise_params
# -----------------------------------------------------------------------------
#' Normalise a database row to a common parameter representation
#'
#' Converts a single-row data frame or named list from `freq_curves` or
#' `bayes_curves` into a unified named list used by all downstream functions.
#'
#' Scale conventions after normalisation:
#' \itemize{
#'   \item \strong{Frequentist} — `c` is on the log10-concentration scale, as
#'     stored in the database; `b` is the scale parameter in log10 units.
#'   \item \strong{Bayesian} — `c` is converted from natural-scale EC50
#'     (concentration units) to log10 via `log10(c)`; `b` is rescaled from
#'     the natural-log domain stored by Stan to log10 via `b / log(10)`.
#'     Both `c` and `x` are therefore on the same log10 scale when passed to
#'     .evaluate_curve()] or .invert_curve().
#' }
#'
#' @param row  Single-row data frame or named list containing at minimum the
#'   columns/elements `model_name`, `a`, `b`, `c`, `d`, `g`.
#' @param method Character scalar. Either `"frequentist"` or `"bayesian"`.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{`family`}{Canonical family string from .canonical_family().}
#'     \item{`a`}{Lower asymptote (response scale).}
#'     \item{`b`}{Scale/slope parameter, log10-concentration domain.}
#'     \item{`c`}{Inflection-point concentration, log10 scale.}
#'     \item{`d`}{Upper asymptote (response scale).}
#'     \item{`g`}{Asymmetry parameter; 1 for symmetric 4-parameter models.}
#'     \item{`source`}{Character, echoes `method`.}
#'   }
#'
#' @keywords internal
.normalise_params <- function(row, method) {
list(
  family = .canonical_family(row$model_name),
  a      = as.numeric(row$a),
  b      = as.numeric(row$b),
  # Frequentist: c is already on the natural concentration scale
  # Bayesian:    c is stored as log_c in the Stan draws; bayes_curves
  #              saves the raw parameter so exp() is required here.
  #              If your pipeline already exponentiates before saving, remove exp().
  c      =  as.numeric(row$c), #if (method == "bayesian") exp(as.numeric(row$c)) else as.numeric(row$c),
  d      = as.numeric(row$d),
  g      = ifelse(is.na(row$g), 1, as.numeric(row$g)), #as.numeric(row$g %||% 1),
  source = method
)
}

# -----------------------------------------------------------------------------
# .pivot_curve_row
# -----------------------------------------------------------------------------

#' Pivot long-format parameter CI rows to a wide named list
#'
#' Takes the subset of `param_ci_df` rows belonging to a single curve and
#' reshapes them into the wide format expected by .normalise_params() and
#' [generate_lods()].
#'
#' @param curve_rows Data frame. Rows from a long-format parameter CI tibble
#'   (e.g. `freq_param_ci_df` or `bayes_param_ci_df`) for a single
#'   `curve_id`, containing columns `parameter`, `estimate`, `conf_lower`,
#'   `conf_upper`, `curve_id`, `method`, and `model_name`.
#'
#' @return Named list with elements `curve_id`, `method`, `model_name`,
#'   `a`, `b`, `c`, `d`, `g` (point estimates), `ci_lower` (named numeric
#'   vector of lower confidence bounds), and `ci_upper` (named numeric vector
#'   of upper confidence bounds).
#'
#' @keywords internal
#' @noRd
.pivot_curve_row <- function(curve_rows) {

  param_rows <- curve_rows[curve_rows$parameter %in% c("a", "b", "c", "d", "g"), ]

  est  <- setNames(param_rows$estimate,   param_rows$parameter)
  low  <- setNames(param_rows$conf_lower, param_rows$parameter)
  high <- setNames(param_rows$conf_upper, param_rows$parameter)

  list(
    curve_id   = curve_rows$curve_id[[1]],
    method     = curve_rows$method[[1]],
    model_name = curve_rows$model_name[[1]],
    a          = as.numeric(est["a"]),
    b          = as.numeric(est["b"]),
    c          = as.numeric(est["c"]),
    d          = as.numeric(est["d"]),
    g          = if ("g" %in% names(est)) as.numeric(est["g"]) else NA_real_,
    ci_lower   = low,
    ci_upper   = high
  )
}

#' Evaluate the forward curve at a given concentration
#'
#' All parameters and x are on the natural scale.
#' Matches functional forms in model_functions.R exactly.
#'
#' @param p Named list with family, a, b, c, d, g
#' @param x Numeric scalar or vector, natural concentration
#' @return Numeric vector of predicted response values
#' @keywords internal
.evaluate_curve <- function(p, x) {
  a <- p$a; b <- p$b; c <- p$c; d <- p$d; g <- p$g
  switch(p$family,
         "logistic4"    = d + (a - d) / (1 + exp((x - c) / b)),
         "logistic5"    = d + (a - d) / (1 + exp((x - c) / b))^g,
         "loglogistic4" = a + (d - a) / (1 + (x / c)^b),
         "loglogistic5" = a + (d - a) * (1 + g * exp(-b * (x - c)))^(-1 / g),
         "gompertz4"    = a + (d - a) * exp(-exp(-b * (x - c))),
         stop(sprintf("[.evaluate_curve] Unsupported family: '%s'", p$family), call. = FALSE)
  )
}

#' Back-calculate concentration from a response value
#'
#' Analytically inverts the forward dose-response curve to find the
#' concentration \eqn{x} corresponding to a given response \eqn{y}.
#' Returns a value on the same scale as \code{p$c} (log10-concentration
#' after normalisation by \code{.normalise_params()}).
#'
#' Use \code{.safe_invert()} rather than calling this function directly;
#' \code{.safe_invert()} clamps \eqn{y} to within the curve's
#' \eqn{(a, d)} range before inversion, preventing the logarithmic
#' blow-up that occurs when \eqn{y} is at or beyond an asymptote.
#'
#' @details
#' Inversion formulae by family (all solved analytically from the
#' corresponding forward model):
#' \describe{
#'   \item{\code{logistic4}}{\eqn{x = c + b \ln\!\left(\frac{a-d}{y-d} - 1\right)}}
#'   \item{\code{logistic5}}{\eqn{x = c + b \ln\!\left(\left(\frac{a-d}{y-d}\right)^{1/g} - 1\right)}}
#'   \item{\code{loglogistic4}}{\eqn{x = c \left(\frac{d-a}{y-a} - 1\right)^{1/b}}}
#'   \item{\code{loglogistic5}}{\eqn{x = c - \frac{\ln\!\left(\left(\frac{y-a}{d-a}\right)^{-g} - 1\right) - \ln g}{b}}}
#'   \item{\code{gompertz4}}{\eqn{x = c - \frac{\ln\!\left(-\ln\frac{y-a}{d-a}\right)}{b}}}
#' }
#' All formulae require \eqn{y} to lie strictly within \eqn{(\min(a,d),\, \max(a,d))};
#' values at or outside the asymptotes produce \code{NaN} or \code{Inf}.
#'
#' @param p Named list as returned by \code{.normalise_params()}, containing
#'   \code{family}, \code{a}, \code{b}, \code{c}, \code{d}, and \code{g}.
#' @param y Numeric scalar. Response value to invert; must lie strictly
#'   within the interval \eqn{(\min(a, d),\, \max(a, d))}.
#'
#' @return Numeric scalar. Concentration on the same scale as \code{p$c},
#'   or \code{NaN} / \code{Inf} if \eqn{y} is outside the asymptote range.
#'
#' @keywords internal
#' @noRd
.invert_curve <- function(p, y) {
  a <- p$a; b <- p$b; c <- p$c; d <- p$d; g <- p$g

  switch(p$family,

         "logistic4" = {
           c + b * log((a - d) / (y - d) - 1)
         },

         "logistic5" = {
           c + b * log(((a - d) / (y - d))^(1 / g) - 1)
         },

         "loglogistic4" = {
           c * ((d - a) / (y - a) - 1)^(1 / b)
         },

         "loglogistic5" = {
           c - (log(((y - a) / (d - a))^(-g) - 1) - log(g)) / b
         },

         "gompertz4" = {
           c - log(-log((y - a) / (d - a))) / b
         },

         stop(sprintf("[.invert_curve] Unsupported family: '%s'", p$family), call. = FALSE)
  )
}


#' Attach computed quality metrics to a curves data frame
#'
#' Left-joins four quality metric tables — sensitivity, limits of detection,
#' detection/reliable detection limits, and limits of quantification — onto
#' a curves data frame by \code{curve_id}.  All merges are left joins so that
#' curves with no matching metrics retain \code{NA} rather than being dropped.
#'
#' @param curves_df  Data frame of fitted curves, typically \code{freq_curves_df}
#'   or \code{bayes_curves_df}.  Must contain a \code{curve_id} column.
#' @param lods        Tibble returned by \code{generate_lods()}, with columns
#'   \code{curve_id}, \code{llod}, and \code{ulod}.
#' @param rdls        Tibble returned by \code{compute_mdc_rdl()}, with columns
#'   \code{curve_id}, \code{mindc}, \code{maxdc}, \code{minrdl}, and
#'   \code{maxrdl}.
#' @param sensitivity Tibble returned by \code{compute_assay_sensitivity()},
#'   with columns \code{curve_id} and \code{dydx_inflect}.
#' @param loqs        Tibble returned by \code{compute_loqs()}, with columns
#'   \code{curve_id}, \code{lloq}, \code{uloq}, \code{lloq_y}, and
#'   \code{uloq_y}.  The \code{method} column is dropped before joining to
#'   avoid a collision with the \code{method} column already present in
#'   \code{curves_df}.
#'
#' @return \code{curves_df} with the following columns appended:
#'   \describe{
#'     \item{\code{dydx_inflect}}{First derivative of response at the
#'       inflection point (assay sensitivity).}
#'     \item{\code{llod}, \code{ulod}}{Lower and upper limits of detection
#'       on the response scale.}
#'     \item{\code{mindc}, \code{maxdc}}{Lower and upper minimum detectable
#'       concentrations.}
#'     \item{\code{minrdl}, \code{maxrdl}}{Lower and upper reliable detection
#'       limits.}
#'     \item{\code{lloq}, \code{uloq}}{Lower and upper limits of
#'       quantification (concentration scale).}
#'     \item{\code{lloq_y}, \code{uloq_y}}{Predicted response values at the
#'       LOQ concentrations.}
#'   }
#'
#' @export
attach_quality_metrics <- function(curves_df, lods, rdls, sensitivity, loqs) {
  curves_df <- merge(
    curves_df,
    sensitivity[, c("curve_id", "dydx_inflect")],
    by = "curve_id",
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    lods,
    by = "curve_id",
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    rdls,
    by = "curve_id",
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    loqs[names(loqs) != "method"],
    by = "curve_id",
    all.x = TRUE
  )

  method <- unique(curves_df$method)

  # if (method == "bayesian") {
  #
  #   cols_to_log <- c(
  #     "inflect_x", "inflect_y", "dydx_inflect",
  #     "llod", "ulod",
  #     "mindc", "maxdc",
  #     "minrdl", "maxrdl",
  #     "lloq", "uloq",
  #     "lloq_y", "uloq_y"
  #   )
  #
  #   curves_df[cols_to_log] <- lapply(
  #     curves_df[cols_to_log],
  #     log10
  #   )
  # }

  return(curves_df)
}
