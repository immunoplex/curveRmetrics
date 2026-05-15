# =============================================================================
# Parameter normalisation layer
# All parameters converted to natural-scale (concentration domain) with
# unified family names before inflection point calculation.
# =============================================================================


#' Resolve a raw model/family name to its canonical form
#'
#' Maps NLS model names (from `model_functions.R`) and Bayesian family codes
#' (from `curve_families.R`) to a single canonical string used throughout this
#' module.
#'
#' @param raw_name A character string; the model or family name to resolve.
#'
#' @return A single canonical family name (character).
#' @keywords internal
.canonical_family <- function(raw_name) {
  lookup <- c(
    "logistic4"    = "logistic4",
    "logistic5"    = "logistic5",
    "loglogistic4" = "loglogistic4",
    "loglogistic5" = "loglogistic5",
    "gompertz4"    = "gompertz4",
    "4pl"          = "logistic4",
    "5pl"          = "logistic5",
    "gompertz"     = "gompertz4"
  )
  canonical <- lookup[raw_name]
  if (is.na(canonical)) {
    stop(sprintf("[.normalise_params] Unknown family/model name: '%s'", raw_name),
         call. = FALSE)
  }
  unname(canonical)
}


#' Transform curve parameters to natural (untransformed) units
#'
#' Converts raw fitted parameters stored in a data frame to their natural-scale
#' equivalents, appending `*_nat` columns.  The transformation rules are:
#'
#' \describe{
#'   \item{Frequentist}{Model was fit on \eqn{\log_{10}(x)}, so
#'     \code{c} (EC50) is back-transformed via \eqn{10^c}, and
#'     \code{b} (Hill slope) is multiplied by \eqn{\ln(10)} to convert from the
#'     \eqn{\log_{10}} scale to the natural-log scale. Asymptotes \code{a} and
#'     \code{d} are also back-transformed via \eqn{10^a}, \eqn{10^d}.}
#'   \item{Bayesian}{Model is parameterised in \eqn{\ln(x)}, with \code{c}
#'     already representing the median effective concentration; \code{a} and
#'     \code{d} are already on the response scale.  \code{c_nat} is set to
#'     \eqn{\ln(c)} and \code{b_nat} to \eqn{1/b} to match the internal Stan
#'     parameterisation.}
#' }
#'
#' The asymmetry parameter \code{g} is dimensionless and unaffected by any
#' scale transformation; it is set to 1 for symmetric (4-parameter) families.
#'
#' @param df A data frame with at minimum the columns \code{curve_id},
#'   \code{method}, \code{model_name}, \code{a}, \code{b}, \code{c}, \code{d},
#'   and \code{g}.
#'
#' @return \code{df} with additional columns \code{a_nat}, \code{b_nat},
#'   \code{c_nat}, \code{d_nat}, and \code{g_nat} on the natural-unit scale.
#'
#' @importFrom dplyr mutate case_when
#' @export
transform_to_natural_units <- function(df) {

  LN10 <- log(10)   # approximately 2.302585 - conversion factor between log bases

  df$model_name <- vapply(
    df$model_name,
    .canonical_family,
    character(1)
  )

  df %>%
    mutate(
      c_nat = case_when(
        method == "frequentist" ~ 10^c,
        method == "bayesian"    ~ log(c),
        TRUE                    ~ NA_real_
      ),
      b_nat = case_when(
        method == "frequentist" ~ b * LN10,
        method == "bayesian"    ~ 1 / b,
        TRUE                    ~ NA_real_
      ),
      a_nat = case_when(
        method == "frequentist" ~ 10^a,
        method == "bayesian"    ~ a
      ),
      d_nat = case_when(
        method == "frequentist" ~ 10^d,
        method == "bayesian"    ~ d
      ),
      g_nat = ifelse(
        model_name %in% c("logistic4", "loglogistic4", "gompertz4"),
        1,
        df$g
      )
    )
}


#' Normalise a single DB row to a natural-scale parameter list
#'
#' Prefers pre-computed \code{*_nat} columns produced by
#' \code{\link{transform_to_natural_units}} when they are present on \code{row}
#' (i.e. when called with \code{curves_nat}).  Falls back to raw DB columns
#' plus inline conversion otherwise.
#'
#' Either way the returned list has an identical structure and scale:
#' \describe{
#'   \item{a, d}{Natural-scale response asymptotes.}
#'   \item{b}{Hill slope on the \eqn{\ln(x)} scale.}
#'   \item{c}{Natural-scale EC50 (concentration units).}
#'   \item{g}{Dimensionless asymmetry/shape parameter.}
#' }
#'
#' @param row A single-row data frame (or named list) with curve parameters.
#' @param method Character; either \code{"frequentist"} or \code{"bayesian"}.
#' @param verbose Logical; if \code{TRUE} (default \code{FALSE}) emits a
#'   message when natural-scale columns are detected and used.
#'
#' @return A named list with elements \code{family}, \code{a}, \code{b},
#'   \code{c}, \code{d}, \code{g}, and \code{source}.
#' @keywords internal
.normalise_params <- function(row, method, verbose = FALSE) {

  family  <- .canonical_family(row$model_name)
  has_nat <- all(c("a_nat", "b_nat", "c_nat", "d_nat", "g_nat") %in% names(row))

  if (has_nat) {
    if (verbose) message("Using natural scale parameters...")
    return(list(
      family = family,
      a      = as.numeric(row$a_nat),
      b      = as.numeric(row$b_nat),
      c      = as.numeric(row$c_nat),
      d      = as.numeric(row$d_nat),
      g      = as.numeric(row$g_nat),
      source = method
    ))
  }

  LN10  <- log(10)
  a_raw <- as.numeric(row$a)
  b_raw <- as.numeric(row$b)
  c_raw <- as.numeric(row$c)
  d_raw <- as.numeric(row$d)
  g_raw <- ifelse(is.na(row$g), 1, as.numeric(row$g))

  if (method == "frequentist") {
    list(
      family = family,
      a      = 10^a_raw,
      b      = b_raw * LN10,
      c      = 10^c_raw,
      d      = 10^d_raw,
      g      = g_raw,
      source = method
    )
  } else {
    list(
      family = family,
      a      = a_raw,
      b      = b_raw,
      c      = c_raw,
      d      = d_raw,
      g      = g_raw,
      source = method
    )
  }
}


#' Pivot a long-format parameter CI data frame to a wide named list
#'
#' Selects rows for parameters \code{a}, \code{b}, \code{c}, \code{d}, and
#' (optionally) \code{g}, then pivots them into named vectors of point
#' estimates and confidence bounds.
#'
#' When the column \code{estimate_nat} is present on \code{curve_rows} the
#' natural-scale estimates and CI bounds are used; otherwise the raw DB
#' columns \code{estimate}, \code{conf_lower}, and \code{conf_upper} are used.
#'
#' @param curve_rows A data frame of CI rows for a single curve, containing at
#'   minimum \code{parameter}, \code{estimate}, \code{conf_lower},
#'   \code{conf_upper}, \code{curve_id}, \code{method}, and \code{model_name}.
#'   May additionally contain \code{estimate_nat}, \code{conf_lower_nat}, and
#'   \code{conf_upper_nat}.
#' @param verbose Logical; if \code{TRUE} (default \code{FALSE}) emits a
#'   message indicating whether natural-scale columns were found.
#'
#' @return A named list with elements \code{curve_id}, \code{method},
#'   \code{model_name}, \code{a}, \code{b}, \code{c}, \code{d}, \code{g},
#'   \code{ci_lower} (named numeric vector), and \code{ci_upper} (named numeric
#'   vector).
#' @keywords internal
.pivot_curve_row <- function(curve_rows, verbose = FALSE) {

  param_rows <- curve_rows[curve_rows$parameter %in% c("a", "b", "c", "d", "g"), ]
  has_nat    <- "estimate_nat" %in% names(param_rows)

  if (verbose) {
    if (has_nat) {
      message("using natural scale parameters...")
    } else {
      message("not using natural parameters...")
    }
  }

  est  <- setNames(if (has_nat) param_rows$estimate_nat   else param_rows$estimate,   param_rows$parameter)
  low  <- setNames(if (has_nat) param_rows$conf_lower_nat else param_rows$conf_lower, param_rows$parameter)
  high <- setNames(if (has_nat) param_rows$conf_upper_nat else param_rows$conf_upper, param_rows$parameter)

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


#' Evaluate a dose-response curve at a given concentration
#'
#' All parameters and \code{x} are expected on the \strong{natural scale}
#' (\code{c} = natural EC50, \code{x} = natural concentration).  The
#' \eqn{\log(x/c)/b} substitution is used throughout so that the formulas
#' accept raw concentration directly.
#'
#' Supported families:
#' \describe{
#'   \item{logistic4}{\eqn{y = d + (a-d) / (1 + \exp((x-c)/b))}}
#'   \item{logistic5}{\eqn{y = d + (a-d) / (1 + \exp((x-c)/b))^g}}
#'   \item{loglogistic4}{\eqn{y = a + (d-a) / (1 + (x/c)^b)}}
#'   \item{loglogistic5}{\eqn{y = a + (d-a)(1 + g\exp(-b(x-c)))^{-1/g}}}
#'   \item{gompertz4}{\eqn{y = a + (d-a)\exp(-\exp(-b(x-c)))}}
#' }
#'
#' @param p A named list of curve parameters with elements \code{family},
#'   \code{a}, \code{b}, \code{c}, \code{d}, and \code{g}.
#' @param x Numeric; concentration value(s) at which to evaluate the curve.
#'
#' @return Numeric; predicted response value(s).
#' @keywords internal
.evaluate_curve <- function(p, x) {
  a <- p$a; b <- p$b; c <- p$c; d <- p$d; g <- p$g
  switch(p$family,
         "logistic4"    = d + (a - d) / (1 + exp((x - c) / b)),
         "logistic5"    = d + (a - d) / (1 + exp((x - c) / b))^g,
         "loglogistic4" = a + (d - a) / (1 + (x / c)^b),
         "loglogistic5" = a + (d - a) * (1 + g * exp(-b * (x - c)))^(-1 / g),
         "gompertz4"    = a + (d - a) * exp(-exp(-b * (x - c))),
         stop(sprintf("[.evaluate_curve] Unsupported family: '%s'", p$family))
  )
}


#' Evaluate a dose-response curve using method-appropriate scaling
#'
#' A thin wrapper around \code{\link{.evaluate_curve}} that applies the
#' correct log transformation depending on whether the model was fit by the
#' frequentist or Bayesian method.
#'
#' \describe{
#'   \item{Frequentist}{Parameters are on the \eqn{\log_{10}} scale.
#'     \code{x} is converted to \eqn{\log_{10}(x)}, the curve is evaluated,
#'     and the result is back-transformed via \eqn{10^y}.}
#'   \item{Bayesian}{Parameters and \code{x} are passed as \eqn{\ln(x)} to
#'     match the Stan model's internal parameterisation.}
#' }
#'
#' @param p A named list of curve parameters (see \code{\link{.evaluate_curve}}).
#' @param x Numeric; natural-scale concentration value(s).
#' @param method Character; \code{"frequentist"} or \code{"bayesian"}.
#'
#' @return Numeric; predicted response on the natural scale.
#' @keywords internal
.evaluate_curve_auto <- function(p, x, method) {
  if (method == "frequentist") {
    x_native <- log10(x)
    y_log10  <- .evaluate_curve(p, x_native)
    return(10^y_log10)
  } else {
    return(.evaluate_curve(p, log(x)))
  }
}


#' Compute the inflection point concentration for a dose-response curve
#'
#' Finds the natural-scale concentration \eqn{x^*} at which the second
#' derivative of the curve is zero (i.e. the point of maximum slope).  All
#' returned values are natural concentrations.
#'
#' Derivations use the substitution \eqn{t = \ln(x/c)}:
#' \describe{
#'   \item{logistic4}{Symmetric \eqn{\Rightarrow x^* = c}.}
#'   \item{logistic5}{\eqn{z^* = (b-1)/(bg+1)},
#'     \eqn{x^* = c \cdot z^{*1/b}}.}
#'   \item{loglogistic4}{Symmetric \eqn{\Rightarrow x^* = c}.}
#'   \item{loglogistic5}{\eqn{u^* = (b-1)/(b+1)},
#'     \eqn{x^* = c \exp((\ln g - \ln u^*)/b)}.}
#'   \item{gompertz4}{\eqn{u^* = (b-1)/b},
#'     \eqn{x^* = c \exp(-\ln u^*/b)}.}
#' }
#'
#' @param p A named list of natural-scale curve parameters (see
#'   \code{\link{.evaluate_curve}}).
#'
#' @return Numeric scalar; natural-scale concentration at the inflection point.
#'   Falls back to \code{c} with a warning when the analytical condition
#'   cannot be satisfied.
#' @keywords internal
.inflection_x <- function(p) {

  b <- p$b
  g <- p$g
  c <- p$c

  switch(p$family,

         "logistic4" = {
           c
         },

         "logistic5" = {
           if (b <= 1) warning(sprintf(
             "[.inflection_x] logistic5: b = %.4f <= 1; inflection point may be unreliable.", b
           ))
           z_star <- (b - 1) / (b * g + 1)
           if (!is.finite(z_star) || z_star <= 0) {
             warning("[.inflection_x] logistic5: z_star <= 0; falling back to c.")
             return(c)
           }
           c * z_star^(1 / b)
         },

         "loglogistic4" = {
           c
         },

         "loglogistic5" = {
           if (b <= 1) warning(sprintf(
             "[.inflection_x] loglogistic5: b = %.4f <= 1; inflection point may be unreliable.", b
           ))
           u_star <- (b - 1) / (b + 1)
           if (!is.finite(u_star) || u_star <= 0) {
             warning("[.inflection_x] loglogistic5: u_star <= 0; falling back to c.")
             return(c)
           }
           c * exp((log(g) - log(u_star)) / b)
         },

         "gompertz4" = {
           if (b <= 1) {
             warning(sprintf(
               "[.inflection_x] gompertz4: b = %.4f <= 1; falling back to c.", b
             ))
             return(c)
           }
           u_star <- (b - 1) / b
           c * exp(-log(u_star) / b)
         },

         stop(sprintf("[.inflection_x] Unsupported family: '%s'", p$family), call. = FALSE)
  )
}


#' Compute the analytical inflection point as the response-range midpoint
#'
#' Finds the natural-scale concentration \eqn{x^*} at which the curve equals
#' the geometric midpoint of its response range, i.e.
#' \eqn{y^* = \sqrt{a \cdot d}} (Bayesian) or the equivalent arithmetic
#' midpoint in \eqn{\log_{10}} space (frequentist).  This definition is
#' unified across both fitting methods and requires no root-finding.
#'
#' @param p A named list of curve parameters.  For frequentist curves the
#'   parameters are on the raw \eqn{\log_{10}} scale; for Bayesian curves they
#'   are on the natural scale (with \code{c = ln(c_raw)}).
#' @param method Character; \code{"frequentist"} or \code{"bayesian"}.
#'
#' @return Numeric scalar; natural-scale concentration at the midpoint
#'   inflection.  Falls back to the back-transformed EC50 if the analytical
#'   condition cannot be evaluated.
#' @keywords internal
.inflect_x_analytical <- function(p, method) {

  a <- p$a; b <- p$b; c <- p$c; d <- p$d; g <- p$g

  if (method == "frequentist") {

    switch(p$family,
           "logistic4" = {
             10^c
           },
           "logistic5" = {
             inner <- 2^(1 / g) - 1
             if (!is.finite(inner) || inner <= 0) return(10^c)
             10^(c + b * log(inner))
           },
           "loglogistic4" = {
             10^c
           },
           "loglogistic5" = {
             denom <- 2^g - 1
             if (!is.finite(denom) || denom <= 0) return(10^c)
             t_star <- log(g / denom) / b
             10^c * exp(t_star)
           },
           "gompertz4" = {
             10^(c - log(log(2)) / b)
           },
           stop(sprintf("[.inflect_x_analytical] Unsupported family: '%s'", p$family))
    )

  } else {

    y_target <- sqrt(a * d)

    switch(p$family,
           "logistic4" = ,
           "logistic5" = {
             ratio <- (a - d) / (y_target - d)
             inner <- ratio^(1 / g) - 1
             if (!is.finite(inner) || inner <= 0) return(exp(c))
             exp(c) * inner^b
           },
           "loglogistic4" = {
             ratio <- (d - y_target) / (y_target - a)
             if (!is.finite(ratio) || ratio <= 0) return(exp(c))
             exp(c) * ratio^(1 / b)
           },
           "loglogistic5" = {
             ratio  <- (y_target - a) / (d - a)
             inner  <- ratio^(-g) - 1
             if (!is.finite(inner) || inner <= 0) return(exp(c))
             t_star <- -log(inner / g) / b
             exp(c + t_star)
           },
           "gompertz4" = {
             ratio <- (y_target - a) / (d - a)
             if (!is.finite(ratio) || ratio <= 0 || ratio >= 1) return(exp(c))
             t_star <- -log(-log(ratio)) / b
             exp(c + t_star)
           },
           stop(sprintf("[.inflect_x_analytical] Unsupported family: '%s'", p$family))
    )
  }
}


#' Compute the inflection point for a single curve row (worker)
#'
#' Normalises parameters via \code{\link{.normalise_params}}, locates the
#' inflection concentration with \code{\link{.inflection_x}}, and evaluates the
#' corresponding response with \code{\link{.evaluate_curve}}.
#'
#' @param row A single-row data frame with curve parameters and a \code{method}
#'   column.
#' @param verbose Logical; if \code{TRUE} (default) emits diagnostic messages
#'   with the resolved parameters and computed inflection coordinates.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{inflect_x}{Natural-scale concentration at the inflection point.}
#'   \item{inflect_y}{Predicted response at the inflection point.}
#'   \item{family}{Canonical family name.}
#'   \item{source}{Fitting method (\code{"frequentist"} or \code{"bayesian"}).}
#'   \item{params}{The normalised parameter list passed to the worker functions.}
#' }
#' @export
compute_inflection_point_worker <- function(row, verbose = TRUE) {

  method <- row$method
  p      <- .normalise_params(row, method)

  if (verbose)
    message(sprintf(
      "[compute_inflection_point] family=%s source=%s a=%.4f b=%.4f c=%.4f d=%.4f g=%.4f",
      p$family, p$source, p$a, p$b, p$c, p$d, p$g
    ))

  inflect_x <- tryCatch(
    .inflection_x(p),
    error = function(e) { warning(e$message); NA_real_ }
  )

  inflect_y <- tryCatch(
    .evaluate_curve(p, inflect_x),
    error = function(e) { warning(e$message); NA_real_ }
  )

  if (verbose)
    message(sprintf(
      "[compute_inflection_point] inflect_x=%.6f inflect_y=%.6f",
      inflect_x, inflect_y
    ))

  list(
    inflect_x = as.numeric(inflect_x),
    inflect_y = as.numeric(inflect_y),
    family    = p$family,
    source    = p$source,
    params    = p
  )
}


#' Compute inflection points for all curves in a data frame
#'
#' Iterates over unique \code{(curve_id, method)} combinations, resolves
#' parameters via \code{\link{.inflect_x_analytical}}, evaluates the response at
#' that concentration, and joins the results back onto the input data frame.
#'
#' The inflection point is defined as the concentration at the midpoint of the
#' response range (\eqn{y^* = \sqrt{a \cdot d}}), computed analytically without
#' root-finding.
#'
#' @param curves_df A data frame with one row per curve containing at minimum
#'   \code{curve_id}, \code{method}, \code{model_name}, \code{a}, \code{b},
#'   \code{c}, \code{d}, and \code{g}.  Frequentist rows hold raw
#'   \eqn{\log_{10}} parameters; Bayesian rows must also carry \code{*_nat}
#'   columns.
#' @param verbose Logical; if \code{TRUE} (default) emits a message per curve
#'   with the computed inflection coordinates.
#'
#' @return \code{curves_df} with two new (or replaced) columns:
#'   \code{inflect_x} and \code{inflect_y}, both on the natural scale.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows left_join select any_of
#' @importFrom tibble tibble
#' @export
compute_inflection_point <- function(curves_df, verbose = TRUE) {

  keys <- unique(curves_df[, c("curve_id", "method")])

  results <- purrr::map(seq_len(nrow(keys)), function(i) {

    cid    <- keys$curve_id[[i]]
    method <- keys$method[[i]]
    row    <- curves_df[curves_df$curve_id == cid & curves_df$method == method, ]
    family <- .canonical_family(row$model_name)

    p <- if (method == "frequentist") {
      list(family = family,
           a = as.numeric(row$a),
           b = as.numeric(row$b),
           c = as.numeric(row$c),
           d = as.numeric(row$d),
           g = if (is.na(row$g)) 1 else as.numeric(row$g))
    } else {
      list(family = family,
           a = as.numeric(row$a_nat),
           b = as.numeric(row$b_nat),
           c = as.numeric(row$c_nat),
           d = as.numeric(row$d_nat),
           g = as.numeric(row$g_nat))
    }

    inflect_x <- tryCatch(
      .inflect_x_analytical(p, method),
      error = function(e) {
        warning(sprintf("[compute_inflection_point] curve_id=%s (%s): %s",
                        cid, method, e$message))
        NA_real_
      }
    )

    inflect_y <- tryCatch({
      if (method == "frequentist")
        10^.evaluate_curve(p, log10(inflect_x))
      else
        .evaluate_curve(p, log(inflect_x))
    }, error = function(e) { warning(e$message); NA_real_ })

    if (verbose)
      message(sprintf(
        "[compute_inflection_point] curve_id=%s (%s)  inflect_x=%.4f  inflect_y=%.4f",
        cid, method, inflect_x, inflect_y
      ))

    tibble::tibble(
      curve_id  = cid, method    = method,
      inflect_x = as.numeric(inflect_x),
      inflect_y = as.numeric(inflect_y)
    )
  })

  dplyr::bind_rows(results) %>%
    { dplyr::left_join(
      dplyr::select(curves_df, -dplyr::any_of(c("inflect_x", "inflect_y"))),
      ., by = c("curve_id", "method")) }
}


#' Compute assay sensitivity at the inflection point
#'
#' Estimates \eqn{dy/dx} at \code{inflect_x} for each curve using a symmetric
#' central finite difference with step size \eqn{h = |x| \times 10^{-5}}.
#' Parameters are normalised via \code{\link{.normalise_params}} before
#' evaluation.
#'
#' @param curves_df A data frame with one row per curve.  Must contain
#'   \code{curve_id}, \code{method}, \code{model_name}, and \code{inflect_x}
#'   (natural-scale, as produced by \code{\link{compute_inflection_point}}),
#'   plus raw or \code{*_nat} parameter columns.
#' @param verbose Logical; if \code{TRUE} (default) emits a message per curve
#'   with \code{inflect_x} and the estimated slope.
#'
#' @return A data frame with columns \code{curve_id}, \code{method},
#'   \code{inflect_x}, and \code{dydx_inflect}.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @export
compute_assay_sensitivity <- function(curves_df, verbose = TRUE) {

  results <- purrr::map(seq_len(nrow(curves_df)), function(i) {

    row  <- curves_df[i, ]
    p    <- .normalise_params(row, row$method)
    ix   <- as.numeric(row$inflect_x)

    dydx <- tryCatch({
      if (!is.na(ix)) {
        h <- if (abs(ix) > 0) abs(ix) * 1e-5 else 1e-8
        (.evaluate_curve(p, ix + h) - .evaluate_curve(p, ix - h)) / (2 * h)
      } else {
        NA_real_
      }
    }, error = function(e) {
      if (verbose) message(sprintf("[compute_assay_sensitivity] curve_id=%s error: %s",
                                   row$curve_id, e$message))
      NA_real_
    })

    if (verbose)
      message(sprintf(
        "[compute_assay_sensitivity] curve_id=%s (%s)  inflect_x=%.4f  dydx=%.6f",
        row$curve_id, row$method, ix, dydx
      ))

    tibble::tibble(
      curve_id     = row$curve_id,
      method       = row$method,
      inflect_x    = ix,
      dydx_inflect = as.numeric(dydx)
    )
  })

  dplyr::bind_rows(results)
}


#' Derive lower and upper limits of detection from asymptote CI bounds
#'
#' Limits of detection are defined on the response scale:
#' \describe{
#'   \item{LLOD}{Upper CI bound of the lower asymptote \code{a}.}
#'   \item{ULOD}{Lower CI bound of the upper asymptote \code{d}.}
#' }
#'
#' When the ULOD is less than the LLOD (indicating an implausible CI) it is
#' set to \code{NA}.  Only rows flagged as \code{is_best_model} are processed.
#'
#' @param param_ci_df A long-format data frame of parameter confidence
#'   intervals with columns \code{curve_id}, \code{method}, \code{parameter},
#'   \code{estimate}, \code{conf_lower}, \code{conf_upper}, and
#'   \code{is_best_model}.  May optionally contain \code{*_nat} variants.
#' @param verbose Logical; if \code{TRUE} emits a message per curve
#'   with the computed LLOD and ULOD.
#'
#' @return A data frame with columns \code{curve_id}, \code{method},
#'   \code{llod}, and \code{ulod}.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @export
generate_lods <- function(param_ci_df, verbose = FALSE) {

  best_df <- param_ci_df[param_ci_df$is_best_model, ]
  keys    <- unique(best_df[, c("curve_id", "method")])

  results <- purrr::map(seq_len(nrow(keys)), function(i) {

    cid    <- keys$curve_id[[i]]
    method <- keys$method[[i]]

    rows <- best_df[best_df$curve_id == cid & best_df$method == method, ]
    row  <- .pivot_curve_row(rows)

    ulod <- as.numeric(row$ci_lower["d"])
    llod <- as.numeric(row$ci_upper["a"])

    if (is.na(ulod) || is.na(llod) || ulod < llod) {
      if (verbose)
        message(sprintf(
          "[generate_lods] curve_id=%s (%s)  ULOD invalid (%.4f vs LLOD %.4f) - ULOD set to NA",
          cid, method, ulod, llod
        ))
      ulod <- NA_real_
    }

    if (verbose)
      message(sprintf(
        "[generate_lods] curve_id=%s (%s)  LLOD=%.4f  ULOD=%s",
        cid, method, llod,
        if (is.na(ulod)) "NA" else sprintf("%.4f", ulod)
      ))

    tibble::tibble(curve_id = cid, method = method, llod = llod, ulod = ulod)
  })

  dplyr::bind_rows(results)
}


#' Back-calculate natural concentration from a response value
#'
#' Analytically inverts each supported dose-response family.  All parameters
#' and the returned concentration are on the natural scale.
#'
#' @param p A named list of natural-scale curve parameters (see
#'   \code{\link{.evaluate_curve}}).
#' @param y Numeric; response value(s) to invert.  Must lie within the range
#'   \eqn{[\min(a, d),\, \max(a, d)]}.
#'
#' @return Numeric; natural-scale concentration(s) corresponding to \code{y}.
#' @keywords internal
.invert_curve <- function(p, y) {
  a <- p$a; b <- p$b; c <- p$c; d <- p$d; g <- p$g

  switch(p$family,

         "logistic4" = {
           c * ((a - d) / (y - d) - 1)^b
         },

         "logistic5" = {
           c * (((a - d) / (y - d))^(1 / g) - 1)^b
         },

         "loglogistic4" = {
           c * ((d - a) / (y - a) - 1)^(1 / b)
         },

         "loglogistic5" = {
           c * exp((log(g) - log(((y - a) / (d - a))^(-g) - 1)) / b)
         },

         "gompertz4" = {
           c * exp(-log(-log((y - a) / (d - a))) / b)
         },

         stop(sprintf("[.invert_curve] Unsupported family: '%s'", p$family), call. = FALSE)
  )
}


#' Compute the Minimum Detectable Concentration and Reliable Detection Limits
#'
#' Derives MDC and RDL by inverting the dose-response curve at the LOD
#' response thresholds:
#' \describe{
#'   \item{MDC}{\code{mindc} = concentration at LLOD;
#'     \code{maxdc} = concentration at ULOD.}
#'   \item{RDL}{\code{minrdl} = concentration at LLOD on the lower CI bound
#'     curve (\code{d} replaced by \code{ci_lower["d"]});
#'     \code{maxrdl} = concentration at ULOD on the upper CI bound curve
#'     (\code{d} replaced by \code{ci_upper["d"]}).}
#' }
#'
#' All returned concentrations are on the natural scale.  Response values
#' outside \eqn{[\min(a, d),\, \max(a, d)]} are silently coerced to
#' \code{NA} before inversion.
#'
#' Only rows flagged as \code{is_best_model} are processed.
#'
#' @param param_ci_df A long-format CI data frame (same structure as required
#'   by \code{\link{generate_lods}}).
#' @param lods A data frame as returned by \code{\link{generate_lods}} with
#'   columns \code{curve_id}, \code{method}, \code{llod}, and \code{ulod}.
#' @param verbose Logical; if \code{TRUE} (default) emits a message per curve
#'   with the four computed concentration limits.
#'
#' @return A data frame with columns \code{curve_id}, \code{method},
#'   \code{mindc}, \code{maxdc}, \code{minrdl}, and \code{maxrdl}.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @export
compute_mdc_rdl <- function(param_ci_df, lods, verbose = TRUE) {

  best_df <- param_ci_df[param_ci_df$is_best_model, ]
  keys    <- unique(best_df[, c("curve_id", "method")])

  results <- purrr::map(seq_len(nrow(keys)), function(i) {

    cid    <- keys$curve_id[[i]]
    method <- keys$method[[i]]

    rows <- best_df[best_df$curve_id == cid & best_df$method == method, ]
    row  <- .pivot_curve_row(rows)

    p <- list(
      family = .canonical_family(row$model_name),
      a      = row$a,
      b      = row$b,
      c      = row$c,
      d      = row$d,
      g      = ifelse(is.na(row$g), 1, row$g),
      source = method
    )

    lod_row <- lods[lods$curve_id == cid & lods$method == method, ]
    llod    <- as.numeric(lod_row$llod)
    ulod    <- as.numeric(lod_row$ulod)

    y_lo <- min(p$a, p$d)
    y_hi <- max(p$a, p$d)

    safe_invert <- function(params, y) {
      if (is.na(y) || y < y_lo || y > y_hi) return(NA_real_)
      tryCatch(.invert_curve(params, y), error = function(e) NA_real_)
    }

    mindc <- safe_invert(p, llod)
    maxdc <- safe_invert(p, ulod)

    p_lo <- modifyList(p, list(d = as.numeric(row$ci_lower["d"])))
    p_hi <- modifyList(p, list(d = as.numeric(row$ci_upper["d"])))

    y_lo_rdl_min <- min(p_lo$a, p_lo$d)
    y_hi_rdl_min <- max(p_lo$a, p_lo$d)
    y_lo_rdl_max <- min(p_hi$a, p_hi$d)
    y_hi_rdl_max <- max(p_hi$a, p_hi$d)

    safe_invert_lo <- function(y) {
      if (is.na(y) || y < y_lo_rdl_min || y > y_hi_rdl_min) return(NA_real_)
      tryCatch(.invert_curve(p_lo, y), error = function(e) NA_real_)
    }

    safe_invert_hi <- function(y) {
      if (is.na(y) || y < y_lo_rdl_max || y > y_hi_rdl_max) return(NA_real_)
      tryCatch(.invert_curve(p_hi, y), error = function(e) NA_real_)
    }

    minrdl <- safe_invert_lo(llod)
    maxrdl <- safe_invert_hi(ulod)

    if (verbose)
      message(sprintf(
        "[compute_mdc_rdl] curve_id=%s (%s)  mindc=%s  maxdc=%s  minrdl=%s  maxrdl=%s",
        cid, method,
        format(mindc,  digits = 4), format(maxdc,  digits = 4),
        format(minrdl, digits = 4), format(maxrdl, digits = 4)
      ))

    tibble::tibble(
      curve_id = cid,
      method   = method,
      mindc    = as.numeric(mindc),
      maxdc    = as.numeric(maxdc),
      minrdl   = as.numeric(minrdl),
      maxrdl   = as.numeric(maxrdl)
    )
  })

  dplyr::bind_rows(results)
}


#' Compute limits of quantification from the second derivative
#'
#' Identifies the lower LOQ (\code{lloq}) and upper LOQ (\code{uloq}) as the
#' concentrations at the maximum and minimum of the second derivative of the
#' fitted curve, respectively.  Local extrema are found via sign changes in
#' the first differences of the pre-computed \code{d2x_y} series; sub-grid
#' vertex positions are refined with a direct 3-point parabola interpolation.
#'
#' The response at each LOQ concentration is evaluated using the
#' method-appropriate evaluator:
#' \itemize{
#'   \item Frequentist: \eqn{10^{f(\log_{10}(x))}}.
#'   \item Bayesian: \eqn{f(\ln(x))}.
#' }
#'
#' @param curves_df A data frame with one row per curve (same requirements as
#'   \code{\link{compute_inflection_point}}), containing both raw and
#'   \code{*_nat} parameter columns.
#' @param second_deriv_df A data frame of pre-computed second derivative values
#'   with columns \code{curve_id}, \code{method}, \code{concentration} (natural
#'   scale), and \code{d2x_y}.
#' @param verbose Logical; if \code{TRUE} (default) emits a message per curve
#'   with the computed LOQ concentrations and response values.
#'
#' @return A data frame with columns \code{curve_id}, \code{method},
#'   \code{lloq}, \code{uloq}, \code{lloq_y}, and \code{uloq_y}.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @export
compute_loqs <- function(curves_df, second_deriv_df, verbose = TRUE) {

  keys <- unique(curves_df[, c("curve_id", "method")])

  results <- purrr::map(seq_len(nrow(keys)), function(i) {

    cid    <- keys$curve_id[[i]]
    method <- keys$method[[i]]

    row    <- curves_df[curves_df$curve_id == cid & curves_df$method == method, ]
    family <- .canonical_family(row$model_name)

    p_raw <- list(family = family,
                  a = as.numeric(row$a), b = as.numeric(row$b),
                  c = as.numeric(row$c), d = as.numeric(row$d),
                  g = if (is.na(row$g)) 1 else as.numeric(row$g))

    p_nat <- list(family = family,
                  a = as.numeric(row$a_nat), b = as.numeric(row$b_nat),
                  c = as.numeric(row$c_nat), d = as.numeric(row$d_nat),
                  g = as.numeric(row$g_nat))

    eval_y <- if (method == "frequentist") {
      function(x) 10^.evaluate_curve(p_raw, log10(x))
    } else {
      function(x) .evaluate_curve(p_nat, log(x))
    }

    sd_sub <- second_deriv_df[second_deriv_df$curve_id == cid &
                                second_deriv_df$method == method, ]

    if (nrow(sd_sub) < 3) {
      if (verbose) message(sprintf(
        "[compute_loqs] curve_id=%s (%s): fewer than 3 second deriv points, skipping",
        cid, method))
      return(tibble::tibble(
        curve_id = cid, method = method,
        lloq = NA_real_, uloq = NA_real_,
        lloq_y = NA_real_, uloq_y = NA_real_
      ))
    }

    x <- as.numeric(sd_sub$concentration)
    y <- as.numeric(sd_sub$d2x_y)

    dy      <- diff(y)
    idx_max <- which(dy[-1] < 0 & dy[-length(dy)] > 0) + 1
    idx_min <- which(dy[-1] > 0 & dy[-length(dy)] < 0) + 1

    interpolate_vertex <- function(idx) {
      xi <- x[(idx - 1):(idx + 1)]
      yi <- y[(idx - 1):(idx + 1)]

      x1 <- xi[1]; x2 <- xi[2]; x3 <- xi[3]
      y1 <- yi[1]; y2 <- yi[2]; y3 <- yi[3]

      denom <- (x1 - x2) * (x1 - x3) * (x2 - x3)
      if (abs(denom) < .Machine$double.eps * 1e8)
        return(list(x = x2, y = y2))

      a_coef <- (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
      b_coef <- (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3)) / denom

      if (abs(a_coef) < .Machine$double.eps)
        return(list(x = x2, y = y2))

      xv <- -b_coef / (2 * a_coef)
      xv <- max(min(xv, max(xi)), min(xi))
      yv <- a_coef * xv^2 + b_coef * xv + (y1 - a_coef * x1^2 - b_coef * x1)

      list(x = as.numeric(xv), y = as.numeric(yv))
    }

    max_df <- if (length(idx_max) > 0) {
      verts <- lapply(idx_max, interpolate_vertex)
      data.frame(x = vapply(verts, `[[`, numeric(1), "x"),
                 y = vapply(verts, `[[`, numeric(1), "y"))
    } else data.frame(x = numeric(0), y = numeric(0))

    min_df <- if (length(idx_min) > 0) {
      verts <- lapply(idx_min, interpolate_vertex)
      data.frame(x = vapply(verts, `[[`, numeric(1), "x"),
                 y = vapply(verts, `[[`, numeric(1), "y"))
    } else data.frame(x = numeric(0), y = numeric(0))

    lloq_x <- if (nrow(max_df) > 0) max_df$x[which.max(max_df$y)] else NA_real_
    uloq_x <- if (nrow(min_df) > 0) min_df$x[which.min(min_df$y)] else NA_real_

    lloq_y <- tryCatch(
      if (!is.na(lloq_x)) as.numeric(eval_y(lloq_x)) else NA_real_,
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_loqs] curve_id=%s (%s) lloq_y eval failed: %s",
          cid, method, e$message))
        NA_real_
      }
    )

    uloq_y <- tryCatch(
      if (!is.na(uloq_x)) as.numeric(eval_y(uloq_x)) else NA_real_,
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_loqs] curve_id=%s (%s) uloq_y eval failed: %s",
          cid, method, e$message))
        NA_real_
      }
    )

    if (verbose) message(sprintf(
      "[compute_loqs] curve_id=%s (%s)  lloq_x=%.4f  uloq_x=%.4f  lloq_y=%.4f  uloq_y=%.4f",
      cid, method, lloq_x, uloq_x, lloq_y, uloq_y))

    tibble::tibble(
      curve_id = cid,
      method   = method,
      lloq     = as.numeric(lloq_x),
      uloq     = as.numeric(uloq_x),
      lloq_y   = as.numeric(lloq_y),
      uloq_y   = as.numeric(uloq_y)
    )
  })

  dplyr::bind_rows(results)
}
