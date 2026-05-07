#' Compute analytical inflection x-coordinate from normalised parameters
#'
#' All inputs are assumed to be on the natural concentration scale.
#' Derivations:
#'
#' logistic4:
#'   y = d + (a-d) / (1 + exp((x-c)/b))
#'   Symmetric => inflection exactly at x = c, y = (a+d)/2
#'
#' logistic5:
#'   y = d + (a-d) / (1 + exp((x-c)/b))^g
#'   Set d²y/dx² = 0: let z = exp((x-c)/b)
#'   Numerator of d²y/dx² proportional to: (g+1)*z - (b-1)/(b*... )
#'   Full derivation => z* = (b-1) / (b*g + 1)  requires b > 1 and z* > 0
#'   => x_infl = c + b * log(z*)
#'
#' loglogistic4:
#'   y = a + (d-a) / (1 + (x/c)^b)
#'   Symmetric on log(x) scale => inflection at x = c
#'
#' loglogistic5:
#'   y = a + (d-a) * (1 + g*exp(-b*(x-c)))^(-1/g)
#'   Set d²y/dx² = 0: let u = g*exp(-b*(x-c))
#'   => u* = (b-1)/(b+1)  ... simplified Richards condition  requires b > 1
#'   => x_infl = c - (log(u*) - log(g)) / b
#'             = c + (log(g) - log((b-1)/(b+1))) / b
#'
#' gompertz4:
#'   y = a + (d-a) * exp(-exp(-b*(x-c)))
#'   Set d²y/dx² = 0: let u = exp(-b*(x-c))
#'   => u*(1 - u) = 1/b ... simplifies to u* = (b-1)/b  (requires b > 1)
#'   => x_infl = c - log(u*) / b = c + log(b/(b-1)) / b
#'   At inflection: y_infl = a + (d-a)*exp(-1)  (always ~36.8% of range from a)
#'
#' @param p Named list with family, a, b, c, d, g (all natural scale)
#' @return Numeric scalar, inflection x on natural concentration scale
#' @keywords internal
.inflection_x <- function(p) {

  b <- p$b
  g <- p$g
  c <- p$c   # natural scale throughout

  switch(p$family,

         "logistic4" = {
           # Symmetric: inflection exactly at x = c
           c
         },

         "logistic5" = {
           # z* = (b-1)/(b*g+1), requires b > 1 and z* > 0
           if (b <= 1) warning(sprintf(
             "[.inflection_x] logistic5: b = %.4f <= 1; inflection point may be unreliable.", b
           ))
           z_star <- (b - 1) / (b * g + 1)
           if (!is.finite(z_star) || z_star <= 0) {
             warning("[.inflection_x] logistic5: z_star <= 0; falling back to c.")
             return(c)
           }
           c + b * log(z_star)
         },

         "loglogistic4" = {
           # Inflection at x = c (symmetric on log-x scale)
           c
         },

         "loglogistic5" = {
           # u* = (b-1)/(b+1), requires b > 1
           if (b <= 1) warning(sprintf(
             "[.inflection_x] loglogistic5: b = %.4f <= 1; inflection point may be unreliable.", b
           ))
           u_star <- (b - 1) / (b + 1)
           if (!is.finite(u_star) || u_star <= 0) {
             warning("[.inflection_x] loglogistic5: u_star <= 0; falling back to c.")
             return(c)
           }
           # x_infl = c + (log(g) - log(u_star)) / b
           c + (log(g) - log(u_star)) / b
         },

         "gompertz4" = {
           # u* = (b-1)/b, requires b > 1
           if (b <= 1) warning(sprintf(
             "[.inflection_x] gompertz4: b = %.4f <= 1; falling back to c.", b
           ))
           if (b <= 1) return(c)
           u_star <- (b - 1) / b
           # x_infl = c + log(b/(b-1)) / b = c - log(u_star) / b
           c - log(u_star) / b
         },

         stop(sprintf("[.inflection_x] Unsupported family: '%s'", p$family), call. = FALSE)
  )
}

#' Compute inflection point for a single curve row
#'
#' Worker called by [compute_inflection_point()]. Normalises parameters,
#' computes the analytical inflection x-coordinate via [.inflection_x()],
#' and evaluates the forward curve at that point.
#'
#' @param row  Single-row data frame from `freq_curves_df` or
#'   `bayes_curves_df`, containing at minimum `method`, `model_name`,
#'   `a`, `b`, `c`, `d`, `g`.
#' @param verbose Logical. If `TRUE`, emits diagnostic messages.
#'
#' @return Named list with elements `inflect_x`, `inflect_y`, `family`,
#'   `source`, and `params` (the full normalised parameter list).
#'
#' @keywords internal
#' @noRd
compute_inflection_point_worker <- function(row, verbose = TRUE) {

  method <- row$method   # "frequentist" or "bayesian"
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

  # if (method == "bayesian") {
  #   inflect_x <- log10(inflect_x)
  #   inflect_y <- log10(inflect_y)
  # }

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

#' Add inflection-point coordinates to a curves data frame
#'
#' Applies compute_inflection_point_worker() to every row of `curves_df`
#' and appends two new columns: `inflect_x` (concentration at the inflection
#' point) and `inflect_y` (predicted response at the inflection point).
#' Both values are on the same scale as `c` after normalisation by
#' .normalise_params().
#'
#' @param curves_df Data frame of fitted curves, typically `freq_curves_df` or
#'   `bayes_curves_df`, containing at minimum `method`, `model_name`,
#'   `a`, `b`, `c`, `d`, `g`.
#'
#' @return `curves_df` with two additional numeric columns: `inflect_x` and
#'   `inflect_y`.
#'
#' @export
compute_inflection_point <- function(curves_df) {
  results <- purrr::map(
    seq_len(nrow(curves_df)),
    ~ compute_inflection_point_worker(curves_df[.x, ], verbose = FALSE)
  )
  curves_df$inflect_x <- vapply(results, `[[`, numeric(1), "inflect_x")
  curves_df$inflect_y <- vapply(results, `[[`, numeric(1), "inflect_y")
  curves_df
}

# -----------------------------------------------------------------------------
# compute_assay_sensitivity
# -----------------------------------------------------------------------------

#' Estimate assay sensitivity as the slope at the inflection point
#'
#' Computes the first derivative of the forward curve at the inflection point
#' using a central-difference approximation.  A larger absolute value indicates
#' greater assay sensitivity (larger response change per unit change in
#' log10-concentration).
#'
#' The inflection x-coordinate is taken directly from `curves_df$inflect_x`;
#' run [compute_inflection_point()] on the data frame first.
#'
#' @param curves_df Data frame of fitted curves with an `inflect_x` column,
#'   as produced by [compute_inflection_point()].
#' @param verbose Logical. If `TRUE`, emits one diagnostic message per curve.
#'
#' @return Tibble with columns:
#'   \describe{
#'     \item{`curve_id`}{Curve identifier.}
#'     \item{`method`}{`"frequentist"` or `"bayesian"`.}
#'     \item{`inflect_x`}{Inflection-point concentration (log10 scale).}
#'     \item{`dydx_inflect`}{First derivative of response with respect to
#'       log10-concentration at the inflection point.}
#'   }
#'
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

    # if (row$method == "bayesian") {
    #  dydx <- dydx * ix * log(10) # apply chain-rule to get to correct log10 scale.
    # }

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
