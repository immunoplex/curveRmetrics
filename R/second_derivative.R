# =============================================================================
# Adaptive, family-aware numerical second derivative of dose-response curves.
#
# Design notes:
#   - All branching is driven by p$family (from .normalise_params()), which
#     already encodes the frequentist / Bayesian scale distinction.
#   - The working range is derived analytically from the inflection location
#     (p$c) and width parameter (p$b); no stored xmin / xmax is required.
#   - Family-aware working axis:
#       logistic4 / logistic5        — natural-x axis (b has x units)
#       loglogistic4 / loglogistic5  — log10(x) axis  (b is log-slope)
#       gompertz4                    — natural-x axis  (asymmetric span)
#   - Non-uniform central differences give correct d²y/dx² regardless of
#     grid spacing, which is critical for log-spaced grids.
#   - Three-pass adaptive strategy:
#       Pass 1 — coarse scan across the full family-derived range
#       Pass 2 — analytical seed brackets (exact for 4-parameter families)
#       Pass 3 — dense refinement inside every detected bracket
#   Typical output: ~300-400 points concentrated around the two LOQ
#   shoulders rather than uniformly spread across the full range.
# =============================================================================


# -----------------------------------------------------------------------------
# .compute_second_deriv
# -----------------------------------------------------------------------------

#' Compute an adaptive numerical second derivative for a single curve
#'
#' Evaluates \eqn{d^2y/dx^2} over a family-appropriate concentration range
#' using a three-pass adaptive grid: a coarse scan to locate candidate
#' extrema, analytical seed brackets to ensure the LOQ shoulders are captured,
#' and a dense refinement pass inside each detected bracket.
#'
#' All output concentrations and second-derivative values are on the
#' natural-concentration axis regardless of whether the curve was evaluated
#' internally on a log10 grid; .d2y_nonuniform() handles the Jacobian
#' correction automatically.
#'
#' @param p        Named list as returned by .normalise_params().  Must
#'   contain `family`, `a`, `b`, `c`, `d`, and `g`.
#' @param n_coarse Integer. Number of points in the initial coarse scan.
#'   Default `80L`.
#' @param n_refine Integer. Number of points per refinement bracket.
#'   Default `80L`.
#'
#' @return Tibble with columns:
#'   \describe{
#'     \item{`concentration`}{Natural-scale concentration values.}
#'     \item{`d2x_y`}{\eqn{d^2y/dx^2} evaluated at each concentration.
#'       `NA` at the two endpoints of each grid segment.}
#'   }
#'
#' @keywords internal
.compute_second_deriv <- function(p, n_coarse = 80L, n_refine = 80L) {

  # ── 1. Family-aware working axis and bounds ──────────────────────────────────
  si <- .sderiv_scale(p)

  # ── 2. Analytical seed positions (u-domain) ──────────────────────────────────
  seeds_u <- .sderiv_seeds(p, si)

  # ── 3. Coarse grid ───────────────────────────────────────────────────────────
  u_coarse <- seq(si$u_lo, si$u_hi, length.out = n_coarse)
  x_coarse <- if (si$log_scale) 10^u_coarse else u_coarse
  y_coarse <- .sderiv_eval_vec(p, x_coarse)

  # ── 4. Detect extremum brackets from coarse scan ────────────────────────────
  d2y_coarse <- .d2y_nonuniform(x_coarse, y_coarse)
  bracks     <- .sderiv_brackets(u_coarse, d2y_coarse)

  # ── 5. Seed-centred brackets ─────────────────────────────────────────────────
  du_c <- (si$u_hi - si$u_lo) / n_coarse
  seed_bracks <- lapply(seeds_u, function(s)
    c(max(si$u_lo, s - 4 * du_c), min(si$u_hi, s + 4 * du_c))
  )
  all_bracks <- c(bracks, seed_bracks)

  # ── 6. Refined points inside each bracket ────────────────────────────────────
  refined_u <- unlist(lapply(all_bracks, function(br)
    seq(br[1], br[2], length.out = n_refine)
  ))

  # ── 7. Merge coarse + refined, evaluate once, compute d2y once ───────────────
  # Coarse points are retained throughout so density ramps up smoothly toward
  # each shoulder rather than stepping sharply at bracket boundaries.
  # A single-pass d2y on the full merged grid eliminates seam NAs entirely.
  u_all <- sort(unique(c(u_coarse, refined_u)))
  x_all <- if (si$log_scale) 10^u_all else u_all

  y_all   <- .sderiv_eval_vec(p, x_all)
  d2y_all <- .d2y_nonuniform(x_all, y_all)

  tibble::tibble(concentration = x_all, d2x_y = d2y_all)
}


# -----------------------------------------------------------------------------
# .sderiv_scale
# -----------------------------------------------------------------------------

#' Determine the family-appropriate working axis and concentration bounds
#'
#' Returns the bounds and axis type used by .compute_second_deriv() to
#' construct its concentration grid.  The working variable \eqn{u} is either
#' the natural concentration or its log10 transform, depending on the model
#' family:
#' \itemize{
#'   \item `logistic4` / `logistic5` — \eqn{u = x}; span is \eqn{\pm 6|b|}
#'     around the inflection point \eqn{c}.
#'   \item `loglogistic4` / `loglogistic5` — \eqn{u = \log_{10}(x)}; span is
#'     \eqn{\pm 6/|b|} around \eqn{\log_{10}(c)}.
#'   \item `gompertz4` — \eqn{u = x}; the left span is doubled relative to
#'     the right to account for the curve's left-skewed shape.
#' }
#'
#' @param p Named list as returned by .normalise_params().
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{`u_lo`}{Lower bound of the working variable.}
#'     \item{`u_hi`}{Upper bound of the working variable.}
#'     \item{`log_scale`}{Logical; `TRUE` when \eqn{u = \log_{10}(x)}.}
#'   }
#'
#' @keywords internal
#' @noRd
.sderiv_scale <- function(p) {
  bw <- abs(p$b)
  cc <- p$c

  switch(p$family,

         "logistic4" = ,
         "logistic5" = {
           span <- max(6 * bw, 1e-9)
           list(u_lo = cc - span, u_hi = cc + span, log_scale = FALSE)
         },

         "loglogistic4" = ,
         "loglogistic5" = {
           u_c  <- log10(max(cc, 1e-300))
           span <- max(6 / bw, 1e-9)
           list(u_lo = u_c - span, u_hi = u_c + span, log_scale = TRUE)
         },

         "gompertz4" = {
           span <- max(6 / bw, 1e-9)
           list(u_lo = cc - 2 * span, u_hi = cc + span, log_scale = FALSE)
         },

         stop(sprintf("[.sderiv_scale] Unknown family '%s'", p$family), call. = FALSE)
  )
}


# -----------------------------------------------------------------------------
# .sderiv_seeds
# -----------------------------------------------------------------------------

#' Compute analytical seed positions for second-derivative extrema
#'
#' Returns the concentration(s) at which \eqn{d^3y/dx^3 = 0}, i.e. the exact
#' positions of the local extrema of \eqn{d^2y/dx^2}.  These seed positions
#' are used by .compute_second_deriv() to ensure refinement brackets are
#' placed accurately even when the coarse scan misses a narrow shoulder.
#'
#' @details
#' Derivations by family:
#' \describe{
#'   \item{`logistic4`}{Set \eqn{d^3y/dx^3 = 0}: \eqn{z^2 - 4z + 1 = 0},
#'     \eqn{z = \exp((x-c)/b)}.  Exact roots \eqn{z^* = 2 \pm \sqrt{3}} give
#'     \eqn{x^* = c + b \ln(2 \pm \sqrt{3})}.}
#'   \item{`logistic5`}{No closed form for arbitrary \eqn{g}.  The 4-parameter
#'     offset is scaled by \eqn{g^{0.3}} as a heuristic (\eqn{g > 1} widens,
#'     \eqn{g < 1} narrows).}
#'   \item{`loglogistic4`}{Exact: \eqn{x^* = c \cdot (2 \pm \sqrt{3})^{1/b}}.}
#'   \item{`loglogistic5`}{Heuristic: exponent adjusted by \eqn{g^{0.3}/b}.}
#'   \item{`gompertz4`}{Set \eqn{d^3y/dx^3 = 0}: \eqn{v^2 - 3v + 1 = 0},
#'     \eqn{v = \exp(-b(x-c))}.  Exact roots \eqn{v^* = \varphi^{\pm 2}}
#'     (\eqn{\varphi} = golden ratio) give \eqn{x^* = c \pm 2\ln\varphi / b}.
#'     The two shoulders are equidistant from \eqn{c} in x; asymmetry appears
#'     only in the heights of the \eqn{d^2y} peaks.}
#' }
#'
#' @param p  Named list as returned by .normalise_params().
#' @param si Named list as returned by .sderiv_scale().
#'
#' @return Numeric vector of seed positions clipped to
#'   `[si$u_lo, si$u_hi]`, in the working \eqn{u}-domain.
#'
#' @keywords internal
#' @noRd
.sderiv_seeds <- function(p, si) {
  bw      <- abs(p$b)
  cc      <- p$c
  gg      <- p$g
  S3      <- sqrt(3)
  PHI     <- (1 + sqrt(5)) / 2
  LN_PHI2 <- 2 * log(PHI)

  nat_seeds <- switch(p$family,

                      "logistic4" = {
                        cc + p$b * log(c(2 - S3, 2 + S3))
                      },

                      "logistic5" = {
                        scale <- bw * max(gg, 1e-9)^0.3
                        cc + sign(p$b) * scale * log(c(2 - S3, 2 + S3))
                      },

                      "loglogistic4" = {
                        cc * c((2 - S3)^(1 / bw), (2 + S3)^(1 / bw))
                      },

                      "loglogistic5" = {
                        fac <- max(gg, 1e-9)^0.3 / bw
                        cc * c((2 - S3)^fac, (2 + S3)^fac)
                      },

                      "gompertz4" = {
                        c(cc - LN_PHI2 / p$b,
                          cc + LN_PHI2 / p$b)
                      }
  )

  u_seeds <- if (si$log_scale) log10(pmax(nat_seeds, 1e-300)) else nat_seeds
  pmax(si$u_lo, pmin(si$u_hi, u_seeds))
}


# -----------------------------------------------------------------------------
# .sderiv_brackets
# -----------------------------------------------------------------------------

#' Locate working-axis intervals that bracket local extrema of d²y
#'
#' Scans the coarse second-derivative vector for sign changes in consecutive
#' first differences, which indicate a local extremum of \eqn{d^2y/dx^2}.
#' Each detected extremum is returned as a bracket \eqn{[u_\text{lo},
#' u_\text{hi}]} widened by one index on each side for numerical safety.
#'
#' @param u   Numeric vector of working-axis positions (sorted, no `NA`).
#' @param d2y Numeric vector of \eqn{d^2y/dx^2} values corresponding to `u`
#'   (may contain `NA`).
#'
#' @return List of two-element numeric vectors `c(u_lo, u_hi)`, one per
#'   detected extremum.  Returns an empty list if fewer than four valid
#'   points are available.
#'
#' @keywords internal
#' @noRd
.sderiv_brackets <- function(u, d2y) {
  ok    <- !is.na(d2y)
  u_v   <- u[ok]; d2y_v <- d2y[ok]
  if (length(d2y_v) < 4L) return(list())

  dd       <- diff(d2y_v)
  sign_chg <- which(dd[-length(dd)] * dd[-1] < 0)

  lapply(sign_chg, function(i) {
    lo <- max(1L, i - 1L)
    hi <- min(length(u_v), i + 2L)
    c(u_v[lo], u_v[hi])
  })
}


# -----------------------------------------------------------------------------
# .d2y_nonuniform
# -----------------------------------------------------------------------------

#' Non-uniform central-difference second derivative
#'
#' Computes \eqn{d^2y/dx^2} at interior points of an arbitrarily-spaced
#' grid using the three-point non-uniform formula:
#'
#' \deqn{
#'   \left.\frac{d^2y}{dx^2}\right|_i \approx
#'   \frac{2}{h_l + h_r}
#'   \left[\frac{y_{i+1} - y_i}{h_r} - \frac{y_i - y_{i-1}}{h_l}\right]
#' }
#'
#' where \eqn{h_l = x_i - x_{i-1}} and \eqn{h_r = x_{i+1} - x_i}.
#'
#' This gives a correct derivative on the natural-concentration axis
#' regardless of whether \eqn{x} was sampled uniformly or on a log scale,
#' making it suitable for both logistic and log-logistic families.
#'
#' @param x Natural-scale concentration vector (sorted ascending, no `NA`).
#' @param y Curve response values at each element of `x`.
#'
#' @return Numeric vector of the same length as `x`.  The two endpoint
#'   elements are always `NA`; interior elements are `NA` where any of the
#'   three required `y` values are `NA` or where consecutive `x` values are
#'   closer than `.Machine$double.eps`.
#'
#' @keywords internal
#' @noRd
.d2y_nonuniform <- function(x, y) {
  n  <- length(x)
  d2 <- rep(NA_real_, n)
  for (i in seq(2L, n - 1L)) {
    if (anyNA(y[c(i - 1L, i, i + 1L)])) next
    h_l <- x[i]      - x[i - 1L]
    h_r <- x[i + 1L] - x[i]
    if (h_l < .Machine$double.eps || h_r < .Machine$double.eps) next
    d2[i] <- 2 * ((y[i + 1L] - y[i]) / h_r - (y[i] - y[i - 1L]) / h_l) /
      (h_l + h_r)
  }
  d2
}


# -----------------------------------------------------------------------------
# .sderiv_eval_vec
# -----------------------------------------------------------------------------

#' Vectorised curve evaluation with per-point error guard
#'
#' Attempts a single vectorised call to [.evaluate_curve()] (fast path) and
#' falls back to a point-by-point loop if the batch call fails — for example,
#' when log-logistic formulas encounter \eqn{x \le 0} at the grid boundary.
#'
#' @param p     Named list as returned by [.normalise_params()].
#' @param x_vec Numeric vector of concentration values to evaluate.
#'
#' @return Numeric vector of predicted response values, same length as
#'   `x_vec`.  Failed evaluations return `NA_real_`.
#'
#' @keywords internal
#' @noRd
.sderiv_eval_vec <- function(p, x_vec) {
  tryCatch(
    as.numeric(.evaluate_curve(p, x_vec)),
    error = function(e) {
      vapply(x_vec,
             function(xv) tryCatch(as.numeric(.evaluate_curve(p, xv)),
                                   error = function(e2) NA_real_),
             numeric(1L))
    }
  )
}


# -----------------------------------------------------------------------------
# compute_second_deriv_df
# -----------------------------------------------------------------------------

#' Compute second derivatives for all curves in a data frame
#'
#' Applies [.compute_second_deriv()] to every row of `curves_df` and
#' returns the results as a single long-format tibble keyed by `curve_id`.
#' Curves for which parameter normalisation or derivative computation fails
#' are included with a single `NA` row so that downstream joins do not
#' silently drop curve identifiers.
#'
#' The returned `concentration` column is always on the natural scale.
#' For Bayesian curves the concentrations therefore reflect actual
#' concentration units; for frequentist curves they are on the log10 scale
#' (matching `p$c` after [.normalise_params()]).
#'
#' @param curves_df Data frame of fitted curves, typically `freq_curves_df`
#'   or `bayes_curves_df`, containing at minimum the columns `curve_id`,
#'   `method`, `model_name`, `a`, `b`, `c`, `d`, and `g`.
#' @param n_coarse  Integer. Coarse-scan grid size passed to
#'   .compute_second_deriv(). Default `80L`.
#' @param n_refine  Integer. Refinement grid size per bracket passed to
#'   .compute_second_deriv(). Default `80L`.
#' @param verbose   Logical. If `TRUE`, emits a diagnostic message for each
#'   curve that fails normalisation or derivative computation.
#'
#' @return Tibble with columns:
#'   \describe{
#'     \item{`curve_id`}{Curve identifier, matching `curves_df$curve_id`.}
#'     \item{`concentration`}{Natural-scale concentration. `NA` for failed
#'       curves.}
#'     \item{`d2x_y`}{\eqn{d^2y/dx^2} at each concentration. `NA` at grid
#'       endpoints and for failed curves.}
#'   }
#'   Endpoint `NA` rows are removed before returning; only interior `NA`
#'   values (from domain errors) are retained.
#'
#' @export
compute_second_deriv_df <- function(curves_df,
                                    n_coarse = 80L,
                                    n_refine = 80L,
                                    verbose  = TRUE) {

  purrr::map_dfr(seq_len(nrow(curves_df)), function(i) {

    row <- curves_df[i, ]
    cid <- row$curve_id

    p <- tryCatch(
      .normalise_params(row, row$method),
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_second_deriv_df] curve_id=%s normalise failed: %s",
          cid, e$message))
        NULL
      }
    )

    if (is.null(p))
      return(tibble::tibble(curve_id      = cid,
                            concentration = NA_real_,
                            d2x_y         = NA_real_))

    result <- tryCatch(
      .compute_second_deriv(p, n_coarse = n_coarse, n_refine = n_refine),
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_second_deriv_df] curve_id=%s deriv failed: %s",
          cid, e$message))
        tibble::tibble(concentration = NA_real_, d2x_y = NA_real_)
      }
    )

    result            <- result[!is.na(result$d2x_y), ]
    result$curve_id   <- cid
    result[, c("curve_id", "concentration", "d2x_y")]
  })
}
