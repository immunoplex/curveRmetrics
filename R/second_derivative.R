# =============================================================================
# second_deriv_natural.R
#
# Adaptive, family-aware numerical second derivative of a dose-response curve.
#
# Key design decisions:
#
#   * ALL families use log_scale = TRUE (working axis u = log10(x)).
#
#   * d²y/du² (u = log10(x)) is computed — stable on log-spaced grids.
#     Using natural x would cause 1/x² blow-up at small concentrations.
#
#   * Curve evaluation uses the actual plotted curve for each method:
#       Frequentist : d²(log10 y)/d(log10 x)²  via  10^.evaluate_curve(p_raw, log10(x))
#       Bayesian    : d²(log10 y)/d(log10 x)²  via  .evaluate_curve(p_nat, log(x))
#     p$eval_fn carries the correct evaluator; p$b and p$c carry the
#     natural-scale params used only for grid sizing (.sderiv_scale /
#     .sderiv_seeds).
#
#   * local() is used in compute_second_deriv_df to force eager capture of
#     p_raw / p_nat inside the purrr loop — avoids R lazy-eval scoping bug.
#
#   * Span formula per family:
#       logistic4/5, gompertz4 : ±6·b_nat/ln(10) log10 units
#       loglogistic4/5         : ±6/b log10 units
#
#   * Seed formula table (all return log10(x*)):
#       logistic4    : u_c ± b_nat · log10(2 ± √3)               [exact]
#       logistic5    : logistic4 seeds scaled by g^0.3            [approx]
#       loglogistic4 : u_c ± log10(2 ± √3) / b                   [exact]
#       loglogistic5 : loglogistic4 seeds scaled by g^0.3         [approx]
#       gompertz4    : u_c ± LN_PHI2 / (b_nat · ln10)            [exact]
# =============================================================================


# Shared numeric constants
LN10    <- log(10)                    # ≈ 2.302585
LN_PHI2 <- log((3 + sqrt(5)) / 2)    # 2·ln(golden ratio) ≈ 0.9624
S3      <- sqrt(3)                    # ≈ 1.7321


#' Compute the working-axis bounds for the second derivative grid
#'
#' Returns the \eqn{\log_{10}(x)} interval that spans the region of interest
#' for the second derivative, centred on \eqn{u_c = \log_{10}(c)} and widened
#' by a family-specific multiple of the slope parameter:
#' \describe{
#'   \item{logistic4/5, gompertz4}{span = \eqn{6 \cdot b_{nat} / \ln(10)},
#'     scaled by \eqn{g^{0.3}} for the 5-parameter variant.}
#'   \item{loglogistic4/5}{span = \eqn{6 / b}.}
#'   \item{gompertz4}{left tail extended by \eqn{2 \times} span to capture
#'     the asymmetric lower extremum.}
#' }
#'
#' @param p A named list of natural-scale curve parameters with elements
#'   \code{family}, \code{b} (natural-log-scale slope for logistic/gompertz;
#'   log10-scale slope for loglogistic), \code{c} (natural EC50), and
#'   optionally \code{g}.
#'
#' @return A named list with elements \code{u_lo}, \code{u_hi} (the
#'   \eqn{\log_{10}(x)} bounds), and \code{log_scale = TRUE}.
#' @keywords internal
.sderiv_scale <- function(p) {

  bw  <- abs(p$b)
  u_c <- log10(max(p$c, 1e-300))

  switch(p$family,

         "logistic4" = ,
         "logistic5" = {
           gg   <- if (p$family == "logistic5") max(p$g, 1e-9) else 1
           span <- max(6 * bw / LN10 * max(gg, 1)^0.3, 1e-9)
           list(u_lo = u_c - span, u_hi = u_c + span, log_scale = TRUE)
         },

         "loglogistic4" = ,
         "loglogistic5" = {
           span <- max(6 / bw, 1e-9)
           list(u_lo = u_c - span, u_hi = u_c + span, log_scale = TRUE)
         },

         "gompertz4" = {
           span <- max(6 * bw / LN10, 1e-9)
           list(u_lo = u_c - 2 * span, u_hi = u_c + span, log_scale = TRUE)
         },

         stop(sprintf("[.sderiv_scale] Unknown family '%s'", p$family), call. = FALSE)
  )
}


#' Analytical seed positions for the extrema of \eqn{d^2y/du^2}
#'
#' Returns a two-element numeric vector of \eqn{\log_{10}(x^*)} values at
#' which \eqn{d^3y/du^3 = 0}, i.e. the expected locations of the maximum and
#' minimum of the second derivative.  These seeds are used as starting points
#' for the dense refinement grid in \code{\link{.compute_second_deriv}}.
#'
#' Derivations use \eqn{t = \ln(x/c)}, \eqn{u = \log_{10}(x)}:
#' \describe{
#'   \item{logistic4}{\eqn{z^2 - 4z + 1 = 0},
#'     \eqn{z^* = 2 \pm \sqrt{3}},
#'     \eqn{u^* = u_c + b \cdot \log_{10}(2 \pm \sqrt{3})} (exact).}
#'   \item{logistic5}{logistic4 seeds scaled by \eqn{g^{0.3}} (approx).}
#'   \item{loglogistic4}{\eqn{u^* = u_c \pm \log_{10}(2 \pm \sqrt{3}) / b} (exact).}
#'   \item{loglogistic5}{loglogistic4 seeds scaled by \eqn{g^{0.3}} (approx).}
#'   \item{gompertz4}{\eqn{v^* = (3 \pm \sqrt{5})/2},
#'     \eqn{u^* = u_c \pm \mathrm{LN\_PHI2} / (b \cdot \ln 10)} (exact).}
#' }
#'
#' @param p A named list of natural-scale curve parameters (see
#'   \code{\link{.sderiv_scale}}).
#'
#' @return A numeric vector of length 2 containing the two seed positions in
#'   \eqn{\log_{10}(x)} units.
#' @keywords internal
.sderiv_seeds <- function(p) {

  bw  <- abs(p$b)
  u_c <- log10(max(p$c, 1e-300))

  switch(p$family,

         "logistic4" = {
           u_c + bw * log10(c(2 - S3, 2 + S3))
         },

         "logistic5" = {
           gg <- max(p$g, 1e-9)
           u_c + bw * gg^0.3 * log10(c(2 - S3, 2 + S3))
         },

         "loglogistic4" = {
           u_c + log10(c(2 - S3, 2 + S3)) / bw
         },

         "loglogistic5" = {
           gg <- max(p$g, 1e-9)
           u_c + (log10(c(2 - S3, 2 + S3)) / bw) * gg^0.3
         },

         "gompertz4" = {
           u_c + c(-LN_PHI2, LN_PHI2) / (bw * LN10)
         },

         stop(sprintf("[.sderiv_seeds] Unknown family '%s'", p$family), call. = FALSE)
  )
}


#' Locate u-intervals that bracket local extrema of \eqn{d^2y/du^2}
#'
#' Scans a discrete second-derivative series for sign changes in consecutive
#' first differences, which indicate the presence of a local maximum or
#' minimum.  For each detected sign change, a bracketing interval of
#' neighbouring \eqn{u = \log_{10}(x)} values is returned.
#'
#' @param u Numeric vector of \eqn{\log_{10}(x)} grid positions.
#' @param d2y Numeric vector of second-derivative values at each position in
#'   \code{u}; \code{NA} values are removed before scanning.
#'
#' @return A list of two-element numeric vectors, each giving the
#'   \eqn{[u_{lo},\, u_{hi}]} bounds of one bracket.  An empty list is
#'   returned when fewer than four non-\code{NA} points are available or no
#'   sign change is detected.
#' @keywords internal
.sderiv_brackets <- function(u, d2y) {
  ok    <- !is.na(d2y)
  u_v   <- u[ok];  d2y_v <- d2y[ok]
  if (length(d2y_v) < 4L) return(list())

  dd       <- diff(d2y_v)
  sign_chg <- which(dd[-length(dd)] * dd[-1] < 0)

  lapply(sign_chg, function(i) {
    lo <- max(1L, i - 1L)
    hi <- min(length(u_v), i + 2L)
    c(u_v[lo], u_v[hi])
  })
}


#' Non-uniform central-difference second derivative in log10-x space
#'
#' Computes \eqn{d^2y/du^2} where \eqn{u = \log_{10}(x)} using a
#' non-uniform central-difference formula that is stable on log-spaced grids.
#' Working in \eqn{\log_{10}(x)} space avoids the \eqn{1/x^2} blow-up that
#' would occur when differentiating in natural concentration.
#'
#' The formula at interior point \eqn{i} is:
#' \deqn{
#'   \frac{d^2y}{du^2}\bigg|_i \approx
#'   \frac{2}{h_l + h_r}
#'   \left[\frac{y_{i+1} - y_i}{h_r} - \frac{y_i - y_{i-1}}{h_l}\right]
#' }
#' where \eqn{h_l = u_i - u_{i-1}} and \eqn{h_r = u_{i+1} - u_i}.
#'
#' @param x Numeric vector of natural-scale concentration values (strictly
#'   positive, log-spaced).
#' @param y Numeric vector of response values corresponding to \code{x}.
#'
#' @return Numeric vector of the same length as \code{x} containing
#'   \eqn{d^2y/du^2}; the first and last elements are always \code{NA}, as
#'   are any interior points where \code{y} contains \code{NA} or where
#'   adjacent grid spacings are below machine epsilon.
#' @keywords internal
.d2y_nonuniform <- function(x, y) {
  u  <- log10(x)
  n  <- length(u)
  d2 <- rep(NA_real_, n)
  for (i in seq(2L, n - 1L)) {
    if (anyNA(y[c(i - 1L, i, i + 1L)])) next
    h_l <- u[i]      - u[i - 1L]
    h_r <- u[i + 1L] - u[i]
    if (h_l < .Machine$double.eps || h_r < .Machine$double.eps) next
    d2[i] <- 2 * ((y[i + 1L] - y[i]) / h_r -
                    (y[i]       - y[i - 1L]) / h_l) / (h_l + h_r)
  }
  d2
}


#' Vectorised curve evaluation with per-point error guard
#'
#' Attempts to evaluate \code{p$eval_fn} (or \code{\link{.evaluate_curve}} as
#' a fallback) over the entire vector \code{x_vec} in a single call.  If that
#' call throws an error, each point is evaluated individually and failures are
#' replaced with \code{NA}.
#'
#' @param p A named list of curve parameters.  If \code{p$eval_fn} is a
#'   function it is used as the evaluator; otherwise
#'   \code{\link{.evaluate_curve}} is called with \code{p} and \code{x}.
#' @param x_vec Numeric vector of concentration values to evaluate.
#'
#' @return Numeric vector of the same length as \code{x_vec} with
#'   \code{NA} substituted for any evaluation failure.
#' @keywords internal
.sderiv_eval_vec <- function(p, x_vec) {
  eval_fn <- if (!is.null(p$eval_fn)) p$eval_fn else function(x) .evaluate_curve(p, x)
  tryCatch(
    as.numeric(eval_fn(x_vec)),
    error = function(e) {
      vapply(x_vec,
             function(xv) tryCatch(as.numeric(eval_fn(xv)),
                                   error = function(e2) NA_real_),
             numeric(1L))
    }
  )
}


#' Three-pass adaptive numerical second derivative for a single curve
#'
#' Computes \eqn{d^2y/du^2} (\eqn{u = \log_{10}(x)}) over a family-derived
#' concentration grid using three passes:
#'
#' \enumerate{
#'   \item \strong{Coarse pass} — \code{n_coarse} evenly spaced
#'     \eqn{\log_{10}(x)} points across the family bounds from
#'     \code{\link{.sderiv_scale}}.
#'   \item \strong{Seed brackets} — \code{\link{.sderiv_seeds}} provides
#'     analytical estimates of the extremum locations; each seed is bracketed
#'     by ±4 coarse grid spacings.
#'   \item \strong{Refinement} — \code{n_refine} dense points are placed
#'     inside every bracket detected from both the coarse pass sign-changes
#'     (\code{\link{.sderiv_brackets}}) and the seed brackets.
#' }
#'
#' All three pass grids are merged, deduplicated, and evaluated in a single
#' call to \code{\link{.sderiv_eval_vec}}, after which
#' \code{\link{.d2y_nonuniform}} computes the final derivative series.
#'
#' @param p A named list of natural-scale curve parameters including a
#'   \code{family} string, \code{b}, \code{c}, \code{g}, and an
#'   \code{eval_fn} closure that maps concentration to response (see
#'   \code{\link{compute_second_deriv_df}}).
#' @param n_coarse Integer; number of points in the initial coarse grid
#'   (default \code{80L}).
#' @param n_refine Integer; number of points added inside each refinement
#'   bracket (default \code{80L}).
#'
#' @return A \code{\link[tibble:tibble]{tibble}} with columns
#'   \code{concentration} (natural scale) and \code{d2x_y}
#'   (\eqn{d^2y/du^2}).
#' @keywords internal
.compute_second_deriv <- function(p, n_coarse = 80L, n_refine = 80L) {

  si      <- .sderiv_scale(p)
  seeds_u <- .sderiv_seeds(p)

  # ── Pass 1: coarse grid ────────────────────────────────────────────────────
  u_coarse <- seq(si$u_lo, si$u_hi, length.out = n_coarse)
  x_coarse <- 10^u_coarse
  y_coarse <- .sderiv_eval_vec(p, x_coarse)

  # ── Pass 2: detect extremum brackets from coarse d²y ──────────────────────
  d2y_coarse <- .d2y_nonuniform(x_coarse, y_coarse)
  bracks     <- .sderiv_brackets(u_coarse, d2y_coarse)

  # ── Seed-centred brackets (±4 coarse spacings around each analytical seed) ─
  du_c        <- (si$u_hi - si$u_lo) / n_coarse
  seed_bracks <- lapply(seeds_u, function(s)
    c(max(si$u_lo, s - 4 * du_c), min(si$u_hi, s + 4 * du_c))
  )
  all_bracks <- c(bracks, seed_bracks)

  # ── Pass 3: dense refinement inside every bracket ─────────────────────────
  refined_u <- unlist(lapply(all_bracks, function(br)
    seq(br[1], br[2], length.out = n_refine)
  ))

  # ── Merge coarse + refined; single-pass evaluate + d2y ────────────────────
  u_all   <- sort(unique(c(u_coarse, refined_u)))
  x_all   <- 10^u_all
  y_all   <- .sderiv_eval_vec(p, x_all)
  d2y_all <- .d2y_nonuniform(x_all, y_all)

  tibble::tibble(concentration = x_all, d2x_y = d2y_all)
}


#' Compute the second derivative series for all curves in a data frame
#'
#' Public entry point.  Iterates over unique \code{(curve_id, method)}
#' combinations and calls \code{\link{.compute_second_deriv}} for each,
#' returning a long-format data frame of second derivative values suitable
#' for use in \code{\link{compute_loqs}}.
#'
#' A method-specific parameter list is built for each curve:
#'
#' \describe{
#'   \item{Frequentist}{Grid parameters use \eqn{b_{nat} = b_{raw} \cdot \ln(10)}
#'     and \eqn{c_{nat} = 10^{c_{raw}}}.  The \code{eval_fn} closure computes
#'     \eqn{d^2(\log_{10} y)/d(\log_{10} x)^2} via
#'     \code{log10(10^.evaluate_curve(p_raw, log10(x)))}.}
#'   \item{Bayesian}{Grid parameters use \eqn{b_{nat} = 1/b_{raw}} and the
#'     raw natural EC50 \eqn{c_{raw}}.  The \code{eval_fn} closure computes
#'     \eqn{d^2(\log_{10} y)/d(\log_{10} x)^2} via
#'     \code{log10(.evaluate_curve(p_nat, log(x)))}.}
#' }
#'
#' \code{local()} is used when constructing each \code{eval_fn} to force
#' eager capture of \code{p_raw} / \code{p_nat}, preventing R's
#' lazy-evaluation scoping from causing all closures to share the same
#' (final-iteration) environment.
#'
#' @param curves_df A data frame with one row per \code{(curve_id, method)}
#'   combination.  Frequentist rows must carry raw \eqn{\log_{10}} parameter
#'   columns (\code{a}, \code{b}, \code{c}, \code{d}, \code{g});  Bayesian
#'   rows must also carry \code{*_nat} columns and a raw \code{c} column for
#'   the natural EC50.
#' @param n_coarse Integer; coarse-grid resolution passed to
#'   \code{\link{.compute_second_deriv}} (default \code{80L}).
#' @param n_refine Integer; refinement-grid resolution per bracket passed to
#'   \code{\link{.compute_second_deriv}} (default \code{80L}).
#' @param verbose Logical; if \code{TRUE} (default) emits a message when the
#'   derivative computation fails for a curve.
#'
#' @return A data frame (row-bound across all curves) with columns
#'   \code{curve_id}, \code{method}, \code{concentration} (natural scale), and
#'   \code{d2x_y} (\eqn{d^2(\log_{10} y)/d(\log_{10} x)^2}).  Rows with
#'   \code{NA} in \code{d2x_y} are dropped.
#'
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @export
compute_second_deriv_df <- function(curves_df,
                                    n_coarse = 80L,
                                    n_refine = 80L,
                                    verbose  = TRUE) {

  keys <- unique(curves_df[, c("curve_id", "method")])

  purrr::map_dfr(seq_len(nrow(keys)), function(i) {

    cid    <- keys$curve_id[[i]]
    method <- keys$method[[i]]
    row    <- curves_df[curves_df$curve_id == cid & curves_df$method == method, ]
    family <- .canonical_family(row$model_name)

    p <- if (method == "frequentist") {

      # Capture raw log10-scale params eagerly (guards against lazy-eval scoping)
      a_ <- as.numeric(row$a);  b_ <- as.numeric(row$b)
      c_ <- as.numeric(row$c);  d_ <- as.numeric(row$d)
      g_ <- if (is.na(row$g)) 1 else as.numeric(row$g)

      p_raw <- list(family = family, a = a_, b = b_, c = c_, d = d_, g = g_)

      list(
        family  = family,
        a       = 10^a_,
        b       = b_ * LN10,     # b_nat → correct grid span in log10-x
        c       = 10^c_,         # c_nat → grid centred at natural EC50
        d       = 10^d_,
        g       = g_,
        eval_fn = local({
          pr <- p_raw
          function(x) {
            y_nat <- 10^.evaluate_curve(pr, log10(x))
            log10(pmax(y_nat, 1e-10))   # d²(log10 y)/d(log10 x)²
          }
        })
      )

    } else {

      # Capture natural-scale params eagerly (guards against lazy-eval scoping)
      a_ <- as.numeric(row$a_nat);  b_ <- as.numeric(row$b_nat)
      c_ <- as.numeric(row$c_nat);  d_ <- as.numeric(row$d_nat)
      g_ <- as.numeric(row$g_nat);  c_raw_ <- as.numeric(row$c)

      p_nat <- list(family = family, a = a_, b = b_, c = c_, d = d_, g = g_)

      list(
        family  = family,
        a       = a_,
        b       = b_,             # b_nat = 1/b_raw → correct grid span
        c       = c_raw_,         # c_raw = natural EC50 → correct grid centre
        d       = d_,
        g       = g_,
        eval_fn = local({
          pn <- p_nat
          function(x) {
            y_nat <- .evaluate_curve(pn, log(x))
            log10(pmax(y_nat, 1e-10))   # d²(log10 y)/d(log10 x)²
          }
        })
      )
    }

    result <- tryCatch(
      .compute_second_deriv(p, n_coarse = n_coarse, n_refine = n_refine),
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_second_deriv_df] curve_id=%s (%s) deriv failed: %s",
          cid, method, e$message))
        tibble::tibble(concentration = NA_real_, d2x_y = NA_real_)
      }
    )

    result <- result[!is.na(result$d2x_y), ]
    result$curve_id <- cid
    result$method   <- method
    result[, c("curve_id", "method", "concentration", "d2x_y")]
  })
}
