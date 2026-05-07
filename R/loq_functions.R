#' Compute limits of quantification from the second derivative of fitted curves
#'
#' Identifies the lower LOQ (LLOQ) and upper LOQ (ULOQ) as the concentration
#' coordinates of the outermost local extrema of the second derivative of the
#' fitted curve, supplied as a pre-computed data frame.  Sub-grid vertex
#' positions are refined by 3-point parabolic interpolation.
#'
#' The second derivative extrema mark the boundaries of the curve's linear
#' (maximum-slope) region: the lower extremum (maximum of \eqn{d^2y/dx^2})
#' gives the LLOQ and the upper extremum (minimum of \eqn{d^2y/dx^2}) gives
#' the ULOQ.
#'
#' @param curves_df  Data frame of fitted curves, e.g. `freq_curves_df` or
#'   `bayes_curves_df`.
#' @param second_deriv_df Data frame with columns `curve_id`, `concentration`
#'   (on the same scale as `c` for that method), and `d2x_y` (second
#'   derivative of response with respect to concentration).
#' @param verbose Logical. If `TRUE`, emits one diagnostic message per curve.
#'
#' @return Tibble with columns:
#'   \describe{
#'     \item{`curve_id`}{Curve identifier.}
#'     \item{`method`}{`"frequentist"` or `"bayesian"`.}
#'     \item{`lloq`}{Lower LOQ concentration.}
#'     \item{`uloq`}{Upper LOQ concentration.}
#'     \item{`lloq_y`}{Predicted response at the LLOQ concentration.}
#'     \item{`uloq_y`}{Predicted response at the ULOQ concentration.}
#'   }
#'
#' @export
compute_loqs <- function(curves_df, second_deriv_df, verbose = TRUE) {

  results <- purrr::map(seq_len(nrow(curves_df)), function(i) {

    row <- curves_df[i, ]
    cid <- row$curve_id
    p   <- .normalise_params(row, row$method)

    # -- pull second derivative data for this curve -----------------
    sd_sub <- second_deriv_df[second_deriv_df$curve_id == cid, ]

    if (nrow(sd_sub) < 3) {
      if (verbose) message(sprintf(
        "[compute_loqs] curve_id=%s: fewer than 3 second deriv points, skipping", cid))
      return(tibble::tibble(
        curve_id = cid, method = row$method,
        lloq = NA_real_, uloq = NA_real_,
        lloq_y = NA_real_, uloq_y = NA_real_
      ))
    }

    x <- as.numeric(sd_sub$concentration)
    y <- as.numeric(sd_sub$d2x_y)

    # ── local extrema via first differences ──────────────────────────────────
    dy      <- diff(y)
    idx_max <- which(dy[-1] < 0 & dy[-length(dy)] > 0) + 1
    idx_min <- which(dy[-1] > 0 & dy[-length(dy)] < 0) + 1

    # ── sub-grid vertex via direct 3-point parabola ───────────────────────────
    # Avoids solve(t(X)%*%X) which squares the condition number and fails
    # when adjacent x points are nearly coincident.
    interpolate_vertex <- function(idx) {
      xi <- x[(idx - 1):(idx + 1)]
      yi <- y[(idx - 1):(idx + 1)]

      x1 <- xi[1]; x2 <- xi[2]; x3 <- xi[3]
      y1 <- yi[1]; y2 <- yi[2]; y3 <- yi[3]

      denom <- (x1 - x2) * (x1 - x3) * (x2 - x3)
      if (abs(denom) < .Machine$double.eps * 1e8)
        return(list(x = x2, y = y2))

      a <- (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
      b <- (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3)) / denom

      if (abs(a) < .Machine$double.eps)
        return(list(x = x2, y = y2))

      xv <- -b / (2 * a)
      xv <- max(min(xv, max(xi)), min(xi))
      yv <- a * xv^2 + b * xv + (y1 - a * x1^2 - b * x1)

      list(x = as.numeric(xv), y = as.numeric(yv))
    }

    # ── collect all local maxima / minima ─────────────────────────────────────
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

    # ── pick global extrema ───────────────────────────────────────────────────
    lloq_x <- if (nrow(max_df) > 0) max_df$x[which.max(max_df$y)] else NA_real_
    uloq_x <- if (nrow(min_df) > 0) min_df$x[which.min(min_df$y)] else NA_real_

    # ── evaluate forward curve at LOQ positions ───────────────────────────────
    lloq_y <- tryCatch(
      if (!is.na(lloq_x)) as.numeric(.evaluate_curve(p, lloq_x)) else NA_real_,
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_loqs] curve_id=%s lloq_y eval failed: %s", cid, e$message))
        NA_real_
      }
    )

    uloq_y <- tryCatch(
      if (!is.na(uloq_x)) as.numeric(.evaluate_curve(p, uloq_x)) else NA_real_,
      error = function(e) {
        if (verbose) message(sprintf(
          "[compute_loqs] curve_id=%s uloq_y eval failed: %s", cid, e$message))
        NA_real_
      }
    )

    if (verbose) message(sprintf(
      "[compute_loqs] curve_id=%s (%s)  lloq_x=%.4f  uloq_x=%.4f  lloq_y=%.4f  uloq_y=%.4f",
      cid, row$method, lloq_x, uloq_x, lloq_y, uloq_y))

    tibble::tibble(
      curve_id = cid,
      method   = row$method,
      lloq     = as.numeric(lloq_x),
      uloq     = as.numeric(uloq_x),
      lloq_y   = as.numeric(lloq_y),
      uloq_y   = as.numeric(uloq_y)
    )
  })

  dplyr::bind_rows(results)
}
