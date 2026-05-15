#' Plot the second derivative of a fitted curve for both methods
#'
#' Renders a log-scale line plot of \eqn{d^2y/dx^2} against concentration for
#' a single curve, overlaying the frequentist and Bayesian traces.  The x-axis
#' is automatically zoomed to the region where the second derivative is
#' meaningfully non-zero (greater than 1 % of its maximum absolute value).
#'
#' @param second_derivative_df A data frame with at minimum the columns
#'   \code{curve_id}, \code{method}, \code{concentration}, and \code{d2x_y}
#'   (the second derivative of the response with respect to concentration).
#' @param curve_id Scalar character or integer; the identifier of the curve to
#'   plot.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_color_manual
#'   scale_x_log10 coord_cartesian labs theme_bw
#' @export
compare_second_derivative <- function(second_derivative_df, curve_id) {

  colors <- list("frequentist" = "#0067a5", "bayesian" = "#f38400")

  d2_sub <- second_derivative_df[second_derivative_df$curve_id == curve_id, ]

  # Zoom to the x range where d2y is meaningfully non-zero (> 1 % of max abs)
  nonzero <- d2_sub[abs(d2_sub$d2x_y) > max(abs(d2_sub$d2x_y), na.rm = TRUE) * 0.01, ]
  x_lo    <- min(nonzero$concentration) * 0.1
  x_hi    <- max(nonzero$concentration) * 10

  ggplot(d2_sub, aes(x = concentration, y = d2x_y, color = method)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_color_manual(values = unlist(colors)) +
    scale_x_log10() +
    coord_cartesian(xlim = c(x_lo, x_hi)) +
    labs(
      title = sprintf("Second Derivative Comparison for curve %s", curve_id),
      x     = "Concentration",
      y     = expression(d^2 * y / dx^2)
    ) +
    theme_bw()
}


#' Attach quality metrics to a curves data frame
#'
#' Left-joins inflection points, assay sensitivity, limits of detection,
#' reliable detection limits, and limits of quantification onto the curves data
#' frame.  The \code{model_name} column is first resolved to its canonical
#' family name via \code{\link{.canonical_family}}.
#'
#' @param curves_df A data frame with one row per \code{(curve_id, method)}
#'   combination, containing at minimum a \code{model_name} column.
#' @param lods A data frame as returned by \code{\link{generate_lods}} with
#'   columns \code{curve_id}, \code{method}, \code{llod}, and \code{ulod}.
#' @param rdls A data frame of reliable detection limits with columns
#'   \code{curve_id}, \code{method}, \code{minrdl}, and \code{maxrdl}.
#' @param sensitivity A data frame as returned by
#'   \code{\link{compute_assay_sensitivity}} with columns \code{curve_id},
#'   \code{method}, and \code{dydx_inflect} (the \code{inflect_x} column, if
#'   present, is dropped before merging to avoid duplication).
#' @param loqs A data frame as returned by \code{\link{compute_loqs}} with
#'   columns \code{curve_id}, \code{method}, \code{lloq}, \code{uloq},
#'   \code{lloq_y}, and \code{uloq_y}.
#' @param inflection A data frame as returned by
#'   \code{\link{compute_inflection_point}} with columns \code{curve_id},
#'   \code{method}, \code{inflect_x}, and \code{inflect_y}.
#'
#' @return \code{curves_df} with additional columns from each quality-metric
#'   data frame, joined on \code{(curve_id, method)}.
#'
#' @export
attach_quality_metrics <- function(curves_df, lods, rdls, sensitivity, loqs, inflection) {

  curves_df$model_name <- vapply(curves_df$model_name, .canonical_family, character(1))

  curves_df <- merge(
    curves_df,
    inflection[, c("curve_id", "method", "inflect_x", "inflect_y")],
    by    = c("curve_id", "method"),
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    sensitivity[, !names(sensitivity) %in% "inflect_x"],
    by    = c("curve_id", "method"),
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    lods,
    by    = c("curve_id", "method"),
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    rdls,
    by    = c("curve_id", "method"),
    all.x = TRUE
  )

  curves_df <- merge(
    curves_df,
    loqs,
    by    = c("curve_id", "method"),
    all.x = TRUE
  )

  curves_df
}


#' Compare fitted curves and quality metrics for a single curve
#'
#' Produces a log-log scatter/line plot showing:
#' \itemize{
#'   \item Observed standard points.
#'   \item Fitted frequentist and Bayesian curves evaluated over the observed
#'     concentration range.
#'   \item Vertical dashed lines at the LLOQ and ULOQ concentrations for each
#'     method.
#'   \item Horizontal lines at the LLOQ response, ULOQ response, LLOD, and
#'     ULOD for each method.
#'   \item A diamond point at the inflection coordinates for each method.
#' }
#'
#' A subtitle in monospaced font summarises the key metric values numerically
#' for both methods.
#'
#' @param curves_nat A data frame of fitted curve parameters (including
#'   \code{*_nat} columns for Bayesian rows) with at minimum \code{curve_id},
#'   \code{method}, \code{model_name}, \code{a}, \code{b}, \code{c}, \code{d},
#'   \code{g}, \code{a_nat}, \code{b_nat}, \code{c_nat}, \code{d_nat}, and
#'   \code{g_nat}.
#' @param standards A data frame of observed calibration standards with at
#'   minimum \code{curve_id}, a concentration column (named by \code{x_var}),
#'   and \code{assay_response}.
#' @param curves_quality_df A data frame as produced by
#'   \code{\link{attach_quality_metrics}} containing \code{curve_id},
#'   \code{method}, and all quality metric columns (\code{lloq}, \code{uloq},
#'   \code{lloq_y}, \code{uloq_y}, \code{llod}, \code{ulod}, \code{inflect_x},
#'   \code{inflect_y}).
#' @param curve_id Scalar character or integer; the identifier of the curve to
#'   plot.
#' @param x_var Character; name of the concentration column in \code{standards}
#'   (default \code{"concentration"}).
#' @param y_var Character; y-axis label string (default \code{"absorbance"}).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_vline geom_hline
#'   scale_colour_manual scale_linetype_manual scale_x_log10 scale_y_log10
#'   labs theme_bw theme element_text .data
#' @importFrom purrr map_dfr
#' @export
compare_quality_metrics <- function(curves_nat, standards, curves_quality_df, curve_id,
                                 x_var = "concentration", y_var = "absorbance") {

  curve_nat_df <- curves_nat[curves_nat$curve_id == curve_id, ]
  standards_df <- standards[standards$curve_id == curve_id, ]
  quality_df   <- curves_quality_df[curves_quality_df$curve_id == curve_id, ]

  # ── Fitted curves ──────────────────────────────────────────────────────────
  conc_min <- min(standards_df[[x_var]][standards_df[[x_var]] > 0])
  conc_seq <- 10^seq(log10(conc_min), log10(max(standards_df[[x_var]])), length.out = 400)

  make_param_list <- function(method) {
    p <- curve_nat_df[curve_nat_df$method == method, ]
    if (method == "frequentist")
      list(family = p$model_name, a = p$a,     b = p$b,     c = p$c,     d = p$d,     g = p$g,     method = method)
    else
      list(family = p$model_name, a = p$a_nat, b = p$b_nat, c = p$c_nat, d = p$d_nat, g = p$g_nat, method = method)
  }

  curve_df <- data.frame(
    concentration = rep(conc_seq, 2),
    y      = c(
      .evaluate_curve_auto(make_param_list("frequentist"), conc_seq, "frequentist"),
      .evaluate_curve_auto(make_param_list("bayesian"),    conc_seq, "bayesian")
    ),
    method = rep(c("frequentist", "bayesian"), each = length(conc_seq))
  )

  method_colors <- c("bayesian" = "#f38400", "frequentist" = "#0067a5")

  # ── Helper: build long-format line data from a list of metric specs ────────
  # Each spec is a list(col, label, linetype); rows with NA or non-positive
  # values are silently dropped.
  make_lines <- function(specs) {
    purrr::map_dfr(specs, function(m) {
      purrr::map_dfr(c("bayesian", "frequentist"), function(meth) {
        row <- quality_df[quality_df$method == meth, ]
        val <- if (m$col %in% names(row)) as.numeric(row[[m$col]]) else NA_real_
        if (!is.na(val) && val > 0)
          data.frame(val = val, method = meth, metric = m$label,
                     linetype = m$linetype, stringsAsFactors = FALSE)
      })
    })
  }

  # ── Subtitle: one line per method summarising key metric values ────────────
  make_metric_label <- function(qdf, meth) {
    r <- qdf[qdf$method == meth, ]
    sprintf(
      "%s: LLOQ=%.1f (y=%.3f)  ULOQ=%.1f (y=%.3f)  LLOD=%.3f  ULOD=%.3f",
      ifelse(meth == "bayesian", "Bayesian", "Frequentist"),
      r$lloq, r$lloq_y, r$uloq, r$uloq_y, r$llod, r$ulod
    )
  }

  subtitle_text <- paste(
    make_metric_label(quality_df, "frequentist"),
    make_metric_label(quality_df, "bayesian"),
    sep = "\n"
  )

  # ── Line data ──────────────────────────────────────────────────────────────
  vlines_df <- make_lines(list(
    list(col = "lloq", label = "LLOQ", linetype = "dashed"),
    list(col = "uloq", label = "ULOQ", linetype = "dashed")
  ))

  hlines_df <- make_lines(list(
    list(col = "lloq_y", label = "LLOQ", linetype = "dashed"),
    list(col = "uloq_y", label = "ULOQ", linetype = "dashed"),
    list(col = "llod",   label = "LLOD", linetype = "dotted"),
    list(col = "ulod",   label = "ULOD", linetype = "dotted")
  ))

  # ── Inflection point markers (one diamond per method) ──────────────────────
  inflect_df <- purrr::map_dfr(c("bayesian", "frequentist"), function(meth) {
    row <- quality_df[quality_df$method == meth, ]
    xv  <- as.numeric(row$inflect_x)
    yv  <- as.numeric(row$inflect_y)
    if (!is.na(xv) && !is.na(yv) && xv > 0 && yv > 0)
      data.frame(x = xv, y = yv, method = meth, stringsAsFactors = FALSE)
  })

  # ── Shared linetype scale ──────────────────────────────────────────────────
  all_linetypes <- c(
    "LLOQ"    = "dashed",
    "ULOQ"    = "dashed",
    "LLOD"    = "dotted",
    "ULOD"    = "dotted",
    "Min DC"  = "longdash",
    "Max DC"  = "longdash",
    "Min RDL" = "twodash",
    "Max RDL" = "twodash"
  )

  ggplot() +
    geom_point(data = standards_df,
               aes(x = .data[[x_var]], y = .data[["assay_response"]]),
               colour = "grey30", alpha = 0.6, size = 2) +
    geom_line(data = curve_df,
              aes(x = concentration, y = y, colour = method),
              linewidth = 0.9) +
    geom_vline(data = vlines_df,
               aes(xintercept = val, colour = method, linetype = metric),
               linewidth = 0.6, alpha = 0.8) +
    geom_hline(data = hlines_df,
               aes(yintercept = val, colour = method, linetype = metric),
               linewidth = 0.6, alpha = 0.8) +
    geom_point(data = inflect_df,
               aes(x = x, y = y, colour = method),
               shape = 18, size = 5) +
    scale_colour_manual(
      name   = "Method",
      values = method_colors,
      labels = c("bayesian" = "Bayesian", "frequentist" = "Frequentist")
    ) +
    scale_linetype_manual(name = "Metric", values = all_linetypes) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title    = sprintf("Quality Metric Comparison for Curve %s", curve_id),
      subtitle = subtitle_text,
      x        = x_var,
      y        = y_var
    ) +
    theme_bw() +
    theme(
      plot.subtitle   = element_text(size = 7, family = "mono", colour = "grey30"),
      legend.position = "right",
      legend.box      = "vertical"
    )
}
