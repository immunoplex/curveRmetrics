#' ##' Compute Limits of  Detection
#'  #' @param param_ci_df  Long-format tibble (freq_param_ci_df / bayes_param_ci_df)
#'  #' @param  verbose Logical. If \code{TRUE}, emit one diagnostic message per curve.
#'  #' @return Tibble: curve_id, llod, ulod'
#'  #' @export
#' generate_lods <- function(param_ci_df, verbose = TRUE) {
#'
#'   best_df   <- param_ci_df[param_ci_df$is_best_model, ]
#'   curve_ids <- unique(best_df$curve_id)
#'
#'   results <- purrr::map(curve_ids, function(cid) {
#'     row  <- .pivot_curve_row(best_df[best_df$curve_id == cid, ])
#'
#'     ulod <- as.numeric(row$ci_lower["d"])
#'     llod <- as.numeric(row$ci_upper["a"])
#'
#'     # if (row$method == "bayesian") {
#'     #   ulod <- log10(ulod)
#'     #   llod <- log10(llod)
#'     # }
#'     ##|| ulod < 0 ||
#'     if (is.na(ulod) || is.na(llod) || ulod < llod) {
#'       if (verbose)
#'         message(sprintf(
#'           "[generate_lods] curve_id=%s  ULOD invalid (%.4f vs LLOD %.4f) - ULOD set to NA",
#'           cid, ulod, llod
#'         ))
#'       ulod <- NA_real_
#'     }
#'
#'     if (verbose)
#'       message(sprintf(
#'         "[generate_lods] curve_id=%s  LLOD=%.4f  ULOD=%s",
#'         cid, llod,
#'         if (is.na(ulod)) "NA" else sprintf("%.4f", ulod)
#'       ))
#'
#'     tibble::tibble(curve_id = cid, llod = llod, ulod = ulod)
#'   })
#'
#'   dplyr::bind_rows(results)
#' }
#'
#' ##' Compute minimum detectable concentrations and reliable detection limits
#' #'
#' #' Back-calculates concentration values from LOD response thresholds using
#' #' curve inversion:
#' #' \describe{
#' #'   \item{MDC (minimum detectable concentration)}{Inversion of the
#' #'     point-estimate curve at LLOD (`mindc`) and ULOD (`maxdc`).}
#' #'   \item{RDL (reliable detection limit)}{Inversion of CI-bound curves:
#' #'     `minrdl` uses a compressed curve (`a = ci_upper["a"]`,
#' #'     `d = ci_lower["d"]`); `maxrdl` uses an expanded curve
#' #'     (`a = ci_lower["a"]`, `d = ci_upper["d"]`).}
#' #' }
#' #'
#' #' For Bayesian curves the inversion result is in log10-concentration (matching
#' #' the normalised parameter scale); it is converted back to natural
#' #' concentration via \eqn{10^x} before storage so that downstream plotting
#' #' code can apply `log10()` uniformly.
#' #'
#' #'
#' #' @param param_ci_df Long-format tibble of parameter estimates and CIs,
#' #'   as passed to [generate_lods()].
#' #' @param lods Tibble returned by [generate_lods()], with columns
#' #'   `curve_id`, `llod`, and `ulod`.
#' #' @param verbose Logical. If `TRUE`, emits one diagnostic message per curve.
#' #'
#' #' @return Tibble with columns:
#' #'   \describe{
#' #'     \item{`curve_id`}{Curve identifier.}
#' #'     \item{`mindc`}{Lower MDC. Log10-concentration for frequentist;
#' #'       natural concentration for Bayesian.}
#' #'     \item{`maxdc`}{Upper MDC. Same scale convention as `mindc`.}
#' #'     \item{`minrdl`}{Lower RDL (conservative bound). Same scale as `mindc`.}
#' #'     \item{`maxrdl`}{Upper RDL (liberal bound). Same scale as `mindc`.}
#' #'   }
#' #'
#' #' @export
#' compute_mdc_rdl <- function(param_ci_df, lods, verbose = TRUE) {
#'
#'   best_df   <- param_ci_df[param_ci_df$is_best_model, ]
#'   curve_ids <- unique(best_df$curve_id)
#'
#'   results <- purrr::map(curve_ids, function(cid) {
#'
#'     row <- .pivot_curve_row(best_df[best_df$curve_id == cid, ])
#'     p   <- .normalise_params(row, row$method)
#'
#'     lod_row <- lods[lods$curve_id == cid, ]
#'     llod    <- as.numeric(lod_row$llod)
#'     ulod    <- as.numeric(lod_row$ulod)
#'
#'     # -- MDC: point-estimate curve --------------------------------------------
#'     mindc <- if (!is.na(llod)) tryCatch(.invert_curve(p, llod), error = function(e) NA_real_) else NA_real_
#'     maxdc <- if (!is.na(ulod)) tryCatch(.invert_curve(p, ulod), error = function(e) NA_real_) else NA_real_
#'
#'     # -- RDL: CI-bound curves ------------------------------------------------
#'     p_lo   <- modifyList(p, list(d = as.numeric(row$ci_lower["d"])))
#'     p_hi   <- modifyList(p, list(#a = as.numeric(row$ci_lower["a"]),
#'       d = as.numeric(row$ci_upper["d"])))
#'
#'     minrdl <- if (!is.na(llod)) tryCatch(.invert_curve(p_lo, llod), error = function(e) NA_real_) else NA_real_
#'     maxrdl <- if (!is.na(ulod)) tryCatch(.invert_curve(p_hi, ulod), error = function(e) NA_real_) else NA_real_
#'
#'     if (verbose)
#'       message(sprintf(
#'         "[compute_mdc_rdl] curve_id=%s  mindc=%s  maxdc=%s  minrdl=%s  maxrdl=%s",
#'         cid,
#'         format(mindc,  digits = 4), format(maxdc,  digits = 4),
#'         format(minrdl, digits = 4), format(maxrdl, digits = 4)
#'       ))
#'
#'     tibble::tibble(
#'       curve_id = cid,
#'       mindc    = as.numeric(mindc),
#'       maxdc    = as.numeric(maxdc),
#'       minrdl   = as.numeric(minrdl),
#'       maxrdl   = as.numeric(maxrdl)
#'     )
#'   })
#'
#'   dplyr::bind_rows(results)
#' }
#'
