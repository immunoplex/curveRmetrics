# =============================================================================
# natural_units.R
#
# Conversion rules (frequentist only; Bayesian already on natural scale):
#   a, d  : 10^x        (log10 response  -> natural response)
#   c     : 10^x        (log10 EC50      -> natural EC50)
#   b     : x * ln(10)  (log10 slope     -> ln slope)
#   g     : unchanged   (dimensionless)
# =============================================================================

LN10 <- log(10)   # ≈ 2.302585


#' Convert a single parameter value to the natural scale
#'
#' Low-level scalar/vector converter applied element-wise by
#' \code{\link{transform_ci_to_natural_units}}.  Conversion rules differ by
#' parameter and fitting method:
#'
#' \tabular{lll}{
#'   \strong{Parameter} \tab \strong{Frequentist} \tab \strong{Bayesian} \cr
#'   \code{a}, \code{d} \tab \eqn{10^x}           \tab unchanged         \cr
#'   \code{c}           \tab \eqn{10^x}            \tab \eqn{\ln(x)}      \cr
#'   \code{b}           \tab \eqn{x \cdot \ln 10}  \tab \eqn{1/x}         \cr
#'   \code{g}           \tab unchanged              \tab unchanged         \cr
#' }
#'
#' Any unrecognised parameter name is passed through unchanged.
#'
#' @param x Numeric value(s) to convert.
#' @param param Character; parameter name — one of \code{"a"}, \code{"b"},
#'   \code{"c"}, \code{"d"}, or \code{"g"}.
#' @param method Character; fitting method — \code{"frequentist"} or
#'   \code{"bayesian"}.
#'
#' @return Numeric value(s) on the natural scale.
#' @keywords internal
.to_nat <- function(x, param, method) {
  switch(param,
         a = , d = switch(method, frequentist = 10^x,      bayesian = x),
         c      = switch(method, frequentist = 10^x,       bayesian = log(x)),
         b      = switch(method, frequentist = x * LN10,   bayesian = 1 / x),
         g      = x,   # dimensionless — unchanged for both methods
         x             # unknown params passed through
  )
}


#' Transform long-format parameter CI data to the natural scale
#'
#' Adds natural-scale counterparts of the estimate and confidence bound
#' columns to a long-format parameter CI data frame, suitable for use with
#' both \code{freq_param_ci_df} and \code{bayes_param_ci_df} (or their
#' row-bound combination).
#'
#' Each value is converted via \code{\link{.to_nat}} according to its
#' \code{parameter} and \code{method} fields.  Because \eqn{10^x} and
#' \eqn{x \cdot \ln(10)} are both monotone increasing, CI ordering is
#' preserved for \code{a}, \code{c}, \code{d}, and \code{b} under the
#' frequentist method.  For the Bayesian \code{b} conversion (\eqn{1/x}),
#' which is monotone \emph{decreasing}, the lower and upper bounds are
#' swapped after transformation so that \code{conf_lower_nat} always holds
#' the smaller value.
#'
#' @param param_ci_df A long-format data frame with at minimum the columns
#'   \code{parameter}, \code{method}, \code{estimate}, \code{conf_lower},
#'   and \code{conf_upper}.
#'
#' @return \code{param_ci_df} with three additional columns:
#' \describe{
#'   \item{estimate_nat}{Point estimate on the natural scale.}
#'   \item{conf_lower_nat}{Lower confidence bound on the natural scale
#'     (\code{min} of the transformed bounds).}
#'   \item{conf_upper_nat}{Upper confidence bound on the natural scale
#'     (\code{max} of the transformed bounds).}
#' }
#'
#' @importFrom dplyr mutate select
#' @export
transform_ci_to_natural_units <- function(param_ci_df) {
  param_ci_df %>%
    mutate(
      estimate_nat   = mapply(.to_nat, estimate,   parameter, method),
      lower_nat      = mapply(.to_nat, conf_lower, parameter, method),
      upper_nat      = mapply(.to_nat, conf_upper, parameter, method),
      conf_lower_nat = pmin(lower_nat, upper_nat),
      conf_upper_nat = pmax(lower_nat, upper_nat)
    ) %>%
    select(-lower_nat, -upper_nat)
}
