utils::globalVariables(c(
  # ggplot2 aesthetics in compare_quality_metrics
  "concentration", "y", "method", "val", "metric", "x",
  # compare_second_derivative
  "d2x_y",
  # transform_ci_to_natural_units
  "estimate", "parameter", "conf_lower", "conf_upper",
  "lower_nat", "upper_nat",
  # transform_to_natural_units
  "model_name",
  # magrittr
  "."
))

#' @importFrom dplyr bind_rows filter
#' @importFrom purrr map map_dfr
#' @importFrom tibble tibble as_tibble
#' @importFrom stats setNames
#' @importFrom utils modifyList
#' @import ggplot2
#' @importFrom magrittr %>%
"_PACKAGE"
