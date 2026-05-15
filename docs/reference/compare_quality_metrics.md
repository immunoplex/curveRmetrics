# Compare fitted curves and quality metrics for a single curve

Produces a log-log scatter/line plot showing:

- Observed standard points.

- Fitted frequentist and Bayesian curves evaluated over the observed
  concentration range.

- Vertical dashed lines at the LLOQ and ULOQ concentrations for each
  method.

- Horizontal lines at the LLOQ response, ULOQ response, LLOD, and ULOD
  for each method.

- A diamond point at the inflection coordinates for each method.

## Usage

``` r
compare_quality_metrics(
  curves_nat,
  standards,
  curves_quality_df,
  curve_id,
  x_var = "concentration",
  y_var = "absorbance"
)
```

## Arguments

- curves_nat:

  A data frame of fitted curve parameters (including `*_nat` columns for
  Bayesian rows) with at minimum `curve_id`, `method`, `model_name`,
  `a`, `b`, `c`, `d`, `g`, `a_nat`, `b_nat`, `c_nat`, `d_nat`, and
  `g_nat`.

- standards:

  A data frame of observed calibration standards with at minimum
  `curve_id`, a concentration column (named by `x_var`), and
  `assay_response`.

- curves_quality_df:

  A data frame as produced by
  [`attach_quality_metrics`](https://immunoplex.github.io/curveRmetrics/reference/attach_quality_metrics.md)
  containing `curve_id`, `method`, and all quality metric columns
  (`lloq`, `uloq`, `lloq_y`, `uloq_y`, `llod`, `ulod`, `inflect_x`,
  `inflect_y`).

- curve_id:

  Scalar character or integer; the identifier of the curve to plot.

- x_var:

  Character; name of the concentration column in `standards` (default
  `"concentration"`).

- y_var:

  Character; y-axis label string (default `"absorbance"`).

## Value

A `ggplot` object.

## Details

A subtitle in monospaced font summarises the key metric values
numerically for both methods.
