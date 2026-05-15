# Estimate assay sensitivity as the slope at the inflection point

Computes the first derivative of the forward curve at the inflection
point using a central-difference approximation. A larger absolute value
indicates greater assay sensitivity (larger response change per unit
change in log10-concentration).

Estimates \\dy/dx\\ at `inflect_x` for each curve using a symmetric
central finite difference with step size \\h = \|x\| \times 10^{-5}\\.
Parameters are normalised via
[`.normalise_params`](https://immunoplex.github.io/curveRmetrics/reference/dot-normalise_params.md)
before evaluation.

## Usage

``` r
compute_assay_sensitivity(curves_df, verbose = TRUE)

compute_assay_sensitivity(curves_df, verbose = TRUE)
```

## Arguments

- curves_df:

  A data frame with one row per curve. Must contain `curve_id`,
  `method`, `model_name`, and `inflect_x` (natural-scale, as produced by
  [`compute_inflection_point`](https://immunoplex.github.io/curveRmetrics/reference/compute_inflection_point.md)),
  plus raw or `*_nat` parameter columns.

- verbose:

  Logical; if `TRUE` (default) emits a message per curve with
  `inflect_x` and the estimated slope.

## Value

Tibble with columns:

- `curve_id`:

  Curve identifier.

- `method`:

  `"frequentist"` or `"bayesian"`.

- `inflect_x`:

  Inflection-point concentration (log10 scale).

- `dydx_inflect`:

  First derivative of response with respect to log10-concentration at
  the inflection point.

A data frame with columns `curve_id`, `method`, `inflect_x`, and
`dydx_inflect`.

## Details

The inflection x-coordinate is taken directly from
`curves_df$inflect_x`; run compute_inflection_point() on the data frame
first.
