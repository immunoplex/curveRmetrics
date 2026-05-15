# Compute limits of quantification from the second derivative

Identifies the lower LOQ (`lloq`) and upper LOQ (`uloq`) as the
concentrations at the maximum and minimum of the second derivative of
the fitted curve, respectively. Local extrema are found via sign changes
in the first differences of the pre-computed `d2x_y` series; sub-grid
vertex positions are refined with a direct 3-point parabola
interpolation.

## Usage

``` r
compute_loqs(curves_df, second_deriv_df, verbose = TRUE)
```

## Arguments

- curves_df:

  A data frame with one row per curve (same requirements as
  [`compute_inflection_point`](https://immunoplex.github.io/curveRmetrics/reference/compute_inflection_point.md)),
  containing both raw and `*_nat` parameter columns.

- second_deriv_df:

  A data frame of pre-computed second derivative values with columns
  `curve_id`, `method`, `concentration` (natural scale), and `d2x_y`.

- verbose:

  Logical; if `TRUE` (default) emits a message per curve with the
  computed LOQ concentrations and response values.

## Value

A data frame with columns `curve_id`, `method`, `lloq`, `uloq`,
`lloq_y`, and `uloq_y`.

## Details

The response at each LOQ concentration is evaluated using the
method-appropriate evaluator:

- Frequentist: \\10^{f(\log\_{10}(x))}\\.

- Bayesian: \\f(\ln(x))\\.
