# Add inflection-point coordinates to a curves data frame

Applies compute_inflection_point_worker() to every row of `curves_df`
and appends two new columns: `inflect_x` (concentration at the
inflection point) and `inflect_y` (predicted response at the inflection
point). Both values are on the same scale as `c` after normalisation by
.normalise_params().

Iterates over unique `(curve_id, method)` combinations, resolves
parameters via
[`.inflect_x_analytical`](https://immunoplex.github.io/curveRmetrics/reference/dot-inflect_x_analytical.md),
evaluates the response at that concentration, and joins the results back
onto the input data frame.

## Usage

``` r
compute_inflection_point(curves_df, verbose = TRUE)

compute_inflection_point(curves_df, verbose = TRUE)
```

## Arguments

- curves_df:

  A data frame with one row per curve containing at minimum `curve_id`,
  `method`, `model_name`, `a`, `b`, `c`, `d`, and `g`. Frequentist rows
  hold raw \\\log\_{10}\\ parameters; Bayesian rows must also carry
  `*_nat` columns.

- verbose:

  Logical; if `TRUE` (default) emits a message per curve with the
  computed inflection coordinates.

## Value

`curves_df` with two additional numeric columns: `inflect_x` and
`inflect_y`.

`curves_df` with two new (or replaced) columns: `inflect_x` and
`inflect_y`, both on the natural scale.

## Details

The inflection point is defined as the concentration at the midpoint of
the response range (\\y^\* = \sqrt{a \cdot d}\\), computed analytically
without root-finding.
