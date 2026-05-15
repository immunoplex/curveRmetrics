# Compute the inflection point for a single curve row (worker)

Normalises parameters via
[`.normalise_params`](https://immunoplex.github.io/curveRmetrics/reference/dot-normalise_params.md),
locates the inflection concentration with
[`.inflection_x`](https://immunoplex.github.io/curveRmetrics/reference/dot-inflection_x.md),
and evaluates the corresponding response with
[`.evaluate_curve`](https://immunoplex.github.io/curveRmetrics/reference/dot-evaluate_curve.md).

## Usage

``` r
compute_inflection_point_worker(row, verbose = TRUE)
```

## Arguments

- row:

  A single-row data frame with curve parameters and a `method` column.

- verbose:

  Logical; if `TRUE` (default) emits diagnostic messages with the
  resolved parameters and computed inflection coordinates.

## Value

A named list with elements:

- inflect_x:

  Natural-scale concentration at the inflection point.

- inflect_y:

  Predicted response at the inflection point.

- family:

  Canonical family name.

- source:

  Fitting method (`"frequentist"` or `"bayesian"`).

- params:

  The normalised parameter list passed to the worker functions.
