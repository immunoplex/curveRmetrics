# Normalise a single DB row to a natural-scale parameter list

Prefers pre-computed `*_nat` columns produced by
[`transform_to_natural_units`](https://immunoplex.github.io/curveRmetrics/reference/transform_to_natural_units.md)
when they are present on `row` (i.e. when called with `curves_nat`).
Falls back to raw DB columns plus inline conversion otherwise.

## Usage

``` r
.normalise_params(row, method, verbose = FALSE)
```

## Arguments

- row:

  A single-row data frame (or named list) with curve parameters.

- method:

  Character; either `"frequentist"` or `"bayesian"`.

- verbose:

  Logical; if `TRUE` (default `FALSE`) emits a message when
  natural-scale columns are detected and used.

## Value

A named list with elements `family`, `a`, `b`, `c`, `d`, `g`, and
`source`.

## Details

Either way the returned list has an identical structure and scale:

- a, d:

  Natural-scale response asymptotes.

- b:

  Hill slope on the \\\ln(x)\\ scale.

- c:

  Natural-scale EC50 (concentration units).

- g:

  Dimensionless asymmetry/shape parameter.
