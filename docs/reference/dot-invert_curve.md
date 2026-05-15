# Back-calculate natural concentration from a response value

Analytically inverts each supported dose-response family. All parameters
and the returned concentration are on the natural scale.

## Usage

``` r
.invert_curve(p, y)
```

## Arguments

- p:

  A named list of natural-scale curve parameters (see
  [`.evaluate_curve`](https://immunoplex.github.io/curveRmetrics/reference/dot-evaluate_curve.md)).

- y:

  Numeric; response value(s) to invert. Must lie within the range
  \\\[\min(a, d),\\ \max(a, d)\]\\.

## Value

Numeric; natural-scale concentration(s) corresponding to `y`.
