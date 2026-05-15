# Evaluate a dose-response curve using method-appropriate scaling

A thin wrapper around
[`.evaluate_curve`](https://immunoplex.github.io/curveRmetrics/reference/dot-evaluate_curve.md)
that applies the correct log transformation depending on whether the
model was fit by the frequentist or Bayesian method.

## Usage

``` r
.evaluate_curve_auto(p, x, method)
```

## Arguments

- p:

  A named list of curve parameters (see
  [`.evaluate_curve`](https://immunoplex.github.io/curveRmetrics/reference/dot-evaluate_curve.md)).

- x:

  Numeric; natural-scale concentration value(s).

- method:

  Character; `"frequentist"` or `"bayesian"`.

## Value

Numeric; predicted response on the natural scale.

## Details

- Frequentist:

  Parameters are on the \\\log\_{10}\\ scale. `x` is converted to
  \\\log\_{10}(x)\\, the curve is evaluated, and the result is
  back-transformed via \\10^y\\.

- Bayesian:

  Parameters and `x` are passed as \\\ln(x)\\ to match the Stan model's
  internal parameterisation.
