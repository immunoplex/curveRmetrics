# Vectorised curve evaluation with per-point error guard

Attempts to evaluate `p$eval_fn` (or
[`.evaluate_curve`](https://immunoplex.github.io/curveRmetrics/reference/dot-evaluate_curve.md)
as a fallback) over the entire vector `x_vec` in a single call. If that
call throws an error, each point is evaluated individually and failures
are replaced with `NA`.

## Usage

``` r
.sderiv_eval_vec(p, x_vec)
```

## Arguments

- p:

  A named list of curve parameters. If `p$eval_fn` is a function it is
  used as the evaluator; otherwise
  [`.evaluate_curve`](https://immunoplex.github.io/curveRmetrics/reference/dot-evaluate_curve.md)
  is called with `p` and `x`.

- x_vec:

  Numeric vector of concentration values to evaluate.

## Value

Numeric vector of the same length as `x_vec` with `NA` substituted for
any evaluation failure.
