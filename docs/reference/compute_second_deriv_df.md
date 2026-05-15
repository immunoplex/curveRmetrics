# Compute the second derivative series for all curves in a data frame

Public entry point. Iterates over unique `(curve_id, method)`
combinations and calls
[`.compute_second_deriv`](https://immunoplex.github.io/curveRmetrics/reference/dot-compute_second_deriv.md)
for each, returning a long-format data frame of second derivative values
suitable for use in
[`compute_loqs`](https://immunoplex.github.io/curveRmetrics/reference/compute_loqs.md).

## Usage

``` r
compute_second_deriv_df(
  curves_df,
  n_coarse = 80L,
  n_refine = 80L,
  verbose = TRUE
)
```

## Arguments

- curves_df:

  A data frame with one row per `(curve_id, method)` combination.
  Frequentist rows must carry raw \\\log\_{10}\\ parameter columns (`a`,
  `b`, `c`, `d`, `g`); Bayesian rows must also carry `*_nat` columns and
  a raw `c` column for the natural EC50.

- n_coarse:

  Integer; coarse-grid resolution passed to
  [`.compute_second_deriv`](https://immunoplex.github.io/curveRmetrics/reference/dot-compute_second_deriv.md)
  (default `80L`).

- n_refine:

  Integer; refinement-grid resolution per bracket passed to
  [`.compute_second_deriv`](https://immunoplex.github.io/curveRmetrics/reference/dot-compute_second_deriv.md)
  (default `80L`).

- verbose:

  Logical; if `TRUE` (default) emits a message when the derivative
  computation fails for a curve.

## Value

A data frame (row-bound across all curves) with columns `curve_id`,
`method`, `concentration` (natural scale), and `d2x_y` (\\d^2(\log\_{10}
y)/d(\log\_{10} x)^2\\). Rows with `NA` in `d2x_y` are dropped.

## Details

A method-specific parameter list is built for each curve:

- Frequentist:

  Grid parameters use \\b\_{nat} = b\_{raw} \cdot \ln(10)\\ and
  \\c\_{nat} = 10^{c\_{raw}}\\. The `eval_fn` closure computes
  \\d^2(\log\_{10} y)/d(\log\_{10} x)^2\\ via
  `log10(10^.evaluate_curve(p_raw, log10(x)))`.

- Bayesian:

  Grid parameters use \\b\_{nat} = 1/b\_{raw}\\ and the raw natural EC50
  \\c\_{raw}\\. The `eval_fn` closure computes \\d^2(\log\_{10}
  y)/d(\log\_{10} x)^2\\ via `log10(.evaluate_curve(p_nat, log(x)))`.

[`local()`](https://rdrr.io/r/base/eval.html) is used when constructing
each `eval_fn` to force eager capture of `p_raw` / `p_nat`, preventing
R's lazy-evaluation scoping from causing all closures to share the same
(final-iteration) environment.
