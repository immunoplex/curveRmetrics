# Three-pass adaptive numerical second derivative for a single curve

Computes \\d^2y/du^2\\ (\\u = \log\_{10}(x)\\) over a family-derived
concentration grid using three passes:

## Usage

``` r
.compute_second_deriv(p, n_coarse = 80L, n_refine = 80L)
```

## Arguments

- p:

  A named list of natural-scale curve parameters including a `family`
  string, `b`, `c`, `g`, and an `eval_fn` closure that maps
  concentration to response (see
  [`compute_second_deriv_df`](https://immunoplex.github.io/curveRmetrics/reference/compute_second_deriv_df.md)).

- n_coarse:

  Integer; number of points in the initial coarse grid (default `80L`).

- n_refine:

  Integer; number of points added inside each refinement bracket
  (default `80L`).

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns `concentration` (natural scale) and `d2x_y` (\\d^2y/du^2\\).

## Details

1.  **Coarse pass** — `n_coarse` evenly spaced \\\log\_{10}(x)\\ points
    across the family bounds from
    [`.sderiv_scale`](https://immunoplex.github.io/curveRmetrics/reference/dot-sderiv_scale.md).

2.  **Seed brackets** —
    [`.sderiv_seeds`](https://immunoplex.github.io/curveRmetrics/reference/dot-sderiv_seeds.md)
    provides analytical estimates of the extremum locations; each seed
    is bracketed by ±4 coarse grid spacings.

3.  **Refinement** — `n_refine` dense points are placed inside every
    bracket detected from both the coarse pass sign-changes
    ([`.sderiv_brackets`](https://immunoplex.github.io/curveRmetrics/reference/dot-sderiv_brackets.md))
    and the seed brackets.

All three pass grids are merged, deduplicated, and evaluated in a single
call to
[`.sderiv_eval_vec`](https://immunoplex.github.io/curveRmetrics/reference/dot-sderiv_eval_vec.md),
after which
[`.d2y_nonuniform`](https://immunoplex.github.io/curveRmetrics/reference/dot-d2y_nonuniform.md)
computes the final derivative series.
