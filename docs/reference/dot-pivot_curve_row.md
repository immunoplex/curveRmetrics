# Pivot a long-format parameter CI data frame to a wide named list

Selects rows for parameters `a`, `b`, `c`, `d`, and (optionally) `g`,
then pivots them into named vectors of point estimates and confidence
bounds.

## Usage

``` r
.pivot_curve_row(curve_rows, verbose = FALSE)
```

## Arguments

- curve_rows:

  A data frame of CI rows for a single curve, containing at minimum
  `parameter`, `estimate`, `conf_lower`, `conf_upper`, `curve_id`,
  `method`, and `model_name`. May additionally contain `estimate_nat`,
  `conf_lower_nat`, and `conf_upper_nat`.

- verbose:

  Logical; if `TRUE` (default `FALSE`) emits a message indicating
  whether natural-scale columns were found.

## Value

A named list with elements `curve_id`, `method`, `model_name`, `a`, `b`,
`c`, `d`, `g`, `ci_lower` (named numeric vector), and `ci_upper` (named
numeric vector).

## Details

When the column `estimate_nat` is present on `curve_rows` the
natural-scale estimates and CI bounds are used; otherwise the raw DB
columns `estimate`, `conf_lower`, and `conf_upper` are used.
