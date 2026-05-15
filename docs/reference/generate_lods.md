# Derive lower and upper limits of detection from asymptote CI bounds

Limits of detection are defined on the response scale:

- LLOD:

  Upper CI bound of the lower asymptote `a`.

- ULOD:

  Lower CI bound of the upper asymptote `d`.

## Usage

``` r
generate_lods(param_ci_df, verbose = FALSE)
```

## Arguments

- param_ci_df:

  A long-format data frame of parameter confidence intervals with
  columns `curve_id`, `method`, `parameter`, `estimate`, `conf_lower`,
  `conf_upper`, and `is_best_model`. May optionally contain `*_nat`
  variants.

- verbose:

  Logical; if `TRUE` emits a message per curve with the computed LLOD
  and ULOD.

## Value

A data frame with columns `curve_id`, `method`, `llod`, and `ulod`.

## Details

When the ULOD is less than the LLOD (indicating an implausible CI) it is
set to `NA`. Only rows flagged as `is_best_model` are processed.
