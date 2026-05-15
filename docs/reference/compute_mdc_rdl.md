# Compute the Minimum Detectable Concentration and Reliable Detection Limits

Derives MDC and RDL by inverting the dose-response curve at the LOD
response thresholds:

- MDC:

  `mindc` = concentration at LLOD; `maxdc` = concentration at ULOD.

- RDL:

  `minrdl` = concentration at LLOD on the lower CI bound curve (`d`
  replaced by `ci_lower["d"]`); `maxrdl` = concentration at ULOD on the
  upper CI bound curve (`d` replaced by `ci_upper["d"]`).

## Usage

``` r
compute_mdc_rdl(param_ci_df, lods, verbose = TRUE)
```

## Arguments

- param_ci_df:

  A long-format CI data frame (same structure as required by
  [`generate_lods`](https://immunoplex.github.io/curveRmetrics/reference/generate_lods.md)).

- lods:

  A data frame as returned by
  [`generate_lods`](https://immunoplex.github.io/curveRmetrics/reference/generate_lods.md)
  with columns `curve_id`, `method`, `llod`, and `ulod`.

- verbose:

  Logical; if `TRUE` (default) emits a message per curve with the four
  computed concentration limits.

## Value

A data frame with columns `curve_id`, `method`, `mindc`, `maxdc`,
`minrdl`, and `maxrdl`.

## Details

All returned concentrations are on the natural scale. Response values
outside \\\[\min(a, d),\\ \max(a, d)\]\\ are silently coerced to `NA`
before inversion.

Only rows flagged as `is_best_model` are processed.
