# Resolve a raw model/family name to its canonical form

Maps NLS model names (from `model_functions.R`) and Bayesian family
codes (from `curve_families.R`) to a single canonical string used
throughout this module.

## Usage

``` r
.canonical_family(raw_name)
```

## Arguments

- raw_name:

  A character string; the model or family name to resolve.

## Value

A single canonical family name (character).
