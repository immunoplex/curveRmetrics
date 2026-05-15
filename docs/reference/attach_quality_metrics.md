# Attach quality metrics to a curves data frame

Left-joins inflection points, assay sensitivity, limits of detection,
reliable detection limits, and limits of quantification onto the curves
data frame. The `model_name` column is first resolved to its canonical
family name via
[`.canonical_family`](https://immunoplex.github.io/curveRmetrics/reference/dot-canonical_family.md).

## Usage

``` r
attach_quality_metrics(curves_df, lods, rdls, sensitivity, loqs, inflection)
```

## Arguments

- curves_df:

  A data frame with one row per `(curve_id, method)` combination,
  containing at minimum a `model_name` column.

- lods:

  A data frame as returned by
  [`generate_lods`](https://immunoplex.github.io/curveRmetrics/reference/generate_lods.md)
  with columns `curve_id`, `method`, `llod`, and `ulod`.

- rdls:

  A data frame of reliable detection limits with columns `curve_id`,
  `method`, `minrdl`, and `maxrdl`.

- sensitivity:

  A data frame as returned by
  [`compute_assay_sensitivity`](https://immunoplex.github.io/curveRmetrics/reference/compute_assay_sensitivity.md)
  with columns `curve_id`, `method`, and `dydx_inflect` (the `inflect_x`
  column, if present, is dropped before merging to avoid duplication).

- loqs:

  A data frame as returned by
  [`compute_loqs`](https://immunoplex.github.io/curveRmetrics/reference/compute_loqs.md)
  with columns `curve_id`, `method`, `lloq`, `uloq`, `lloq_y`, and
  `uloq_y`.

- inflection:

  A data frame as returned by
  [`compute_inflection_point`](https://immunoplex.github.io/curveRmetrics/reference/compute_inflection_point.md)
  with columns `curve_id`, `method`, `inflect_x`, and `inflect_y`.

## Value

`curves_df` with additional columns from each quality-metric data frame,
joined on `(curve_id, method)`.
