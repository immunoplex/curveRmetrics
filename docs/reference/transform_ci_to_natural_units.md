# Transform long-format parameter CI data to the natural scale

Adds natural-scale counterparts of the estimate and confidence bound
columns to a long-format parameter CI data frame, suitable for use with
both `freq_param_ci_df` and `bayes_param_ci_df` (or their row-bound
combination).

## Usage

``` r
transform_ci_to_natural_units(param_ci_df)
```

## Arguments

- param_ci_df:

  A long-format data frame with at minimum the columns `parameter`,
  `method`, `estimate`, `conf_lower`, and `conf_upper`.

## Value

`param_ci_df` with three additional columns:

- estimate_nat:

  Point estimate on the natural scale.

- conf_lower_nat:

  Lower confidence bound on the natural scale (`min` of the transformed
  bounds).

- conf_upper_nat:

  Upper confidence bound on the natural scale (`max` of the transformed
  bounds).

## Details

Each value is converted via
[`.to_nat`](https://immunoplex.github.io/curveRmetrics/reference/dot-to_nat.md)
according to its `parameter` and `method` fields. Because \\10^x\\ and
\\x \cdot \ln(10)\\ are both monotone increasing, CI ordering is
preserved for `a`, `c`, `d`, and `b` under the frequentist method. For
the Bayesian `b` conversion (\\1/x\\), which is monotone *decreasing*,
the lower and upper bounds are swapped after transformation so that
`conf_lower_nat` always holds the smaller value.
