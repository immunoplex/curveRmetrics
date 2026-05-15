# Transform curve parameters to natural (untransformed) units

Converts raw fitted parameters stored in a data frame to their
natural-scale equivalents, appending `*_nat` columns. The transformation
rules are:

## Usage

``` r
transform_to_natural_units(df)
```

## Arguments

- df:

  A data frame with at minimum the columns `curve_id`, `method`,
  `model_name`, `a`, `b`, `c`, `d`, and `g`.

## Value

`df` with additional columns `a_nat`, `b_nat`, `c_nat`, `d_nat`, and
`g_nat` on the natural-unit scale.

## Details

- Frequentist:

  Model was fit on \\\log\_{10}(x)\\, so `c` (EC50) is back-transformed
  via \\10^c\\, and `b` (Hill slope) is multiplied by \\\ln(10)\\ to
  convert from the \\\log\_{10}\\ scale to the natural-log scale.
  Asymptotes `a` and `d` are also back-transformed via \\10^a\\,
  \\10^d\\.

- Bayesian:

  Model is parameterised in \\\ln(x)\\, with `c` already representing
  the median effective concentration; `a` and `d` are already on the
  response scale. `c_nat` is set to \\\ln(c)\\ and `b_nat` to \\1/b\\ to
  match the internal Stan parameterisation.

The asymmetry parameter `g` is dimensionless and unaffected by any scale
transformation; it is set to 1 for symmetric (4-parameter) families.
