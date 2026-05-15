# Compute the working-axis bounds for the second derivative grid

Returns the \\\log\_{10}(x)\\ interval that spans the region of interest
for the second derivative, centred on \\u_c = \log\_{10}(c)\\ and
widened by a family-specific multiple of the slope parameter:

- logistic4/5, gompertz4:

  span = \\6 \cdot b\_{nat} / \ln(10)\\, scaled by \\g^{0.3}\\ for the
  5-parameter variant.

- loglogistic4/5:

  span = \\6 / b\\.

- gompertz4:

  left tail extended by \\2 \times\\ span to capture the asymmetric
  lower extremum.

## Usage

``` r
.sderiv_scale(p)
```

## Arguments

- p:

  A named list of natural-scale curve parameters with elements `family`,
  `b` (natural-log-scale slope for logistic/gompertz; log10-scale slope
  for loglogistic), `c` (natural EC50), and optionally `g`.

## Value

A named list with elements `u_lo`, `u_hi` (the \\\log\_{10}(x)\\
bounds), and `log_scale = TRUE`.
