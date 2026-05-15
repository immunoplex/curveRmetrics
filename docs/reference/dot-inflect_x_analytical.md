# Compute the analytical inflection point as the response-range midpoint

Finds the natural-scale concentration \\x^\*\\ at which the curve equals
the geometric midpoint of its response range, i.e. \\y^\* = \sqrt{a
\cdot d}\\ (Bayesian) or the equivalent arithmetic midpoint in
\\\log\_{10}\\ space (frequentist). This definition is unified across
both fitting methods and requires no root-finding.

## Usage

``` r
.inflect_x_analytical(p, method)
```

## Arguments

- p:

  A named list of curve parameters. For frequentist curves the
  parameters are on the raw \\\log\_{10}\\ scale; for Bayesian curves
  they are on the natural scale (with `c = ln(c_raw)`).

- method:

  Character; `"frequentist"` or `"bayesian"`.

## Value

Numeric scalar; natural-scale concentration at the midpoint inflection.
Falls back to the back-transformed EC50 if the analytical condition
cannot be evaluated.
