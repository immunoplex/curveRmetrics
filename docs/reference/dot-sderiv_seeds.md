# Analytical seed positions for the extrema of \\d^2y/du^2\\

Returns a two-element numeric vector of \\\log\_{10}(x^\*)\\ values at
which \\d^3y/du^3 = 0\\, i.e. the expected locations of the maximum and
minimum of the second derivative. These seeds are used as starting
points for the dense refinement grid in
[`.compute_second_deriv`](https://immunoplex.github.io/curveRmetrics/reference/dot-compute_second_deriv.md).

## Usage

``` r
.sderiv_seeds(p)
```

## Arguments

- p:

  A named list of natural-scale curve parameters (see
  [`.sderiv_scale`](https://immunoplex.github.io/curveRmetrics/reference/dot-sderiv_scale.md)).

## Value

A numeric vector of length 2 containing the two seed positions in
\\\log\_{10}(x)\\ units.

## Details

Derivations use \\t = \ln(x/c)\\, \\u = \log\_{10}(x)\\:

- logistic4:

  \\z^2 - 4z + 1 = 0\\, \\z^\* = 2 \pm \sqrt{3}\\, \\u^\* = u_c + b
  \cdot \log\_{10}(2 \pm \sqrt{3})\\ (exact).

- logistic5:

  logistic4 seeds scaled by \\g^{0.3}\\ (approx).

- loglogistic4:

  \\u^\* = u_c \pm \log\_{10}(2 \pm \sqrt{3}) / b\\ (exact).

- loglogistic5:

  loglogistic4 seeds scaled by \\g^{0.3}\\ (approx).

- gompertz4:

  \\v^\* = (3 \pm \sqrt{5})/2\\, \\u^\* = u_c \pm \mathrm{LN\\PHI2} / (b
  \cdot \ln 10)\\ (exact).
