# Locate u-intervals that bracket local extrema of \\d^2y/du^2\\

Scans a discrete second-derivative series for sign changes in
consecutive first differences, which indicate the presence of a local
maximum or minimum. For each detected sign change, a bracketing interval
of neighbouring \\u = \log\_{10}(x)\\ values is returned.

## Usage

``` r
.sderiv_brackets(u, d2y)
```

## Arguments

- u:

  Numeric vector of \\\log\_{10}(x)\\ grid positions.

- d2y:

  Numeric vector of second-derivative values at each position in `u`;
  `NA` values are removed before scanning.

## Value

A list of two-element numeric vectors, each giving the \\\[u\_{lo},\\
u\_{hi}\]\\ bounds of one bracket. An empty list is returned when fewer
than four non-`NA` points are available or no sign change is detected.
