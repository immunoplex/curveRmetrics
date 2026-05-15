# Plot the second derivative of a fitted curve for both methods

Renders a log-scale line plot of \\d^2y/dx^2\\ against concentration for
a single curve, overlaying the frequentist and Bayesian traces. The
x-axis is automatically zoomed to the region where the second derivative
is meaningfully non-zero (greater than 1 % of its maximum absolute
value).

## Usage

``` r
compare_second_derivative(second_derivative_df, curve_id)
```

## Arguments

- second_derivative_df:

  A data frame with at minimum the columns `curve_id`, `method`,
  `concentration`, and `d2x_y` (the second derivative of the response
  with respect to concentration).

- curve_id:

  Scalar character or integer; the identifier of the curve to plot.

## Value

A `ggplot` object.
