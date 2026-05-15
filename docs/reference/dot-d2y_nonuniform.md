# Non-uniform central-difference second derivative in log10-x space

Computes \\d^2y/du^2\\ where \\u = \log\_{10}(x)\\ using a non-uniform
central-difference formula that is stable on log-spaced grids. Working
in \\\log\_{10}(x)\\ space avoids the \\1/x^2\\ blow-up that would occur
when differentiating in natural concentration.

## Usage

``` r
.d2y_nonuniform(x, y)
```

## Arguments

- x:

  Numeric vector of natural-scale concentration values (strictly
  positive, log-spaced).

- y:

  Numeric vector of response values corresponding to `x`.

## Value

Numeric vector of the same length as `x` containing \\d^2y/du^2\\; the
first and last elements are always `NA`, as are any interior points
where `y` contains `NA` or where adjacent grid spacings are below
machine epsilon.

## Details

The formula at interior point \\i\\ is: \$\$ \frac{d^2y}{du^2}\bigg\|\_i
\approx \frac{2}{h_l + h_r} \left\[\frac{y\_{i+1} - y_i}{h_r} -
\frac{y_i - y\_{i-1}}{h_l}\right\] \$\$ where \\h_l = u_i - u\_{i-1}\\
and \\h_r = u\_{i+1} - u_i\\.
