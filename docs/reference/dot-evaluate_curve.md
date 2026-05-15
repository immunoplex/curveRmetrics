# Evaluate a dose-response curve at a given concentration

All parameters and `x` are expected on the **natural scale** (`c` =
natural EC50, `x` = natural concentration). The \\\log(x/c)/b\\
substitution is used throughout so that the formulas accept raw
concentration directly.

## Usage

``` r
.evaluate_curve(p, x)
```

## Arguments

- p:

  A named list of curve parameters with elements `family`, `a`, `b`,
  `c`, `d`, and `g`.

- x:

  Numeric; concentration value(s) at which to evaluate the curve.

## Value

Numeric; predicted response value(s).

## Details

Supported families:

- logistic4:

  \\y = d + (a-d) / (1 + \exp((x-c)/b))\\

- logistic5:

  \\y = d + (a-d) / (1 + \exp((x-c)/b))^g\\

- loglogistic4:

  \\y = a + (d-a) / (1 + (x/c)^b)\\

- loglogistic5:

  \\y = a + (d-a)(1 + g\exp(-b(x-c)))^{-1/g}\\

- gompertz4:

  \\y = a + (d-a)\exp(-\exp(-b(x-c)))\\
