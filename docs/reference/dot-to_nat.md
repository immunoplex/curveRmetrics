# Convert a single parameter value to the natural scale

Low-level scalar/vector converter applied element-wise by
[`transform_ci_to_natural_units`](https://immunoplex.github.io/curveRmetrics/reference/transform_ci_to_natural_units.md).
Conversion rules differ by parameter and fitting method:

## Usage

``` r
.to_nat(x, param, method)
```

## Arguments

- x:

  Numeric value(s) to convert.

- param:

  Character; parameter name — one of `"a"`, `"b"`, `"c"`, `"d"`, or
  `"g"`.

- method:

  Character; fitting method — `"frequentist"` or `"bayesian"`.

## Value

Numeric value(s) on the natural scale.

## Details

|               |                    |              |
|---------------|--------------------|--------------|
| **Parameter** | **Frequentist**    | **Bayesian** |
| `a`, `d`      | \\10^x\\           | unchanged    |
| `c`           | \\10^x\\           | \\\ln(x)\\   |
| `b`           | \\x \cdot \ln 10\\ | \\1/x\\      |
| `g`           | unchanged          | unchanged    |

Any unrecognised parameter name is passed through unchanged.
