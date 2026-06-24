# Project 3-d compositional data onto x-y plane and vice versa

Points in the 3-d simplex space with coordinates (x, y ,z) such that x +
y + z = 1 are projected into the 2-d plane they reside in. This function
can be used to convert the 3-d compositional data into 2-d and then be
overlayed on the plots output by
[`ternary_plot`](https://rishvish.github.io/DImodelsVis/reference/ternary_plot.md),
[`conditional_ternary_plot`](https://rishvish.github.io/DImodelsVis/reference/conditional_ternary_plot.md)
and
[`grouped_ternary_plot`](https://rishvish.github.io/DImodelsVis/reference/grouped_ternary_plot.md).

## Usage

``` r
prop_to_tern_proj(data, prop, x = ".x", y = ".y")

tern_to_prop_proj(data, x, y, prop = c("p1", "p2", "p3"))
```

## Arguments

- data:

  A data-frame containing the x-y coordinates of the points.

- prop:

  A character vector specifying the columns names of variable containing
  the projected compositions. Default is "p1", "p2", and "p3".

- x:

  A character string specifying the name for the column containing the x
  component of the x-y projection of the simplex.

- y:

  A character string specifying the name for the column containing the y
  component of the x-y projection of the simplex.

## Value

A data-frame with the following two columns appended (when transforming
to x-y projection)

- .x (or value specified in "x"):

  The x component of the x-y projection of the simplex point.

- .y (or value specified in "y"):

  The y component of the x-y projection of the simplex point.

A data-frame with the following three columns appended (when
transforming to compositional projection)

- p1 (or first value specified in "prop"):

  The first component of the 3-d simplex point.

- p2 (or second value specified in "prop"):

  The second component of the 3-d simplex point.

- p3 (or third value specified in "prop"):

  The third component of the 3-d simplex point.

## Examples

``` r
## Convert proportions to x-y co-ordinates
library(DImodels)
data(sim0)
sim0 <- sim0[1:16, ]

prop_to_tern_proj(data = sim0, prop = c("p1", "p2", "p3"))
#>           p1        p2        p3   .x        .y community richness response
#> 1  1.0000000 0.0000000 0.0000000 0.50 0.8660254         1        1   24.855
#> 2  0.0000000 1.0000000 0.0000000 0.00 0.0000000         2        1   19.049
#> 3  0.0000000 0.0000000 1.0000000 1.00 0.0000000         3        1   16.292
#> 4  0.8000000 0.2000000 0.0000000 0.40 0.6928203         4        2   31.529
#> 5  0.2000000 0.8000000 0.0000000 0.10 0.1732051         5        2   25.102
#> 6  0.8000000 0.0000000 0.2000000 0.60 0.6928203         6        2   24.615
#> 7  0.2000000 0.0000000 0.8000000 0.90 0.1732051         7        2   18.654
#> 8  0.0000000 0.8000000 0.2000000 0.20 0.0000000         8        2   24.697
#> 9  0.0000000 0.2000000 0.8000000 0.80 0.0000000         9        2   25.017
#> 10 0.5000000 0.5000000 0.0000000 0.25 0.4330127        10        2   32.743
#> 11 0.5000000 0.0000000 0.5000000 0.75 0.4330127        11        2   25.320
#> 12 0.0000000 0.5000000 0.5000000 0.50 0.0000000        12        2   30.214
#> 13 0.6000000 0.2000000 0.2000000 0.50 0.5196152        13        3   30.511
#> 14 0.2000000 0.6000000 0.2000000 0.30 0.1732051        14        3   31.005
#> 15 0.2000000 0.2000000 0.6000000 0.70 0.1732051        15        3   28.964
#> 16 0.3333333 0.3333333 0.3333333 0.50 0.2886751        16        3   32.091

# Change names of the x and y projections
prop_to_tern_proj(data = sim0, prop = c("p1", "p2", "p3"),
                  x = "x-proj", y = "y-proj")
#>           p1        p2        p3 x-proj    y-proj community richness response
#> 1  1.0000000 0.0000000 0.0000000   0.50 0.8660254         1        1   24.855
#> 2  0.0000000 1.0000000 0.0000000   0.00 0.0000000         2        1   19.049
#> 3  0.0000000 0.0000000 1.0000000   1.00 0.0000000         3        1   16.292
#> 4  0.8000000 0.2000000 0.0000000   0.40 0.6928203         4        2   31.529
#> 5  0.2000000 0.8000000 0.0000000   0.10 0.1732051         5        2   25.102
#> 6  0.8000000 0.0000000 0.2000000   0.60 0.6928203         6        2   24.615
#> 7  0.2000000 0.0000000 0.8000000   0.90 0.1732051         7        2   18.654
#> 8  0.0000000 0.8000000 0.2000000   0.20 0.0000000         8        2   24.697
#> 9  0.0000000 0.2000000 0.8000000   0.80 0.0000000         9        2   25.017
#> 10 0.5000000 0.5000000 0.0000000   0.25 0.4330127        10        2   32.743
#> 11 0.5000000 0.0000000 0.5000000   0.75 0.4330127        11        2   25.320
#> 12 0.0000000 0.5000000 0.5000000   0.50 0.0000000        12        2   30.214
#> 13 0.6000000 0.2000000 0.2000000   0.50 0.5196152        13        3   30.511
#> 14 0.2000000 0.6000000 0.2000000   0.30 0.1732051        14        3   31.005
#> 15 0.2000000 0.2000000 0.6000000   0.70 0.1732051        15        3   28.964
#> 16 0.3333333 0.3333333 0.3333333   0.50 0.2886751        16        3   32.091
## Convert x-y co-ordinates to proportions
library(DImodels)
data(sim0)
sim0 <- sim0[1:16, ]

proj_data <- prop_to_tern_proj(data = sim0, prop = c("p1", "p2", "p3"))

tern_to_prop_proj(data = proj_data, x = ".x", y = ".y")
#>      .x        .y        p1        p2        p3 community richness response
#> 1  0.50 0.8660254 1.0000000 0.0000000 0.0000000         1        1   24.855
#> 2  0.00 0.0000000 0.0000000 1.0000000 0.0000000         2        1   19.049
#> 3  1.00 0.0000000 0.0000000 0.0000000 1.0000000         3        1   16.292
#> 4  0.40 0.6928203 0.8000000 0.2000000 0.0000000         4        2   31.529
#> 5  0.10 0.1732051 0.2000000 0.8000000 0.0000000         5        2   25.102
#> 6  0.60 0.6928203 0.8000000 0.0000000 0.2000000         6        2   24.615
#> 7  0.90 0.1732051 0.2000000 0.0000000 0.8000000         7        2   18.654
#> 8  0.20 0.0000000 0.0000000 0.8000000 0.2000000         8        2   24.697
#> 9  0.80 0.0000000 0.0000000 0.2000000 0.8000000         9        2   25.017
#> 10 0.25 0.4330127 0.5000000 0.5000000 0.0000000        10        2   32.743
#> 11 0.75 0.4330127 0.5000000 0.0000000 0.5000000        11        2   25.320
#> 12 0.50 0.0000000 0.0000000 0.5000000 0.5000000        12        2   30.214
#> 13 0.50 0.5196152 0.6000000 0.2000000 0.2000000        13        3   30.511
#> 14 0.30 0.1732051 0.2000000 0.6000000 0.2000000        14        3   31.005
#> 15 0.70 0.1732051 0.2000000 0.2000000 0.6000000        15        3   28.964
#> 16 0.50 0.2886751 0.3333333 0.3333333 0.3333333        16        3   32.091

# Change prop names
tern_to_prop_proj(data = proj_data, x = ".x", y = ".y",
                  prop = c("prop1", "prop2", "prop3"))
#>      .x        .y     prop1     prop2     prop3        p1        p2        p3
#> 1  0.50 0.8660254 1.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000
#> 2  0.00 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000 1.0000000 0.0000000
#> 3  1.00 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000 1.0000000
#> 4  0.40 0.6928203 0.8000000 0.2000000 0.0000000 0.8000000 0.2000000 0.0000000
#> 5  0.10 0.1732051 0.2000000 0.8000000 0.0000000 0.2000000 0.8000000 0.0000000
#> 6  0.60 0.6928203 0.8000000 0.0000000 0.2000000 0.8000000 0.0000000 0.2000000
#> 7  0.90 0.1732051 0.2000000 0.0000000 0.8000000 0.2000000 0.0000000 0.8000000
#> 8  0.20 0.0000000 0.0000000 0.8000000 0.2000000 0.0000000 0.8000000 0.2000000
#> 9  0.80 0.0000000 0.0000000 0.2000000 0.8000000 0.0000000 0.2000000 0.8000000
#> 10 0.25 0.4330127 0.5000000 0.5000000 0.0000000 0.5000000 0.5000000 0.0000000
#> 11 0.75 0.4330127 0.5000000 0.0000000 0.5000000 0.5000000 0.0000000 0.5000000
#> 12 0.50 0.0000000 0.0000000 0.5000000 0.5000000 0.0000000 0.5000000 0.5000000
#> 13 0.50 0.5196152 0.6000000 0.2000000 0.2000000 0.6000000 0.2000000 0.2000000
#> 14 0.30 0.1732051 0.2000000 0.6000000 0.2000000 0.2000000 0.6000000 0.2000000
#> 15 0.70 0.1732051 0.2000000 0.2000000 0.6000000 0.2000000 0.2000000 0.6000000
#> 16 0.50 0.2886751 0.3333333 0.3333333 0.3333333 0.3333333 0.3333333 0.3333333
#>    community richness response
#> 1          1        1   24.855
#> 2          2        1   19.049
#> 3          3        1   16.292
#> 4          4        2   31.529
#> 5          5        2   25.102
#> 6          6        2   24.615
#> 7          7        2   18.654
#> 8          8        2   24.697
#> 9          9        2   25.017
#> 10        10        2   32.743
#> 11        11        2   25.320
#> 12        12        2   30.214
#> 13        13        3   30.511
#> 14        14        3   31.005
#> 15        15        3   28.964
#> 16        16        3   32.091
```
