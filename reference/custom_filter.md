# Special custom filtering for compositional data

A handy wrapper around the dplyr
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html) function
enabling the user to filter rows which satisfy specific conditions for
compositional data like all equi-proportional communities, or
communities with a given value of richness without having to make any
changes to the data or adding any additional columns. All other
functionalities are same as the dplyr
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html)
function.

## Usage

``` r
custom_filter(data, ..., prop = NULL, special = NULL)
```

## Arguments

- data:

  A data frame containing the compositional variables which should be
  used to perform the filtering.

- ...:

  Any additional arguments specified to the dplyr
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  function. Filtering conditions for any additional variables can also
  be specified here.

- prop:

  A character/numeric vector indicating the columns containing the
  compositional variables in \`data\`.

- special:

  A character string specifying the filtering condition. Four special
  keywords can be specified here for filtering 1. richness: A positive
  integer value to filter communities with a specific number of
  compositional variables (variables with non-zero values). 2. evenness:
  A numeric value between 0 and 1, to filter rows based on the relative
  abundances of the compositional variables where a higher value
  signifies a more even community with equal proportions of all
  variables. 3. equi: A boolean variable indicating whether to filter
  rows containing equi-proportional communities, i.e., communities where
  all variables have the same non-zero proportion. 4. monos: A boolean
  value indicating whether to filter communities containing a single
  compositional variable, i.e., richness == 1. These keywords can be
  combined using any logical operators and can even be combined with any
  other variables in the data. Please use the exact keywords
  (case-sensitive) in the query to get appropriate results. See examples
  for more details.

## Value

A subset of the original data which matches the specified filtering
conditions.

## Examples

``` r
library(DImodels)
library(dplyr)

## Load data
data(sim3)

# The special filter keywords should be specified as a string
# Filter communities containing 3 species
head(custom_filter(data = sim3, prop = 4:12,
                   special = "richness == 3"))
#>   community richness treatment p1        p2 p3        p4 p5 p6 p7        p8
#> 1        46        3         A  0 0.0000000  0 0.3333333  0  0  0 0.3333333
#> 2        46        3         B  0 0.0000000  0 0.3333333  0  0  0 0.3333333
#> 3        46        3         A  0 0.0000000  0 0.3333333  0  0  0 0.3333333
#> 4        46        3         B  0 0.0000000  0 0.3333333  0  0  0 0.3333333
#> 5        47        3         A  0 0.3333333  0 0.0000000  0  0  0 0.3333333
#> 6        47        3         B  0 0.3333333  0 0.0000000  0  0  0 0.3333333
#>          p9 response
#> 1 0.3333333   16.812
#> 2 0.3333333    9.902
#> 3 0.3333333   14.462
#> 4 0.3333333    8.845
#> 5 0.3333333   14.575
#> 6 0.3333333   10.924

# Filter communities at richness 6 OR evenness 0
head(custom_filter(data = sim3, prop = 4:12,
                   special = "richness == 6 | evenness == 0"), 12)
#>    community richness treatment p1 p2 p3 p4 p5 p6 p7 p8 p9 response
#> 1          1        1         A  0  0  0  0  0  0  0  0  1   10.265
#> 2          1        1         B  0  0  0  0  0  0  0  0  1    7.740
#> 3          1        1         A  0  0  0  0  0  0  0  0  1   12.173
#> 4          1        1         B  0  0  0  0  0  0  0  0  1    8.497
#> 5          2        1         A  0  0  0  0  0  0  0  1  0   10.763
#> 6          2        1         B  0  0  0  0  0  0  0  1  0    8.989
#> 7          2        1         A  0  0  0  0  0  0  0  1  0   10.161
#> 8          2        1         B  0  0  0  0  0  0  0  1  0    7.193
#> 9          3        1         A  0  0  0  0  0  0  1  0  0   10.171
#> 10         3        1         B  0  0  0  0  0  0  1  0  0    6.053
#> 11         3        1         A  0  0  0  0  0  0  1  0  0    8.383
#> 12         3        1         B  0  0  0  0  0  0  1  0  0    4.569

# Filter all monoculture AND treatment "A" (treatment is column present in data)
head(custom_filter(data = sim3, prop = 4:12,
                   special = "monos == TRUE & treatment == 'A'"), 10)
#>    community richness treatment p1 p2 p3 p4 p5 p6 p7 p8 p9 response
#> 1          1        1         A  0  0  0  0  0  0  0  0  1   10.265
#> 2          1        1         A  0  0  0  0  0  0  0  0  1   12.173
#> 3          2        1         A  0  0  0  0  0  0  0  1  0   10.763
#> 4          2        1         A  0  0  0  0  0  0  0  1  0   10.161
#> 5          3        1         A  0  0  0  0  0  0  1  0  0   10.171
#> 6          3        1         A  0  0  0  0  0  0  1  0  0    8.383
#> 7          4        1         A  0  0  0  0  0  1  0  0  0   10.182
#> 8          4        1         A  0  0  0  0  0  1  0  0  0    8.240
#> 9          5        1         A  0  0  0  0  1  0  0  0  0   15.467
#> 10         5        1         A  0  0  0  0  1  0  0  0  0   14.790

# Filter all equi proportional communities but NOT monocultures
head(custom_filter(data = sim3, prop = 4:12,
                   special = "equi == TRUE & monos == FALSE"))
#>   community richness treatment p1 p2 p3 p4 p5 p6  p7  p8  p9 response
#> 1        10        2         A  0  0  0  0  0  0 0.0 0.5 0.5   10.734
#> 2        10        2         B  0  0  0  0  0  0 0.0 0.5 0.5    7.465
#> 3        10        2         A  0  0  0  0  0  0 0.0 0.5 0.5   10.054
#> 4        10        2         B  0  0  0  0  0  0 0.0 0.5 0.5    8.615
#> 5        11        2         A  0  0  0  0  0  0 0.5 0.0 0.5    8.614
#> 6        11        2         B  0  0  0  0  0  0 0.5 0.0 0.5    7.851

# Can also use normal filter
sim3 %>% custom_filter(p1 == 1)
#>   community richness treatment p1 p2 p3 p4 p5 p6 p7 p8 p9 response
#> 1         9        1         A  1  0  0  0  0  0  0  0  0   13.539
#> 2         9        1         B  1  0  0  0  0  0  0  0  0   11.768
#> 3         9        1         A  1  0  0  0  0  0  0  0  0   11.869
#> 4         9        1         B  1  0  0  0  0  0  0  0  0   10.393

# Both special filtering and normal filtering can be combined as well
sim3 %>% custom_filter(prop = paste0("p", 1:9),
                       special = "richness == 1",
                       community %in% c(7, 9))
#>   community richness treatment p1 p2 p3 p4 p5 p6 p7 p8 p9 response
#> 1         7        1         A  0  0  1  0  0  0  0  0  0   12.413
#> 2         7        1         B  0  0  1  0  0  0  0  0  0    8.376
#> 3         7        1         A  0  0  1  0  0  0  0  0  0   12.264
#> 4         7        1         B  0  0  1  0  0  0  0  0  0    8.371
#> 5         9        1         A  1  0  0  0  0  0  0  0  0   13.539
#> 6         9        1         B  1  0  0  0  0  0  0  0  0   11.768
#> 7         9        1         A  1  0  0  0  0  0  0  0  0   11.869
#> 8         9        1         B  1  0  0  0  0  0  0  0  0   10.393
```
