# Get all equi-proportional communities at specific levels of richness

Get all equi-proportional communities at specific levels of richness

## Usage

``` r
get_equi_comms(
  nvars,
  richness_lvl = 1:nvars,
  variables = paste0("Var", 1:nvars),
  threshold = 1e+06
)
```

## Arguments

- nvars:

  Number of variables in the design

- richness_lvl:

  The richness levels (number of non-zero compositional variables in a
  community) at which to return the equi-proportional communities.
  Defaults to each richness level from 1 up to \`nvars\` (both
  inclusive).

- variables:

  Names for the variables. Will be used as column names for the final
  result. Default is "Var" followed by column number.

- threshold:

  The maximum number of communities to select for each level of richness
  for situations when there are too many equi-proportional communities.
  Default value is a million.  
  Note: if threshold \< \`number of possible equi-proportional
  communities\` at a given level of richness, a random selection of
  communities equal to the number specified in threshold would be
  returned.

## Value

A dataframe consisting all or a random selection of equi-proportional
communities at each level of richness

## Examples

``` r
## Get all equi-proportional communities for each level of richness upto 10
data10 <- get_equi_comms(10)
head(data10, 12)
#>    Var1 Var2 Var3 Var4 Var5 Var6 Var7 Var8 Var9 Var10 Richness
#> 1   1.0  0.0  0.0    0    0    0    0    0    0     0        1
#> 2   0.0  1.0  0.0    0    0    0    0    0    0     0        1
#> 3   0.0  0.0  1.0    0    0    0    0    0    0     0        1
#> 4   0.0  0.0  0.0    1    0    0    0    0    0     0        1
#> 5   0.0  0.0  0.0    0    1    0    0    0    0     0        1
#> 6   0.0  0.0  0.0    0    0    1    0    0    0     0        1
#> 7   0.0  0.0  0.0    0    0    0    1    0    0     0        1
#> 8   0.0  0.0  0.0    0    0    0    0    1    0     0        1
#> 9   0.0  0.0  0.0    0    0    0    0    0    1     0        1
#> 10  0.0  0.0  0.0    0    0    0    0    0    0     1        1
#> 11  0.5  0.5  0.0    0    0    0    0    0    0     0        2
#> 12  0.5  0.0  0.5    0    0    0    0    0    0     0        2

## Change variable names
data4 <- get_equi_comms(4, variables = c("Lollium perenne", "Chichorum intybus",
                                         "Trifolium repens", "Trifolium pratense"))
head(data4)
#>   Lollium perenne Chichorum intybus Trifolium repens Trifolium pratense
#> 1             1.0               0.0              0.0                  0
#> 2             0.0               1.0              0.0                  0
#> 3             0.0               0.0              1.0                  0
#> 4             0.0               0.0              0.0                  1
#> 5             0.5               0.5              0.0                  0
#> 6             0.5               0.0              0.5                  0
#>   Richness
#> 1        1
#> 2        1
#> 3        1
#> 4        1
#> 5        2
#> 6        2

## Get equi-proportional communities at specific levels of richness
## Get all equi-proportional communities of four variables at richness
## levels 1 and 3
data4_13 <- get_equi_comms(nvars = 4, richness = c(1, 3))
data4_13
#>        Var1      Var2      Var3      Var4 Richness
#> 1 1.0000000 0.0000000 0.0000000 0.0000000        1
#> 2 0.0000000 1.0000000 0.0000000 0.0000000        1
#> 3 0.0000000 0.0000000 1.0000000 0.0000000        1
#> 4 0.0000000 0.0000000 0.0000000 1.0000000        1
#> 5 0.3333333 0.3333333 0.3333333 0.0000000        3
#> 6 0.3333333 0.3333333 0.0000000 0.3333333        3
#> 7 0.3333333 0.0000000 0.3333333 0.3333333        3
#> 8 0.0000000 0.3333333 0.3333333 0.3333333        3

## If threshold is specified and it is less than the number of possible
## equi-proportional communites at a given level of richness, then a
## random selection of communities from the total possible would be returned
## Return only 2 random equi-proportional communities at the chosen richness
## levels
data4_13_2 <- get_equi_comms(nvars = 4, richness = c(1, 3), threshold = 2)
data4_13_2
#>        Var1      Var2      Var3      Var4 Richness
#> 1 0.0000000 1.0000000 0.0000000 0.0000000        1
#> 2 0.0000000 1.0000000 0.0000000 0.0000000        1
#> 3 0.3333333 0.3333333 0.3333333 0.0000000        3
#> 4 0.3333333 0.0000000 0.3333333 0.3333333        3

## Set threshold to a very high positive number to ensure
## random selection is never performed
data_no_random <- get_equi_comms(nvars = 15,
                                 threshold = .Machine$integer.max)
head(data_no_random)
#>   Var1 Var2 Var3 Var4 Var5 Var6 Var7 Var8 Var9 Var10 Var11 Var12 Var13 Var14
#> 1    1    0    0    0    0    0    0    0    0     0     0     0     0     0
#> 2    0    1    0    0    0    0    0    0    0     0     0     0     0     0
#> 3    0    0    1    0    0    0    0    0    0     0     0     0     0     0
#> 4    0    0    0    1    0    0    0    0    0     0     0     0     0     0
#> 5    0    0    0    0    1    0    0    0    0     0     0     0     0     0
#> 6    0    0    0    0    0    1    0    0    0     0     0     0     0     0
#>   Var15 Richness
#> 1     0        1
#> 2     0        1
#> 3     0        1
#> 4     0        1
#> 5     0        1
#> 6     0        1
```
