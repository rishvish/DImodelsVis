# Combine variable proportions into groups

Combine variable proportions into groups

## Usage

``` r
group_prop(data, prop, FG = NULL)
```

## Arguments

- data:

  A data frame containing the compositional variables which need to be
  grouped.

- prop:

  A character/numeric vector indicating the columns containing the
  compositional variables in \`data\`.

- FG:

  A character vector of same length as \`prop\` specifying the group
  each variable belongs to.

## Value

A data-frame with additional columns appended to the end that contain
the grouped variable proportions.

## Examples

``` r
library(DImodels)

data(sim1)

head(group_prop(data = sim1, prop = 3:6,
                FG = c("Gr1", "Gr1", "Gr1", "Gr2")))
#>   community block  p1  p2  p3  p4 response Gr1 Gr2
#> 1         1     1 0.7 0.1 0.1 0.1   10.815 0.9 0.1
#> 2         1     2 0.7 0.1 0.1 0.1   11.232 0.9 0.1
#> 3         1     3 0.7 0.1 0.1 0.1   10.192 0.9 0.1
#> 4         1     4 0.7 0.1 0.1 0.1    8.157 0.9 0.1
#> 5         2     1 0.1 0.7 0.1 0.1    6.724 0.9 0.1
#> 6         2     2 0.1 0.7 0.1 0.1   11.093 0.9 0.1

head(group_prop(data = sim1, prop = 3:6,
                FG = c("Group1", "Group2", "Group1", "Group3")))
#>   community block  p1  p2  p3  p4 response Group1 Group2 Group3
#> 1         1     1 0.7 0.1 0.1 0.1   10.815    0.8    0.1    0.1
#> 2         1     2 0.7 0.1 0.1 0.1   11.232    0.8    0.1    0.1
#> 3         1     3 0.7 0.1 0.1 0.1   10.192    0.8    0.1    0.1
#> 4         1     4 0.7 0.1 0.1 0.1    8.157    0.8    0.1    0.1
#> 5         2     1 0.1 0.7 0.1 0.1    6.724    0.2    0.7    0.1
#> 6         2     2 0.1 0.7 0.1 0.1   11.093    0.2    0.7    0.1

## Data is returned as is, if no groups are specified in FG
head(group_prop(data = sim1, prop = 3:6))
#>   community block  p1  p2  p3  p4 response
#> 1         1     1 0.7 0.1 0.1 0.1   10.815
#> 2         1     2 0.7 0.1 0.1 0.1   11.232
#> 3         1     3 0.7 0.1 0.1 0.1   10.192
#> 4         1     4 0.7 0.1 0.1 0.1    8.157
#> 5         2     1 0.1 0.7 0.1 0.1    6.724
#> 6         2     2 0.1 0.7 0.1 0.1   11.093
```
