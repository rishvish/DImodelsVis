# Return colour-blind friendly colours

Utility function to return either a distinct colour-blind friendly
colour for each variable or if a functional grouping is specified, then
shades of the same colour for variables within a functional group

## Usage

``` r
get_colours(vars, FG = NULL)
```

## Arguments

- vars:

  Either a numeric value \`n\` to get n colours, or a character vector
  of values where each value will be mapped to a colour.

- FG:

  A character vector describing the functional grouping to which each
  variable belongs. Variables within the same group will have different
  shades of the same colour.

## Value

A named vector containing the hex codes of colours

## Examples

``` r
## Get n colours
get_colours(vars = 4)
#> [1] "#009E73" "#AA4499" "#0072B2" "#F0E442"

# Get a color-map for each value specified in vars
get_colours(vars = c("p1", "p2", "p3", "p4"))
#> [1] "#009E73" "#AA4499" "#0072B2" "#F0E442"

# Group values of vars using FG. Variables in the same group
# will have same shades of a colour
get_colours(vars = 4, FG = c("G1", "G1", "G2", "G2"))
#> [1] "#2D6852" "#5FCBA0" "#7A406F" "#B26CA5"
```
