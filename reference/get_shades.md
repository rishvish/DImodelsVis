# Returns shades of colours

Returns shades of colours

## Usage

``` r
get_shades(colours = c("#808080"), shades = 3)
```

## Arguments

- colours:

  A character vector of colours recognizable by R, to produces shades of

- shades:

  A numeric vector giving the number of shades for each colour

## Value

A list consisting of hex codes describing the shades of each colour

## Examples

``` r
## Shades for a single colour
get_shades(c("red"))
#> $red
#> [1] "#A32929" "#EA3E3E" "#EB6161"
#> 

## Shades for multiple colours
get_shades(c("red", "blue" ,"#A5F8E3", "#808080"), shades = c(2, 3, 4, 5))
#> $red
#> [1] "#BB3030" "#EA5050"
#> 
#> $blue
#> [1] "#222286" "#3535C2" "#5959D5"
#> 
#> $`#A5F8E3`
#> [1] "#72E8C0" "#8DEED1" "#B9F4E4" "#EBFCF8"
#> 
#> $`#808080`
#> [1] "#4D4D4D" "#666666" "#808080" "#9A9A9A" "#B3B3B3"
#> 

## A single value for shade would imply all colours get the same number of shades
get_shades(c("red", "blue" ,"#A5F8E3", "#808080"), shades = 2)
#> $red
#> [1] "#BB3030" "#EA5050"
#> 
#> $blue
#> [1] "#28289A" "#4747CD"
#> 
#> $`#A5F8E3`
#> [1] "#99EFD7" "#DCFAF2"
#> 
#> $`#808080`
#> [1] "#666666" "#9A9A9A"
#> 
```
