# Default theme for DImodelsVis

Default theme for DImodelsVis

## Usage

``` r
theme_DI(
  font_size = 14,
  font_family = "",
  legend = c("top", "bottom", "left", "right", "none")
)
```

## Arguments

- font_size:

  Base font size for text across the plot

- font_family:

  Font family for text across the plot

- legend:

  One of c("top", "bottom", "left", "right", "none") specifying the
  position of the legend. The legend position can also be specified as a
  numeric vector of form c(x, y) with x and y having values between 0
  and 1. If specified as a numeric vector the legend within the plotting
  region where c(0,0) corresponds to the "bottom left" and c(1,1)
  corresponds to the "top right" position. The default position is
  "top".

## Value

A ggplot theme object

## Examples

``` r
library(ggplot2)

plot_data <- mtcars
plot_data$gear <- as.factor(plot_data$gear)
ggplot(data = plot_data,
       aes(x = mpg, y = disp, colour = gear))+
   geom_point(size = 3)+
   facet_wrap(~cyl) +
   theme_DI()
```
