# Data preparation for visualising change in response over diversity gradient

Helper function for creating the data to visualise a scatter-plot of the
response over a diversity gradient. The "richness" and "evenness"
diversity gradients are currently supported. The average (predicted)
response is calculated from all communities present at a given level of
the chosen diversity gradient in \`data\`. The output of this function
can be passed to the
[`gradient_change_plot`](https://rishvish.github.io/DImodelsVis/reference/gradient_change_plot.md)
function to visualise results.

## Usage

``` r
gradient_change_data(
  data,
  prop,
  add_var = list(),
  gradient = c("richness", "evenness"),
  prediction = TRUE,
  ...
)
```

## Arguments

- data:

  A data-frame consisting of variable proportions and any other
  necessary variables to make predictions from \`model\` or
  \`coefficients\`.

- prop:

  A vector identifying the column-names or indices of the columns
  containing the variable proportions in \`data\`.

- add_var:

  A list specifying values for additional predictor variables in the
  model independent of the compositional predictor variables. This could
  be useful for comparing the predictions across different values for a
  non-compositional variable. If specified as a list, it will be
  expanded to show a plot for each unique combination of values
  specified, while if specified as a data-frame, one plot would be
  generated for each row in the data and they will be arranged in a grid
  according to the value specified in \`nrow\` and \`ncol\`.

- gradient:

  Diversity gradient to show on the X-axis, one of "richness" or
  "evenness". Defaults to "richness". See \`Details\` for more
  information.

- prediction:

  A logical value indicating whether to pass the final data to the
  \`add_prediction\` function and append the predictions to the data.
  Default value is TRUE, but often it would be desirable to make
  additional changes to the data before making any predictions, so the
  user can set this to FALSE and manually call the \`add_prediction\`
  function.

- ...:

  Arguments passed on to
  [`add_prediction`](https://rishvish.github.io/DImodelsVis/reference/add_prediction.md)

  `model`

  :   A regression model object which will be used to make predictions
      for the observations in \`data\`. Will override \`coefficients\`
      if specified.

  `coefficients`

  :   If a regression model is not available (or can't be fit in R), the
      regression coefficients from a model fit in some other language
      can be used to calculate predictions. However, the user would have
      to ensure there's an appropriate one-to-one positional mapping
      between the data columns and the coefficient values. Further, they
      would also have to provide a variance-covariance matrix of the
      coefficients in the \`vcov\` parameter if they want the associated
      CI for the prediction or it would not be possible to calculate
      confidence/prediction intervals using this method.

  `vcov`

  :   If regression coefficients are specified, then the
      variance-covariance matrix of the coefficients can be specified
      here to calculate the associated confidence interval around each
      prediction. Failure to do so would result in no confidence
      intervals being returned. Ensure \`coefficients\` and \`vcov\`
      have the same positional mapping with the data.

  `coeff_cols`

  :   If \`coefficients\` are specified and a one-to-one positional
      mapping between the data-columns and coefficient vector is not
      present. A character string or numeric index can be specified here
      to reorder the data columns and match the corresponding
      coefficient value to the respective data column. See the "Use
      model coefficients for prediction" section in examples.

  `conf.level`

  :   The confidence level for calculating confidence/prediction
      intervals. Default is 0.95.

  `interval`

  :   Type of interval to calculate:

      "none" (default)

      :   No interval to be calculated.

      "confidence"

      :   Calculate a confidence interval.

      "prediction"

      :   Calculate a prediction interval.

## Value

The data-frame with the following columns appended at the end

- .Richness:

  The richness (number of non-zero compositional variables) within each
  observation.

- .Evenness:

  The evenness (metric quantifying the relative abundance of each
  compositional variable) within each observation.

- .Gradient:

  An character string defining the diversity gradient used for averaging
  the response.

- .add_str_ID:

  An identifier column for grouping the cartesian product of all
  additional columns specified in \`add_var\` parameter (if \`add_var\`
  is specified).

- .Pred:

  The predicted response for each obsvervation.

- .Lower:

  The lower limit of the prediction/confidence interval for each
  observation.

- .Upper:

  The upper limit of the prediction/confidence interval for each
  observation.

- .Avg:

  The averaged value of the predicted response for each unique value of
  the selected diversity gradient.

## Details

Currently two diversity gradients are supported

- **Richness**: A metric describing the number of non-zero compositional
  variables in an observation.

- **Evenness**: A metric quantifying the relative abundances of all
  compositional variables in an observation. Defined as \$\$(2s/(s-1))
  \sum\_{i, j = 1; i \< j}^{s}{p_i \* p_j}\$\$ where \\s\\ is the total
  number of compositional variables and \\p_i\\ and \\p_j\\ are the
  proportions of the variables \\i\\ and \\j\\. See Kirwan et al., 2007
  \<[doi:10.1890/08-1684.1](https://doi.org/10.1890/08-1684.1) \> and
  Kirwan et al., 2009
  \<[doi:10.1890/08-1684.1](https://doi.org/10.1890/08-1684.1) \> for
  more information.

Here's a small example of how these metrics are calculated for a few
observations. Suppose we have four compositional variables (i.e. \\s =
4\\) and have the following three observations

- A = (0.5, 0.5, 0, 0)

- B = (0.25, 0.25, 0.25, 0.25)

- C = (1, 0, 0, 0)

The richness values for these three observations would be as follows

- A = 2 (Since two of the four compositional variables were non-zero)

- B = 4 (Since all four compositional variables were non-zero)

- C = 1 (Since one of the four compositional variables were non-zero)

The evenness values would be calculated as follows

- A = \\2\*4/(4-1)\*(0.5\*0.5+0.5\*0+0.5\*0+0.5\*0+0.5\*0+0\*0) = 0.67\\

- B =
  \\2\*4/(4-1)\*(0.25\*0.25+0.25\*0.25+0..25\*0.25+0.25\*0.25+0.25\*0.25+0.25\*0)
  = 1\\

- C = \\2\*4/(4-1)\*(1\*0+1\*0+1\*0+0\*0+0\*0+0\*0) = 0\\

## Examples

``` r
library(DImodels)
library(dplyr)

## Load data
data(sim2)

## Fit model
mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4)^2, data = sim2)

## Create data
## By default response would be averaged on the basis of richness
head(gradient_change_data(data = sim2,
                          prop = c("p1", "p2", "p3", "p4"),
                          model = mod))
#> ✔ Finished data preparation
#> # A tibble: 6 × 12
#>   community block    p1    p2    p3    p4 response .Richness .Evenness .Gradient
#>       <int> <fct> <dbl> <dbl> <dbl> <dbl>    <dbl>     <dbl>     <dbl> <chr>    
#> 1         1 1       0.7   0.1   0.1   0.1     20.2         4      0.64 .Richness
#> 2         1 2       0.7   0.1   0.1   0.1     20.1         4      0.64 .Richness
#> 3         1 3       0.7   0.1   0.1   0.1     20.9         4      0.64 .Richness
#> 4         1 4       0.7   0.1   0.1   0.1     17.0         4      0.64 .Richness
#> 5         2 1       0.1   0.7   0.1   0.1     17.2         4      0.64 .Richness
#> 6         2 2       0.1   0.7   0.1   0.1     19.9         4      0.64 .Richness
#> # ℹ 2 more variables: .Pred <dbl>, .Avg <dbl>

## Average response with respect to evenness
head(gradient_change_data(data = sim2,
                          prop = c("p1", "p2", "p3", "p4"),
                          model = mod,
                          gradient = "evenness"))
#> ✔ Finished data preparation
#> # A tibble: 6 × 12
#>   community block    p1    p2    p3    p4 response .Richness .Evenness .Gradient
#>       <int> <fct> <dbl> <dbl> <dbl> <dbl>    <dbl>     <dbl>     <dbl> <chr>    
#> 1         1 1       0.7   0.1   0.1   0.1     20.2         4      0.64 .Evenness
#> 2         1 2       0.7   0.1   0.1   0.1     20.1         4      0.64 .Evenness
#> 3         1 3       0.7   0.1   0.1   0.1     20.9         4      0.64 .Evenness
#> 4         1 4       0.7   0.1   0.1   0.1     17.0         4      0.64 .Evenness
#> 5         2 1       0.1   0.7   0.1   0.1     17.2         4      0.64 .Evenness
#> 6         2 2       0.1   0.7   0.1   0.1     19.9         4      0.64 .Evenness
#> # ℹ 2 more variables: .Pred <dbl>, .Avg <dbl>

## Additional variables can also be added to the data by either specifying
## them directly in the `data` or by using the `add_var` argument
## Refit model
sim2$block <- as.numeric(sim2$block)
new_mod <- update(mod, ~. + block, data = sim2)
## This model has block so we can either specify block in the data
subset_data <- sim2[c(1,5,9,11), 2:6]
subset_data
#>    block  p1  p2  p3  p4
#> 1      1 0.7 0.1 0.1 0.1
#> 5      1 0.1 0.7 0.1 0.1
#> 9      1 0.1 0.1 0.7 0.1
#> 11     3 0.1 0.1 0.7 0.1
head(gradient_change_data(data = subset_data,
                          prop = c("p1", "p2", "p3", "p4"),
                          model = mod,
                          gradient = "evenness"))
#> ✔ Finished data preparation
#> # A tibble: 4 × 10
#>   block    p1    p2    p3    p4 .Richness .Evenness .Gradient .Pred  .Avg
#>   <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>     <dbl> <chr>     <dbl> <dbl>
#> 1     1   0.7   0.1   0.1   0.1         4      0.64 .Evenness  18.4  17.3
#> 2     1   0.1   0.7   0.1   0.1         4      0.64 .Evenness  17.5  17.3
#> 3     1   0.1   0.1   0.7   0.1         4      0.64 .Evenness  16.6  17.3
#> 4     3   0.1   0.1   0.7   0.1         4      0.64 .Evenness  16.6  17.3
## Or we could add the variable using `add_var`
subset_data <- sim2[c(1,5,9,11), 3:6]
subset_data
#>     p1  p2  p3  p4
#> 1  0.7 0.1 0.1 0.1
#> 5  0.1 0.7 0.1 0.1
#> 9  0.1 0.1 0.7 0.1
#> 11 0.1 0.1 0.7 0.1
head(gradient_change_data(data = subset_data,
                          prop = c("p1", "p2", "p3", "p4"),
                          model = new_mod,
                          gradient = "evenness",
                          add_var = list(block = c(1, 2))))
#> ✔ Finished data preparation
#> # A tibble: 6 × 11
#>      p1    p2    p3    p4 block .add_str_ID .Richness .Evenness .Gradient .Pred
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <chr>           <dbl>     <dbl> <chr>     <dbl>
#> 1   0.7   0.1   0.1   0.1     1 block: 1            4      0.64 .Evenness  19.1
#> 2   0.1   0.7   0.1   0.1     1 block: 1            4      0.64 .Evenness  18.2
#> 3   0.1   0.1   0.7   0.1     1 block: 1            4      0.64 .Evenness  17.3
#> 4   0.1   0.1   0.7   0.1     1 block: 1            4      0.64 .Evenness  17.3
#> 5   0.7   0.1   0.1   0.1     2 block: 2            4      0.64 .Evenness  18.7
#> 6   0.1   0.7   0.1   0.1     2 block: 2            4      0.64 .Evenness  17.8
#> # ℹ 1 more variable: .Avg <dbl>
## The benefit of specifying the variable this way is we have an ID
## columns now called `.add_str_ID` which could be used to create a
## separate plot for each value of the additional variable


## Model coefficients can also be used, but then user would have
## to specify the data with all columns corresponding to each coefficient
coef_data <- sim2 %>%
               mutate(`p1:p2` = p1*p2, `p1:p3` = p1*p2, `p1:p4` = p1*p4,
                      `p2:p3` = p2*p3, `p2:p4` = p2*p4, `p3:p4` = p3*p4) %>%
               select(p1, p2, p3, p4,
                      `p1:p2`, `p1:p3`, `p1:p4`,
                      `p2:p3`, `p2:p4`, `p3:p4`) %>%
               slice(1,5,9,11)
print(coef_data)
#>     p1  p2  p3  p4 p1:p2 p1:p3 p1:p4 p2:p3 p2:p4 p3:p4
#> 1  0.7 0.1 0.1 0.1  0.07  0.07  0.07  0.01  0.01  0.01
#> 5  0.1 0.7 0.1 0.1  0.07  0.07  0.01  0.07  0.07  0.01
#> 9  0.1 0.1 0.7 0.1  0.01  0.01  0.01  0.07  0.01  0.07
#> 11 0.1 0.1 0.7 0.1  0.01  0.01  0.01  0.07  0.01  0.07
print(mod$coefficients)
#>        p1        p2        p3        p4     p1:p2     p1:p3     p1:p4     p2:p3 
#> 10.699426 10.228917  8.939289  8.532857 33.894874 37.552444 32.720996 26.739691 
#>     p2:p4     p3:p4 
#> 33.188799 27.771368 
gradient_change_data(data = coef_data,
                     prop = c("p1", "p2", "p3", "p4"),
                     gradient = "evenness",
                     coefficients = mod$coefficients,
                     interval = "none")
#> ✔ Finished data preparation
#> # A tibble: 4 × 15
#>      p1    p2    p3    p4 `p1:p2` `p1:p3` `p1:p4` `p2:p3` `p2:p4` `p3:p4`
#>   <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1   0.7   0.1   0.1   0.1    0.07    0.07    0.07    0.01    0.01    0.01
#> 2   0.1   0.7   0.1   0.1    0.07    0.07    0.01    0.07    0.07    0.01
#> 3   0.1   0.1   0.7   0.1    0.01    0.01    0.01    0.07    0.01    0.07
#> 4   0.1   0.1   0.7   0.1    0.01    0.01    0.01    0.07    0.01    0.07
#> # ℹ 5 more variables: .Richness <dbl>, .Evenness <dbl>, .Gradient <chr>,
#> #   .Pred <dbl>, .Avg <dbl>
```
