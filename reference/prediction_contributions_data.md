# Model term contributions to predicted response

The helper function for preparing the data to split the predicted
response from a regression model into contributions (predictor
coefficient \* predictor value) by the terms in the model. The output of
this function can be passed to the
\`[prediction_contributions_plot](https://rishvish.github.io/DImodelsVis/reference/prediction_contributions_plot.md)\`
function to visualise the results.

## Usage

``` r
prediction_contributions_data(
  data,
  model = NULL,
  coefficients = NULL,
  coeff_cols = NULL,
  vcov = NULL,
  add_var = list(),
  groups = list(),
  conf.level = 0.95,
  interval = c("confidence", "prediction", "none"),
  bar_labs = rownames(data)
)
```

## Arguments

- data:

  A user-defined data-frame containing values for compositional
  variables along with any additional variables that the user wishes to
  predict for. If left blank, a selection of observations (2 from each
  level of richness) from the original data used to fit the model would
  be selected.

- model:

  A Diversity Interactions model object fit by using the
  [`DI()`](https://rdrr.io/pkg/DImodels/man/DI.html) function from the
  [`DImodels`](https://rdrr.io/pkg/DImodels/man/DImodels-package.html)
  package.

- coefficients:

  If a regression model is not available (or can't be fit in R), the
  regression coefficients from a model fit in some other language can be
  used to calculate predictions. However, the user would have to ensure
  there's an appropriate one-to-one positional mapping between the data
  columns and the coefficient values. Further, they would also have to
  provide a variance-covariance matrix of the coefficients in the
  \`vcov\` parameter if they want the associated CI for the prediction
  or it would not be possible to calculate confidence/prediction
  intervals using this method.

- coeff_cols:

  If \`coefficients\` are specified and a one-to-one positional mapping
  between the data-columns and coefficient vector is not present. A
  character string or numeric index can be specified here to reorder the
  data columns and match the corresponding coefficient value to the
  respective data column. See the "Use model coefficients for
  prediction" section in examples.

- vcov:

  If regression coefficients are specified, then the variance-covariance
  matrix of the coefficients can be specified here to calculate the
  associated confidence interval around each prediction. Failure to do
  so would result in no confidence intervals being returned. Ensure
  \`coefficients\` and \`vcov\` have the same positional mapping with
  the data.

- add_var:

  A list specifying values for additional predictor variables in the
  model independent of the compositional predictor variables. This could
  be useful for comparing the predictions across different values for a
  non-compositional variable. If specified as a list, it will be
  expanded to show a plot for each unique combination of values
  specified, while if specified as a data-frame, one plot would be
  generated for each row in the data and they will be arranged in a grid
  according to the value specified in \`nrow\` and \`ncol\`.

- groups:

  A list specifying groupings to arrange coefficients into. The
  coefficients within a group will be added together and shown as a
  single component on the respective bars in the plot. This could be
  useful for grouping multiple similar terms into a single term for
  better visibility.

- conf.level:

  The confidence level for calculating confidence or prediction
  intervals.

- interval:

  Type of interval to calculate:

  "none"

  :   No interval to be calculated.

  "confidence" (default)

  :   Calculate a confidence interval.

  "prediction"

  :   Calculate a prediction interval.

- bar_labs:

  The labels to be shown for each bar in the plot. The user has three
  options: - By default, the row-names in the data would be used as
  labels for the bars. - A character string or numeric index indicating
  an ID column in data. - A character vector of same length as the
  number of rows in the data, which manually specifies the names for
  each bar. If none of the three options are available, the function
  would assign a unique ID for each bar.

## Value

A data-frame with the following columns. Any additional columns which
weren't used when fitting the model would also be present.

- .Community:

  An identifier column to discern each observation in the data. These
  are the labels which will be displayed for the bars in the plot.

- .add_str_ID:

  An identifier column for grouping the cartesian product of all
  additional columns specified in \`add_var\` parameter (if \`add_var\`
  is specified).

- .Pred:

  The predicted repsonse for each observation.

- .Lower:

  The lower limit of the prediction interval for each observation.

- .Upper:

  The lower limit of the prediction interval for each observation.

- .Contributions:

  An identifier describing the name of the coefficient contributing to
  the response.

- .Value:

  The contributed value of the respective coefficient/group to the total
  prediction.

## Examples

``` r
library(DImodels)
library(dplyr)

## Load data
data(sim2)

## Fit model
mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4)^2, data = sim2)

prediction_contributions_data(data = sim2[c(1,5,9,11), ],
                              model = mod)
#> ✔ Finished data preparation.
#> # A tibble: 40 × 10
#>    .x_labs .Community  community block response .Pred .Lower .Upper
#>    <chr>   <fct>           <int> <fct>    <dbl> <dbl>  <dbl>  <dbl>
#>  1 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  2 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  3 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  4 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  5 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  6 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  7 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  8 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  9 1       Community 1         1 1         20.2  18.4   17.6   19.3
#> 10 1       Community 1         1 1         20.2  18.4   17.6   19.3
#> # ℹ 30 more rows
#> # ℹ 2 more variables: .Contributions <chr>, .Value <dbl>

## Specific coefficients can also be grouped together
## Either by their indices in the model coefficient vector
prediction_contributions_data(data = sim2[c(1,5,9,11), ],
                              model = mod,
                              groups = list("Interactions" = 5:10))
#> ✔ Finished data preparation.
#> # A tibble: 20 × 10
#>    .x_labs .Community  community block response .Pred .Lower .Upper
#>    <chr>   <fct>           <int> <fct>    <dbl> <dbl>  <dbl>  <dbl>
#>  1 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  2 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  3 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  4 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  5 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  6 5       Community 2         2 1         17.2  17.5   16.7   18.4
#>  7 5       Community 2         2 1         17.2  17.5   16.7   18.4
#>  8 5       Community 2         2 1         17.2  17.5   16.7   18.4
#>  9 5       Community 2         2 1         17.2  17.5   16.7   18.4
#> 10 5       Community 2         2 1         17.2  17.5   16.7   18.4
#> 11 9       Community 3         3 1         17.9  16.6   15.8   17.5
#> 12 9       Community 3         3 1         17.9  16.6   15.8   17.5
#> 13 9       Community 3         3 1         17.9  16.6   15.8   17.5
#> 14 9       Community 3         3 1         17.9  16.6   15.8   17.5
#> 15 9       Community 3         3 1         17.9  16.6   15.8   17.5
#> 16 11      Community 4         3 3         17.6  16.6   15.8   17.5
#> 17 11      Community 4         3 3         17.6  16.6   15.8   17.5
#> 18 11      Community 4         3 3         17.6  16.6   15.8   17.5
#> 19 11      Community 4         3 3         17.6  16.6   15.8   17.5
#> 20 11      Community 4         3 3         17.6  16.6   15.8   17.5
#> # ℹ 2 more variables: .Contributions <chr>, .Value <dbl>
## Or by specifying the coefficient names as character strings
prediction_contributions_data(data = sim2[c(1,5,9,11), ],
                              model = mod,
                              groups = list("p1_Ints" = c("p1:p2",
                                                          "p1:p3",
                                                          "p1:p4")))
#> ✔ Finished data preparation.
#> # A tibble: 32 × 10
#>    .x_labs .Community  community block response .Pred .Lower .Upper
#>    <chr>   <fct>           <int> <fct>    <dbl> <dbl>  <dbl>  <dbl>
#>  1 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  2 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  3 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  4 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  5 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  6 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  7 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  8 1       Community 1         1 1         20.2  18.4   17.6   19.3
#>  9 5       Community 2         2 1         17.2  17.5   16.7   18.4
#> 10 5       Community 2         2 1         17.2  17.5   16.7   18.4
#> # ℹ 22 more rows
#> # ℹ 2 more variables: .Contributions <chr>, .Value <dbl>

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
head(prediction_contributions_data(data = subset_data,
                                   model = new_mod))
#> ✔ Finished data preparation.
#> # A tibble: 6 × 7
#>   .x_labs .Community  .Pred .Lower .Upper .Contributions .Value
#>   <chr>   <fct>       <dbl>  <dbl>  <dbl> <chr>           <dbl>
#> 1 1       Community 1  19.1   18.2   20.1 p1              8.29 
#> 2 1       Community 1  19.1   18.2   20.1 p2              1.14 
#> 3 1       Community 1  19.1   18.2   20.1 p3              1.01 
#> 4 1       Community 1  19.1   18.2   20.1 p4              0.967
#> 5 1       Community 1  19.1   18.2   20.1 block          -0.456
#> 6 1       Community 1  19.1   18.2   20.1 p1:p2           2.37 
## Or we could add the variable using `add_var`
subset_data <- sim2[c(1,5,9,11), 3:6]
subset_data
#>     p1  p2  p3  p4
#> 1  0.7 0.1 0.1 0.1
#> 5  0.1 0.7 0.1 0.1
#> 9  0.1 0.1 0.7 0.1
#> 11 0.1 0.1 0.7 0.1
head(prediction_contributions_data(data = subset_data,
                                   model = new_mod,
                                   add_var = list(block = c(1, 2))))
#> ✔ Finished data preparation.
#> # A tibble: 6 × 8
#>   .x_labs .Community  .add_str_ID .Pred .Lower .Upper .Contributions .Value
#>   <chr>   <fct>       <chr>       <dbl>  <dbl>  <dbl> <chr>           <dbl>
#> 1 1       Community 1 block: 1     19.1   18.2   20.1 p1              8.29 
#> 2 1       Community 1 block: 1     19.1   18.2   20.1 p2              1.14 
#> 3 1       Community 1 block: 1     19.1   18.2   20.1 p3              1.01 
#> 4 1       Community 1 block: 1     19.1   18.2   20.1 p4              0.967
#> 5 1       Community 1 block: 1     19.1   18.2   20.1 block          -0.456
#> 6 1       Community 1 block: 1     19.1   18.2   20.1 p1:p2           2.37 
## The benefit of specifying the variable this way is we have an ID
## columns now called `.add_str_ID` which would be used to create a
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
prediction_contributions_data(data = coef_data,
                              coefficients = mod$coefficients,
                              interval = "none")
#> ✔ Finished data preparation.
#> # A tibble: 40 × 5
#>    .x_labs .Community  .Pred .Contributions .Value
#>    <chr>   <fct>       <dbl> <chr>           <dbl>
#>  1 1       Community 1  18.4 p1              7.49 
#>  2 1       Community 1  18.4 p2              1.02 
#>  3 1       Community 1  18.4 p3              0.894
#>  4 1       Community 1  18.4 p4              0.853
#>  5 1       Community 1  18.4 p1:p2           2.37 
#>  6 1       Community 1  18.4 p1:p3           2.63 
#>  7 1       Community 1  18.4 p1:p4           2.29 
#>  8 1       Community 1  18.4 p2:p3           0.267
#>  9 1       Community 1  18.4 p2:p4           0.332
#> 10 1       Community 1  18.4 p3:p4           0.278
#> # ℹ 30 more rows
## To get uncertainity using coefficients vcov matrix would have to specified
prediction_contributions_data(data = coef_data,
                              coefficients = mod$coefficients,
                              vcov = vcov(mod))
#> ✔ Finished data preparation.
#> # A tibble: 40 × 7
#>    .x_labs .Community  .Pred .Lower .Upper .Contributions .Value
#>    <chr>   <fct>       <dbl>  <dbl>  <dbl> <chr>           <dbl>
#>  1 1       Community 1  18.4   17.6   19.3 p1              7.49 
#>  2 1       Community 1  18.4   17.6   19.3 p2              1.02 
#>  3 1       Community 1  18.4   17.6   19.3 p3              0.894
#>  4 1       Community 1  18.4   17.6   19.3 p4              0.853
#>  5 1       Community 1  18.4   17.6   19.3 p1:p2           2.37 
#>  6 1       Community 1  18.4   17.6   19.3 p1:p3           2.63 
#>  7 1       Community 1  18.4   17.6   19.3 p1:p4           2.29 
#>  8 1       Community 1  18.4   17.6   19.3 p2:p3           0.267
#>  9 1       Community 1  18.4   17.6   19.3 p2:p4           0.332
#> 10 1       Community 1  18.4   17.6   19.3 p3:p4           0.278
#> # ℹ 30 more rows

## Specifying `bar_labs`
## Our data has four rows so we'd need four labels in bar_labs
prediction_contributions_data(data = coef_data,
                              coefficients = mod$coefficients,
                              vcov = vcov(mod),
                              bar_labs = c("p1 Domm", "p2 Domm",
                                           "p3 Domm", "p4 Domm"))
#> ✔ Finished data preparation.
#> # A tibble: 40 × 7
#>    .x_labs .Community  .Pred .Lower .Upper .Contributions .Value
#>    <chr>   <fct>       <dbl>  <dbl>  <dbl> <chr>           <dbl>
#>  1 p1 Domm Community 1  18.4   17.6   19.3 p1              7.49 
#>  2 p1 Domm Community 1  18.4   17.6   19.3 p2              1.02 
#>  3 p1 Domm Community 1  18.4   17.6   19.3 p3              0.894
#>  4 p1 Domm Community 1  18.4   17.6   19.3 p4              0.853
#>  5 p1 Domm Community 1  18.4   17.6   19.3 p1:p2           2.37 
#>  6 p1 Domm Community 1  18.4   17.6   19.3 p1:p3           2.63 
#>  7 p1 Domm Community 1  18.4   17.6   19.3 p1:p4           2.29 
#>  8 p1 Domm Community 1  18.4   17.6   19.3 p2:p3           0.267
#>  9 p1 Domm Community 1  18.4   17.6   19.3 p2:p4           0.332
#> 10 p1 Domm Community 1  18.4   17.6   19.3 p3:p4           0.278
#> # ℹ 30 more rows
```
