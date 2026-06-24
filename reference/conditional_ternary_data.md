# Data preparation for conditional ternary diagrams

The helper function for preparing the underlying data for creating
conditional ternary diagrams, where we fix \\n-3\\ variables to have a
constant value \\p_1, p_2, p_3, ..., p\_{n-3}\\ such that \\P = p_1 +
p_2 + p_3 + ... p\_{n - 3}\\ and \\0 \le P \le 1\\ and vary the
proportion of the remaining three variables between \\0\\ and \\1-P\\ to
visualise the change in the predicted response as a contour map within a
ternary diagram. The output of this function can be passed to the
[`conditional_ternary_plot`](https://rishvish.github.io/DImodelsVis/reference/conditional_ternary_plot.md)
function to plot the results. Viewing multiple 2-d slices across
multiple variables should allow to create an approximation of how the
response varies across the n-dimensional simplex.

## Usage

``` r
conditional_ternary_data(
  prop,
  FG = NULL,
  values = NULL,
  tern_vars = NULL,
  conditional = NULL,
  add_var = list(),
  resolution = 3,
  prediction = TRUE,
  ...
)
```

## Arguments

- prop:

  A character vector indicating the model coefficients corresponding to
  variable proportions. These variables should be compositional in
  nature (i.e., proportions should sum to 1).

- FG:

  A character vector specifying the grouping of the variables specified
  in \`prop\`. Specifying this parameter would call the
  grouped_ternary_data function internally. See
  [`grouped_ternary`](https://rishvish.github.io/DImodelsVis/reference/grouped_ternary.md)
  or
  [`grouped_ternary_data`](https://rishvish.github.io/DImodelsVis/reference/grouped_ternary_data.md)
  for more information.

- values:

  A numeric vector specifying the proportional split of the variables
  within a group. The default is to split the group proportion equally
  between each variable in the group.

- tern_vars:

  A character vector giving the names of the three variables to be shown
  in the ternary diagram.

- conditional:

  A data-frame describing the names of the compositional variables and
  their respective values at which to slice the simplex space. The
  format should be, for example, as follows:  
  `data.frame("p1" = c(0, 0.5), "p2" = c(0.2, 0.1))`  
  One figure would be created for each row in \`conditional\` with the
  respective values of all specified variables. Any compositional
  variables not specified in \`conditional\` will be assumed to be 0.

- add_var:

  A list or data-frame specifying values for additional variables in the
  model other than the proportions (i.e. not part of the simplex
  design). This could be useful for comparing the predictions across
  different values for a non-compositional variable. If specified as a
  list, it will be expanded to show a plot for each unique combination
  of values specified, while if specified as a data-frame, one plot
  would be generated for each row in the data.

- resolution:

  A number between 1 and 10 describing the resolution of the resultant
  graph. A high value would result in a higher definition figure but at
  the cost of being computationally expensive.

- prediction:

  A logical value indicating whether to pass the final data to the
  \`[add_prediction](https://rishvish.github.io/DImodelsVis/reference/add_prediction.md)\`
  function and append the predictions to the data. Default value is
  `TRUE`, but often it would be desirable to make additional changes to
  the data before making any predictions, so the user can set this to
  `FALSE` and manually call the
  \`[add_prediction](https://rishvish.github.io/DImodelsVis/reference/add_prediction.md)\`
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

A data-frame containing compositional columns with names specified in
\`prop\` parameter along with any additional columns specified in
\`add_var\` parameter. The first five columns of the data contain the
three variables (specified in \`tern_vars\`) shown in the ternary along
with their 2-d projection and should not be modified. The following
additional columns could also be present in the data.

- .x:

  The x-projection of the points within the ternary.

- .y:

  The y-projection of the points within the ternary.

- .add_str_ID:

  An identifier column for grouping the cartesian product of all
  additional columns specified in \`add_var\` parameter (if \`add_var\`
  is specified).

- .Sp:

  An identifier column specifying the variable(s) along which the high
  dimensional simplex is sliced.

- .Value:

  The value(s) (between 0 and 1) along the direction of variable(s) in
  \`.Sp\` at which the high dimensional simplex is sliced.

- .Facet:

  An identifier column formed by combining \`.Sp\` and \`.value\` to
  group observations within a specific slice of the high dimensional
  simplex.

- .Pred:

  The predicted response for each observation (if \`prediction\` is
  `TRUE`).

- .Lower:

  The lower limit of the prediction/confidence interval for each
  observation.

- .Upper:

  The upper limit of the prediction/confidence interval for each
  observation.

## Examples

``` r
library(DImodels)

## Load data
data(sim4)

## Fit model
mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4 + p5 + p6)^2, data = sim4)

## Create data
## Any species not specified in `tern_vars` or conditional will be assumed
## to be 0, for example p5 and p6 here.
head(conditional_ternary_data(prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
                              tern_vars = c("p1", "p2", "p3"),
                              conditional = data.frame("p4" = c(0, 0.2, 0.5)),
                              model = mod,
                              resolution = 1))
#> ✔ Finished data preparation.
#>   p1        p2          p3          .x .y p4 p5 p6 .Sp .Value .Facet    .Pred
#> 1  0 1.0000000 0.000000000 0.000000000  0  0  0  0  p4      0 p4 = 0 19.82837
#> 2  0 0.9949749 0.005025126 0.005025126  0  0  0  0  p4      0 p4 = 0 19.88595
#> 3  0 0.9899497 0.010050251 0.010050251  0  0  0  0  p4      0 p4 = 0 19.94305
#> 4  0 0.9849246 0.015075377 0.015075377  0  0  0  0  p4      0 p4 = 0 19.99966
#> 5  0 0.9798995 0.020100503 0.020100503  0  0  0  0  p4      0 p4 = 0 20.05579
#> 6  0 0.9748744 0.025125628 0.025125628  0  0  0  0  p4      0 p4 = 0 20.11143

## Can also condition on multiple species
cond <- data.frame(p4 = c(0, 0.2), p5 = c(0.5, 0.1), p6 = c(0, 0.3))
cond
#>    p4  p5  p6
#> 1 0.0 0.5 0.0
#> 2 0.2 0.1 0.3
head(conditional_ternary_data(prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
                              tern_vars = c("p1", "p2", "p3"),
                              conditional = cond,
                              model = mod,
                              resolution = 1))
#> ✔ Finished data preparation.
#>   p1        p2          p3          .x .y p4  p5 p6        .Sp    .Value
#> 1  0 0.5000000 0.000000000 0.000000000  0  0 0.5  0 p4, p5, p6 0, 0.5, 0
#> 2  0 0.4974874 0.002512563 0.005025126  0  0 0.5  0 p4, p5, p6 0, 0.5, 0
#> 3  0 0.4949749 0.005025126 0.010050251  0  0 0.5  0 p4, p5, p6 0, 0.5, 0
#> 4  0 0.4924623 0.007537688 0.015075377  0  0 0.5  0 p4, p5, p6 0, 0.5, 0
#> 5  0 0.4899497 0.010050251 0.020100503  0  0 0.5  0 p4, p5, p6 0, 0.5, 0
#> 6  0 0.4874372 0.012562814 0.025125628  0  0 0.5  0 p4, p5, p6 0, 0.5, 0
#>                     .Facet    .Pred
#> 1 p4 = 0; p5 = 0.5; p6 = 0 21.46473
#> 2 p4 = 0; p5 = 0.5; p6 = 0 21.50009
#> 3 p4 = 0; p5 = 0.5; p6 = 0 21.53533
#> 4 p4 = 0; p5 = 0.5; p6 = 0 21.57045
#> 5 p4 = 0; p5 = 0.5; p6 = 0 21.60544
#> 6 p4 = 0; p5 = 0.5; p6 = 0 21.64031

## Fit model
mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4 + p5 + p6)^2 + treatment,
           data = sim4)

## Can also add any additional variables independent of the simplex
## Notice the additional `.add_str_ID` column
head(conditional_ternary_data(prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
                              tern_vars = c("p1", "p2", "p3"),
                              conditional = data.frame("p4" = c(0, 0.2, 0.5)),
                              add_var = list("treatment" = c(50, 150)),
                              model = mod,
                              resolution = 1))
#> ✔ Finished data preparation.
#>   p1        p2          p3          .x .y p4 treatment   .add_str_ID p5 p6 .Sp
#> 1  0 1.0000000 0.000000000 0.000000000  0  0        50 treatment: 50  0  0  p4
#> 2  0 0.9949749 0.005025126 0.005025126  0  0        50 treatment: 50  0  0  p4
#> 3  0 0.9899497 0.010050251 0.010050251  0  0        50 treatment: 50  0  0  p4
#> 4  0 0.9849246 0.015075377 0.015075377  0  0        50 treatment: 50  0  0  p4
#> 5  0 0.9798995 0.020100503 0.020100503  0  0        50 treatment: 50  0  0  p4
#> 6  0 0.9748744 0.025125628 0.025125628  0  0        50 treatment: 50  0  0  p4
#>   .Value .Facet    .Pred
#> 1      0 p4 = 0 16.90156
#> 2      0 p4 = 0 16.95914
#> 3      0 p4 = 0 17.01624
#> 4      0 p4 = 0 17.07285
#> 5      0 p4 = 0 17.12898
#> 6      0 p4 = 0 17.18462

## It could be desirable to take the output of this function and add
## additional variables to the data before making predictions
## Use `prediction = FALSE` to get data without any predictions
cond_data <- conditional_ternary_data(prop = c("p1", "p2", "p3", "p4", "p5", "p6"),
                                      tern_vars = c("p1", "p2", "p3"),
                                      conditional = data.frame("p4" = c(0, 0.2, 0.5)),
                                      prediction = FALSE,
                                      resolution = 1)
#> ✔ Finished data preparation.
## The data can then be modified and the `add_prediction` function can be
## called manually using either the model object or model coefficients
cond_data$treatment <- 50
head(add_prediction(data = cond_data, model = mod))
#>   p1        p2          p3          .x .y p4 p5 p6 .Sp .Value .Facet treatment
#> 1  0 1.0000000 0.000000000 0.000000000  0  0  0  0  p4      0 p4 = 0        50
#> 2  0 0.9949749 0.005025126 0.005025126  0  0  0  0  p4      0 p4 = 0        50
#> 3  0 0.9899497 0.010050251 0.010050251  0  0  0  0  p4      0 p4 = 0        50
#> 4  0 0.9849246 0.015075377 0.015075377  0  0  0  0  p4      0 p4 = 0        50
#> 5  0 0.9798995 0.020100503 0.020100503  0  0  0  0  p4      0 p4 = 0        50
#> 6  0 0.9748744 0.025125628 0.025125628  0  0  0  0  p4      0 p4 = 0        50
#>      .Pred
#> 1 16.90156
#> 2 16.95914
#> 3 17.01624
#> 4 17.07285
#> 5 17.12898
#> 6 17.18462
```
