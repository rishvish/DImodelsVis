# Prepare data for visualising change in response across points in the simplex space

This is the helper function to prepare the underlying data for
visualising the change in a response variable between two points in a
simplex space. The two points specified by the \`starts\` and \`ends\`
parameters are joined by a straight line across the simplex space and
the response is predicted for the starting, ending and intermediate
communities along this line. The associated uncertainty along this
prediction is also returned. The output of this function can be passed
to the
[`simplex_path_plot`](https://rishvish.github.io/DImodelsVis/reference/simplex_path_plot.md)
function to visualise the change in response.

## Usage

``` r
simplex_path_data(starts, ends, prop, add_var = list(), prediction = TRUE, ...)
```

## Arguments

- starts:

  A data-frame specifying the starting proportions of the compositional
  variables. If a model object is specified then this data should
  contain all the variables present in the model object including any
  additional non-compositional variables. If a coefficient vector is
  specified then data should contain same number of columns as the
  number of elements in the coefficient vector and a one-to-one
  positional mapping would be assumed between the data columns and the
  elements of the coefficient vector.

- ends:

  A data-frame specifying the ending proportions of the compositional
  variables. If a model object is specified then this data should
  contain all the variables present in the model object including any
  additional non-compositional variables. If a coefficient vector is
  specified then data should contain same number of columns as the
  number of elements in the coefficient vector and a one-to-one
  positional mapping would be assumed between the data columns and the
  elements of the coefficient vector.

- prop:

  A vector of column names identifying the columns containing the
  variable proportions (i.e., compositional columns) in the data.

- add_var:

  A list or data-frame specifying values for additional variables in the
  model other than the proportions (i.e. not part of the simplex
  design). This could be useful for comparing the predictions across
  different values for a non-compositional variable. If specified as a
  list, it will be expanded to show a plot for each unique combination
  of values specified, while if specified as a data-frame, one plot
  would be generated for each row in the data.

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

A data frame with the following columns appended at the end

- .InterpConst:

  The value of the interpolation constant for creating the intermediate
  compositions between the start and end compositions.

- .Group:

  An identifier column to discern between the different curves.

- .add_str_ID:

  An identifier column for grouping the cartesian product of all
  additional columns specified in \`add_var\` parameter (if \`add_var\`
  is specified).

- .Pred:

  The predicted response for each observation.

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
data(sim2)

## Fit model
mod <- glm(response ~ (p1 + p2 + p3 + p4)^2 + 0, data = sim2)

## Create data for visualising change in response as we move from
## a species dominated by 70% of one species to a monoculture of
## same species
head(simplex_path_data(starts = sim2[c(1, 5, 9, 13), 3:6],
                       ends = sim2[c(48, 52, 56, 60), 3:6],
                       prop = c("p1", "p2", "p3", "p4"),
                       model = mod))
#> ✔ Finished data preparation.
#>      p1    p2    p3    p4 .InterpConst .Group    .Pred   .Lower   .Upper
#> 1 0.700 0.100 0.100 0.100         0.00      1 18.42849 17.58232 19.27465
#> 2 0.703 0.099 0.099 0.099         0.01      1 18.37345 17.52663 19.22027
#> 3 0.706 0.098 0.098 0.098         0.02      1 18.31797 17.47045 19.16548
#> 4 0.709 0.097 0.097 0.097         0.03      1 18.26203 17.41378 19.11028
#> 5 0.712 0.096 0.096 0.096         0.04      1 18.20565 17.35661 19.05468
#> 6 0.715 0.095 0.095 0.095         0.05      1 18.14881 17.29895 18.99868

## Create data for visualising change in response as we move from
## the centroid mixture to each monoculture
## If either of starts or ends have only row, then they'll be recycled
## to match the number of rows in the other
## Notice starts has only one row here, but will be recycled to have 4
## since ends has 4 four rows
head(simplex_path_data(starts = sim2[c(18),3:6],
                       ends = sim2[c(48, 52, 56, 60),3:6],
                       prop = c("p1", "p2", "p3", "p4"),
                       model = mod))
#> ✔ Finished data preparation.
#>       p1     p2     p3     p4 .InterpConst .Group    .Pred   .Lower   .Upper
#> 1 0.2500 0.2500 0.2500 0.2500         0.00      1 21.59188 20.99819 22.18557
#> 2 0.2575 0.2475 0.2475 0.2475         0.01      1 21.62206 21.02781 22.21630
#> 3 0.2650 0.2450 0.2450 0.2450         0.02      1 21.64942 21.05356 22.24528
#> 4 0.2725 0.2425 0.2425 0.2425         0.03      1 21.67397 21.07553 22.27242
#> 5 0.2800 0.2400 0.2400 0.2400         0.04      1 21.69572 21.09380 22.29763
#> 6 0.2875 0.2375 0.2375 0.2375         0.05      1 21.71465 21.10848 22.32082

## Changing the confidence level for the prediction interval
## Use `conf.level` parameter
head(simplex_path_data(starts = sim2[c(18), 3:6],
                       ends = sim2[c(48, 52, 56, 60),3:6],
                       prop = c("p1", "p2", "p3", "p4"),
                       model = mod, conf.level = 0.99))
#> ✔ Finished data preparation.
#>       p1     p2     p3     p4 .InterpConst .Group    .Pred   .Lower   .Upper
#> 1 0.2500 0.2500 0.2500 0.2500         0.00      1 21.59188 20.80038 22.38338
#> 2 0.2575 0.2475 0.2475 0.2475         0.01      1 21.62206 20.82981 22.41430
#> 3 0.2650 0.2450 0.2450 0.2450         0.02      1 21.64942 20.85502 22.44382
#> 4 0.2725 0.2425 0.2425 0.2425         0.03      1 21.67397 20.87613 22.47182
#> 5 0.2800 0.2400 0.2400 0.2400         0.04      1 21.69572 20.89325 22.49818
#> 6 0.2875 0.2375 0.2375 0.2375         0.05      1 21.71465 20.90651 22.52279

## Adding additional variables to the data using `add_var`
## Notice the new .add_str_ID column in the output
sim2$block <- as.numeric(sim2$block)
new_mod <- update(mod, ~ . + block, data = sim2)
head(simplex_path_data(starts = sim2[c(18), 3:6],
                       ends = sim2[c(48, 52, 56, 60), 3:6],
                       prop = c("p1", "p2", "p3", "p4"),
                       model = new_mod, conf.level = 0.99,
                       add_var = list("block" = c(1, 2))))
#> ✔ Finished data preparation.
#>       p1     p2     p3     p4 .InterpConst .Group block .add_str_ID    .Pred
#> 1 0.2500 0.2500 0.2500 0.2500         0.00      1     1    block: 1 22.27542
#> 2 0.2575 0.2475 0.2475 0.2475         0.01      1     1    block: 1 22.30560
#> 3 0.2650 0.2450 0.2450 0.2450         0.02      1     1    block: 1 22.33296
#> 4 0.2725 0.2425 0.2425 0.2425         0.03      1     1    block: 1 22.35751
#> 5 0.2800 0.2400 0.2400 0.2400         0.04      1     1    block: 1 22.37926
#> 6 0.2875 0.2375 0.2375 0.2375         0.05      1     1    block: 1 22.39819
#>     .Lower   .Upper
#> 1 21.27206 23.27879
#> 2 21.30171 23.30948
#> 3 21.32756 23.33836
#> 4 21.34970 23.36533
#> 5 21.36819 23.39032
#> 6 21.38312 23.41326

## Use predict = FALSE to get raw data structure
out_data <- simplex_path_data(starts = sim2[c(18), 3:6],
                              ends = sim2[c(48, 52, 56, 60), 3:6],
                              prop = c("p1", "p2", "p3", "p4"),
                              model = new_mod,
                              prediction = FALSE)
#> ✔ Finished data preparation.
head(out_data)
#>       p1     p2     p3     p4 .InterpConst .Group
#> 1 0.2500 0.2500 0.2500 0.2500         0.00      1
#> 2 0.2575 0.2475 0.2475 0.2475         0.01      1
#> 3 0.2650 0.2450 0.2450 0.2450         0.02      1
#> 4 0.2725 0.2425 0.2425 0.2425         0.03      1
#> 5 0.2800 0.2400 0.2400 0.2400         0.04      1
#> 6 0.2875 0.2375 0.2375 0.2375         0.05      1
## Manually add block
out_data$block = 3
## Call `add_prediction` to get prediction
head(add_prediction(data = out_data, model = new_mod, interval = "conf"))
#>       p1     p2     p3     p4 .InterpConst .Group block    .Pred   .Lower
#> 1 0.2500 0.2500 0.2500 0.2500         0.00      1     3 21.36404 20.78032
#> 2 0.2575 0.2475 0.2475 0.2475         0.01      1     3 21.39421 20.80999
#> 3 0.2650 0.2450 0.2450 0.2450         0.02      1     3 21.42157 20.83590
#> 4 0.2725 0.2425 0.2425 0.2425         0.03      1     3 21.44613 20.85812
#> 5 0.2800 0.2400 0.2400 0.2400         0.04      1     3 21.46787 20.87674
#> 6 0.2875 0.2375 0.2375 0.2375         0.05      1     3 21.48680 20.89183
#>     .Upper
#> 1 21.94775
#> 2 21.97843
#> 3 22.00725
#> 4 22.03413
#> 5 22.05900
#> 6 22.08178
```
