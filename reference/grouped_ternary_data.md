# Prepare data for grouped ternary diagrams

The helper function for preparing the underlying data for creating
grouped ternary diagrams where the proportions of the compositional
variables are combined into groups and visualised on a ternary diagram.
These are very useful when we have multiple compositional variables that
can be grouped together by some hierarchical grouping structure. For
example, grouping species in a ecosystem based on the functions they
perform, or grouping political parties based on their national
alliances. Grouping variables this way allows us to reduce the
dimensionality of the compositional data and visualise it. This is akin
to looking at a 2-d slice of the high dimensional simplex. The relative
proportions of each variable within a group can be adjust to look at
different slices of the simplex. Looking at multiple such slices would
enable us to create an approximation of how the response varies across
the original n-dimensional simplex. The output of this function can be
passed to the
[`grouped_ternary_plot`](https://rishvish.github.io/DImodelsVis/reference/grouped_ternary_plot.md)
function to plot the results.

## Usage

``` r
grouped_ternary_data(
  prop,
  FG,
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

  A character vector specifying the groupings of the variables specified
  in \`prop\`.

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
\`FG\` and \`prop\` parameters along with any additional columns
specified in \`add_var\` parameter and the following columns appended at
the end.

- .x:

  The x-projection of the points within the ternary.

- .y:

  The y-projection of the points within the ternary.

- .add_str_ID:

  An identifier column for grouping the cartesian product of all
  additional columns specified in \`add_var\` parameter (if \`add_var\`
  is specified).

- .Sp:

  An identifier column specifying the functional group along which the
  high dimensional simplex is sliced (if there are more than 3 groups).

- .Value:

  The value (between 0 and 1) along the direction of functional group in
  \`.Sp\` at which the high dimensional simplex is sliced.

- .Facet:

  An identifier column formed by combining \`.Sp\` and \`.value\` to
  group observations within a specific slice of the high dimensional
  simplex.

- .Pred:

  The predicted response for each observation. (if \`prediction\` is
  `TRUE`)

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
data(sim3)

## Fit model
mod <- glm(response ~ 0 + (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9)^2,
           data = sim3)

## Create data
## We have nine (p1 to p9) variables here and using \code{\link{conditional_ternary}}
## to visualise the simplex space won't be very helpful as there are too
## variables to condition on
## Instead we group the nine-variables into three groups called "G", "L" and "H"
head(grouped_ternary_data(model = mod,
                          prop = paste0("p", 1:9),
                          FG = c("G","G","G","G","G","L","L","H","H"),
                          resolution = 1))
#> Warning: The proportional split of species in the groups was not specified in `values`.
#> Assuming an equal split for species in each group.
#> ✔ Finished data preparation.
#>   G         L           H          .x .y p1 p2 p3 p4 p5        p6        p7
#> 1 0 1.0000000 0.000000000 0.000000000  0  0  0  0  0  0 0.5000000 0.5000000
#> 2 0 0.9949749 0.005025126 0.005025126  0  0  0  0  0  0 0.4974874 0.4974874
#> 3 0 0.9899497 0.010050251 0.010050251  0  0  0  0  0  0 0.4949749 0.4949749
#> 4 0 0.9849246 0.015075377 0.015075377  0  0  0  0  0  0 0.4924623 0.4924623
#> 5 0 0.9798995 0.020100503 0.020100503  0  0  0  0  0  0 0.4899497 0.4899497
#> 6 0 0.9748744 0.025125628 0.025125628  0  0  0  0  0  0 0.4874372 0.4874372
#>            p8          p9    .Pred
#> 1 0.000000000 0.000000000 7.493868
#> 2 0.002512563 0.002512563 7.517016
#> 3 0.005025126 0.005025126 7.540053
#> 4 0.007537688 0.007537688 7.562977
#> 5 0.010050251 0.010050251 7.585788
#> 6 0.012562814 0.012562814 7.608488

## By default the variables within a group take up an equal share of the
## group proportion. So for example, each point along the above ternary
## would have a 50:50 split of the variables in the group "L" or "H", thus
## the vertex where "L" is 1, would mean that p6 and p7 are 0.5 each,
## similarly, the vertex "H" is made up of 0.5 of p8 and p9 while the "G"
## vertex is comprised of 0.2 of each of p1, p2, p3, p4, and p5. The concepts
## also extend to points along the edges and interior of the ternary.

## Change the proportional split of species within an FG by using `values`
## `values` takes a numeric vector where the position of each element
## describes the proportion of the corresponding species within the
## corresponding FG
## For examples this vector describes, 2-% each of p1, p2, p3, p4 and p5,
## in G, 0% and 100% of p6 and p7, respectively in G2 and 30% and 70% of
## p8 and p9, respectively in G3.
vals <- c(0.2, 0.2, 0.2, 0.2, 0.2,
          0, 1,
          0.3, 0.7)
head(grouped_ternary_data(prop = paste0("p", 1:9),
                          FG = c("G","G","G","G","G","L","L","H","H"),
                          values = vals,
                          resolution = 1,
                          model = mod))
#> ✔ Finished data preparation.
#>   G         L           H          .x .y p1 p2 p3 p4 p5 p6        p7
#> 1 0 1.0000000 0.000000000 0.000000000  0  0  0  0  0  0  0 1.0000000
#> 2 0 0.9949749 0.005025126 0.005025126  0  0  0  0  0  0  0 0.9949749
#> 3 0 0.9899497 0.010050251 0.010050251  0  0  0  0  0  0  0 0.9899497
#> 4 0 0.9849246 0.015075377 0.015075377  0  0  0  0  0  0  0 0.9849246
#> 5 0 0.9798995 0.020100503 0.020100503  0  0  0  0  0  0  0 0.9798995
#> 6 0 0.9748744 0.025125628 0.025125628  0  0  0  0  0  0  0 0.9748744
#>            p8          p9    .Pred
#> 1 0.000000000 0.000000000 7.110611
#> 2 0.001507538 0.003517588 7.135017
#> 3 0.003015075 0.007035176 7.159316
#> 4 0.004522613 0.010552764 7.183509
#> 5 0.006030151 0.014070352 7.207595
#> 6 0.007537688 0.017587940 7.231576

## Can also add any additional experimental structures
## Notice .add_str_ID in the data
head(grouped_ternary_data(prop = paste0("p", 1:9),
                          FG = c("G","G","G","G","G","L","L","H","H"),
                          add_var = list("treatment" = c("50", "150")),
                          values = vals,
                          model = mod,
                          resolution = 1))
#> ✔ Finished data preparation.
#>   G         L           H          .x .y p1 p2 p3 p4 p5 p6        p7
#> 1 0 1.0000000 0.000000000 0.000000000  0  0  0  0  0  0  0 1.0000000
#> 2 0 0.9949749 0.005025126 0.005025126  0  0  0  0  0  0  0 0.9949749
#> 3 0 0.9899497 0.010050251 0.010050251  0  0  0  0  0  0  0 0.9899497
#> 4 0 0.9849246 0.015075377 0.015075377  0  0  0  0  0  0  0 0.9849246
#> 5 0 0.9798995 0.020100503 0.020100503  0  0  0  0  0  0  0 0.9798995
#> 6 0 0.9748744 0.025125628 0.025125628  0  0  0  0  0  0  0 0.9748744
#>            p8          p9 treatment   .add_str_ID    .Pred
#> 1 0.000000000 0.000000000        50 treatment: 50 7.110611
#> 2 0.001507538 0.003517588        50 treatment: 50 7.135017
#> 3 0.003015075 0.007035176        50 treatment: 50 7.159316
#> 4 0.004522613 0.010552764        50 treatment: 50 7.183509
#> 5 0.006030151 0.014070352        50 treatment: 50 7.207595
#> 6 0.007537688 0.017587940        50 treatment: 50 7.231576

## It could be desirable to take the output of this function and add
## additional variables to the data before making predictions
## Use `prediction = FALSE` to get data without any predictions
grouped_data <- grouped_ternary_data(prop = paste0("p", 1:9),
                                     FG = c("G","G","G","G","G","L","L","H","H"),
                                     values = vals,
                                     resolution = 1,
                                     prediction = FALSE)
#> ✔ Finished data preparation.
grouped_data$treatment <- 250
# Add predictions
head(add_prediction(data = grouped_data, model = mod))
#>   G         L           H          .x .y p1 p2 p3 p4 p5 p6        p7
#> 1 0 1.0000000 0.000000000 0.000000000  0  0  0  0  0  0  0 1.0000000
#> 2 0 0.9949749 0.005025126 0.005025126  0  0  0  0  0  0  0 0.9949749
#> 3 0 0.9899497 0.010050251 0.010050251  0  0  0  0  0  0  0 0.9899497
#> 4 0 0.9849246 0.015075377 0.015075377  0  0  0  0  0  0  0 0.9849246
#> 5 0 0.9798995 0.020100503 0.020100503  0  0  0  0  0  0  0 0.9798995
#> 6 0 0.9748744 0.025125628 0.025125628  0  0  0  0  0  0  0 0.9748744
#>            p8          p9 treatment    .Pred
#> 1 0.000000000 0.000000000       250 7.110611
#> 2 0.001507538 0.003517588       250 7.135017
#> 3 0.003015075 0.007035176       250 7.159316
#> 4 0.004522613 0.010552764       250 7.183509
#> 5 0.006030151 0.014070352       250 7.207595
#> 6 0.007537688 0.017587940       250 7.231576
```
