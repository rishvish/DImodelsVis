# Prepare data for effects plots of compositional predictors

The helper function to create the underlying data for visualising the
effect of increasing or decreasing (or both) the proportion of a
variable from a set of compositional predictors. This is a special case
of the
[`simplex_path`](https://rishvish.github.io/DImodelsVis/reference/simplex_path.md)
function where the end points are either the monoculture (i.e. variable
of interest = 1, while all others equal 0) of the variable of interest
(when increasing the proportion) or a community without the variable of
interest (when decreasing the proportion). The observations specified in
\`data\` are connected to the respective communities (monoculture of the
variable of interest or the community without the variable of interest)
by a straight line across the simplex; This has the effect of changing
the proportion of the variable of interest whilst adjusting the
proportion of the other variables but keeping the ratio of their
relative proportions unchanged, thereby preserving the compositional
nature of the data. See examples for more information. The output of
this function can be passed to the
[`visualise_effects_plot`](https://rishvish.github.io/DImodelsVis/reference/visualise_effects_plot.md)
function to visualise the results.

## Usage

``` r
visualise_effects_data(
  data,
  prop,
  var_interest = NULL,
  effect = c("increase", "decrease", "both"),
  add_var = list(),
  prediction = TRUE,
  ...
)
```

## Arguments

- data:

  A dataframe specifying the initial communities of interest for which
  to visualise the effect of increasing/decreasing a variable. If a
  model object is specified then this data should contain all the
  variables present in the model object including any additional
  variables not part of the simplex design. If a coefficient vector is
  specified then data should contain same number of columns as the
  number of elements in the coefficient vector and a one-to-one
  positional mapping would be assumed between the data columns and the
  elements of the coefficient vector.

- prop:

  A vector of column names or indices identifying the columns containing
  the variable proportions (i.e., compositional columns) in the data.

- var_interest:

  A character vector specifying the variable for which to visualise the
  effect of change on the response. If left blank, all variables would
  be assumed to be of interest.

- effect:

  One of "increase", "decrease" or "both" to indicate whether to look at
  the effect of increasing the proportion, decreasing the proportion or
  doing both simultaneously, respectively on the response. The default
  in "increasing".

- add_var:

  A list specifying values for additional variables in the model other
  than the proportions (i.e. not part of the simplex design). This would
  be useful to compare the predictions across different values for a
  categorical variable. One plot will be generated for each unique
  combination of values specified here.

- prediction:

  A logical value indicating whether to pass the final data to
  \`[add_prediction](https://rishvish.github.io/DImodelsVis/reference/add_prediction.md)\`
  and add predictions to the data. Default value is `TRUE`, but often it
  would be desirable to make additional changes to the data before
  making any predictions, so the user can set this to `FALSE` and
  manually call the
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

- .Sp:

  An identifier column to discern the variable of interest being
  modified in each curve.

- .Proportion:

  The value of the variable of interest within the community.

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

- .Marginal:

  The marginal change in the response (first derivative) with respect to
  the gradual change in the proportion of the species of interest.

- .Threshold:

  A numeric value indicating the maximum proportion of the species of
  interest within a particular community which has a positive marginal
  effect on the response.

- .MarEffect:

  A character string entailing whether the increase/decrease of the
  species of interest from the particular community would result in a
  positive or negative marginal effect on the response.

- .Effect:

  An identifier column signifying whether considering the effect of
  species addition or species decrease.

## Examples

``` r
library(DImodels)

## Load data
data(sim1)

## Fit model
mod <- glm(response ~ p1 + p2 + p3 + p4 + 0, data = sim1)

## Create data for visualising effect of increasing the proportion of
## variable p1 in data
## Notice how the proportion of `p1` increases while the proportion of
## the other variables decreases whilst maintaining their relative proportions
head(visualise_effects_data(data = sim1, prop = c("p1", "p2", "p3", "p4"),
                            var_interest = "p1", effect = "increase",
                            model = mod))
#> ✔ Finished data preparation.
#>      p1    p2    p3    p4 community block response .Sp .Proportion .Group
#> 1 0.700 0.100 0.100 0.100         1     1   10.815  p1       0.700      1
#> 2 0.703 0.099 0.099 0.099         1     1   10.815  p1       0.703      1
#> 3 0.706 0.098 0.098 0.098         1     1   10.815  p1       0.706      1
#> 4 0.709 0.097 0.097 0.097         1     1   10.815  p1       0.709      1
#> 5 0.712 0.096 0.096 0.096         1     1   10.815  p1       0.712      1
#> 6 0.715 0.095 0.095 0.095         1     1   10.815  p1       0.715      1
#>      .Pred   .Lower   .Upper .Marginal .Threshold .MarEffect  .Effect
#> 1 10.42269 9.789966 11.05541  1.559416      0.826   Negative increase
#> 2 10.42737 9.791585 11.06315  1.559416      0.826   Negative increase
#> 3 10.43204 9.793199 11.07089  1.559416      0.826   Negative increase
#> 4 10.43672 9.794807 11.07864  1.559416      0.826   Negative increase
#> 5 10.44140 9.796410 11.08639  1.559416      0.826   Negative increase
#> 6 10.44608 9.798008 11.09415  1.559416      0.826   Negative increase

## Create data for visualising the effect of decreasing the proportion
## variable p1 in data using `effect = "decrease"`
head(visualise_effects_data(data = sim1, prop = c("p1", "p2", "p3", "p4"),
                            var_interest = "p1", effect = "decrease",
                            model = mod))
#> ✔ Finished data preparation.
#>      p1        p2        p3        p4 community block response .Sp .Proportion
#> 1 0.000 0.3333333 0.3333333 0.3333333         1     1   10.815  p1       0.000
#> 2 0.007 0.3310000 0.3310000 0.3310000         1     1   10.815  p1       0.007
#> 3 0.014 0.3286667 0.3286667 0.3286667         1     1   10.815  p1       0.014
#> 4 0.021 0.3263333 0.3263333 0.3263333         1     1   10.815  p1       0.021
#> 5 0.028 0.3240000 0.3240000 0.3240000         1     1   10.815  p1       0.028
#> 6 0.035 0.3216667 0.3216667 0.3216667         1     1   10.815  p1       0.035
#>   .Group    .Pred   .Lower   .Upper .Marginal .Threshold .MarEffect  .Effect
#> 1      1 9.331096 8.884095 9.778097  1.559416      0.287   Negative decrease
#> 2      1 9.342012 8.900575 9.783450  1.559416      0.287   Negative decrease
#> 3      1 9.352928 8.916964 9.788892  1.559416      0.287   Negative decrease
#> 4      1 9.363844 8.933260 9.794427  1.559416      0.287   Negative decrease
#> 5      1 9.374760 8.949459 9.800060  1.559416      0.287   Negative decrease
#> 6      1 9.385676 8.965558 9.805793  1.559416      0.287   Negative decrease

## Create data for visualising the effect of increasing and decreasing the
## proportion variable p3 in data using `effect = "both"`
head(visualise_effects_data(data = sim1, prop = c("p1", "p2", "p3", "p4"),
                            var_interest = "p3", effect = "decrease",
                            model = mod))
#> ✔ Finished data preparation.
#>          p1        p2    p3        p4 community block response .Sp .Proportion
#> 1 0.7777778 0.1111111 0.000 0.1111111         1     1   10.815  p3       0.000
#> 2 0.7770000 0.1110000 0.001 0.1110000         1     1   10.815  p3       0.001
#> 3 0.7762222 0.1108889 0.002 0.1108889         1     1   10.815  p3       0.002
#> 4 0.7754444 0.1107778 0.003 0.1107778         1     1   10.815  p3       0.003
#> 5 0.7746667 0.1106667 0.004 0.1106667         1     1   10.815  p3       0.004
#> 6 0.7738889 0.1105556 0.005 0.1105556         1     1   10.815  p3       0.005
#>   .Group    .Pred   .Lower   .Upper  .Marginal .Threshold .MarEffect  .Effect
#> 1      1 10.51542 9.795624 11.23522 -0.9273268      0.066   Negative decrease
#> 2      1 10.51449 9.795636 11.23335 -0.9273268      0.066   Negative decrease
#> 3      1 10.51357 9.795646 11.23148 -0.9273268      0.066   Negative decrease
#> 4      1 10.51264 9.795656 11.22962 -0.9273268      0.066   Negative decrease
#> 5      1 10.51171 9.795664 11.22776 -0.9273268      0.066   Negative decrease
#> 6      1 10.51078 9.795671 11.22590 -0.9273268      0.066   Negative decrease

## Getting prediction intervals at a 99% confidence level
head(visualise_effects_data(data = sim1, prop = c("p1", "p2", "p3", "p4"),
                            var_interest = "p1", effect = "decrease",
                            model = mod, conf.level = 0.99,
                            interval = "prediction"))
#> ✔ Finished data preparation.
#>      p1        p2        p3        p4 community block response .Sp .Proportion
#> 1 0.000 0.3333333 0.3333333 0.3333333         1     1   10.815  p1       0.000
#> 2 0.007 0.3310000 0.3310000 0.3310000         1     1   10.815  p1       0.007
#> 3 0.014 0.3286667 0.3286667 0.3286667         1     1   10.815  p1       0.014
#> 4 0.021 0.3263333 0.3263333 0.3263333         1     1   10.815  p1       0.021
#> 5 0.028 0.3240000 0.3240000 0.3240000         1     1   10.815  p1       0.028
#> 6 0.035 0.3216667 0.3216667 0.3216667         1     1   10.815  p1       0.035
#>   .Group    .Pred   .Lower   .Upper .Marginal .Threshold .MarEffect  .Effect
#> 1      1 9.331096 5.855633 12.80656  1.559416      0.287   Negative decrease
#> 2      1 9.342012 5.867809 12.81621  1.559416      0.287   Negative decrease
#> 3      1 9.352928 5.879950 12.82591  1.559416      0.287   Negative decrease
#> 4      1 9.363844 5.892055 12.83563  1.559416      0.287   Negative decrease
#> 5      1 9.374760 5.904125 12.84539  1.559416      0.287   Negative decrease
#> 6      1 9.385676 5.916160 12.85519  1.559416      0.287   Negative decrease

## Adding additional variables to the data using `add_var`
## Notice the new .add_str_ID column in the output
sim1$block <- as.numeric(sim1$block)
new_mod <- update(mod, ~ . + block, data = sim1)
head(visualise_effects_data(data = sim1[, 3:6], prop = c("p1", "p2", "p3", "p4"),
                            var_interest = "p1", effect = "both",
                            model = new_mod,
                            add_var = list("block" = c(1, 2))))
#> ✔ Finished data preparation.
#>      p1        p2        p3        p4 block .add_str_ID .Sp .Proportion .Group
#> 1 0.000 0.3333333 0.3333333 0.3333333     1    block: 1  p1       0.000      1
#> 2 0.014 0.3286667 0.3286667 0.3286667     1    block: 1  p1       0.014      1
#> 3 0.028 0.3240000 0.3240000 0.3240000     1    block: 1  p1       0.028      1
#> 4 0.042 0.3193333 0.3193333 0.3193333     1    block: 1  p1       0.042      1
#> 5 0.056 0.3146667 0.3146667 0.3146667     1    block: 1  p1       0.056      1
#> 6 0.070 0.3100000 0.3100000 0.3100000     1    block: 1  p1       0.070      1
#>      .Pred   .Lower   .Upper .Marginal .Threshold .MarEffect .Effect
#> 1 9.847626 9.241980 10.45327  1.559416      0.946   Negative    both
#> 2 9.869458 9.271266 10.46765  1.559416      0.946   Negative    both
#> 3 9.891290 9.300209 10.48237  1.559416      0.946   Negative    both
#> 4 9.913122 9.328796 10.49745  1.559416      0.946   Negative    both
#> 5 9.934953 9.357015 10.51289  1.559416      0.946   Negative    both
#> 6 9.956785 9.384853 10.52872  1.559416      0.946   Negative    both

## Create data for visualising effect of decreasing variable p2 from
## the original communities in the data but using model coefficients
## When specifying coefficients the data should have a one-to-one
## positional mapping with specified coefficients.
init_comms <- sim1[, c("p1", "p2", "p3", "p4")]
head(visualise_effects_data(data = init_comms, prop = 1:4,
                            var_interest = "p2",
                            effect = "decrease",
                            interval = "none",
                            coefficients = mod$coefficients))
#> ✔ Finished data preparation.
#>          p1    p2        p3        p4 .Sp .Proportion .Group    .Pred .Marginal
#> 1 0.7777778 0.000 0.1111111 0.1111111  p2       0.000      1 10.41016 0.1252706
#> 2 0.7770000 0.001 0.1110000 0.1110000  p2       0.001      1 10.41029 0.1252706
#> 3 0.7762222 0.002 0.1108889 0.1108889  p2       0.002      1 10.41041 0.1252706
#> 4 0.7754444 0.003 0.1107778 0.1107778  p2       0.003      1 10.41054 0.1252706
#> 5 0.7746667 0.004 0.1106667 0.1106667  p2       0.004      1 10.41066 0.1252706
#> 6 0.7738889 0.005 0.1105556 0.1105556  p2       0.005      1 10.41079 0.1252706
#>   .Threshold .MarEffect  .Effect
#> 1      0.067   Negative decrease
#> 2      0.067   Negative decrease
#> 3      0.067   Negative decrease
#> 4      0.067   Negative decrease
#> 5      0.067   Negative decrease
#> 6      0.067   Negative decrease

## Note that to get confidence interval when specifying
## model coefficients we'd also need to provide a variance covariance
## matrix using the `vcov` argument
head(visualise_effects_data(data = init_comms, prop = 1:4,
                            var_interest = "p2",
                            effect = "decrease",
                            interval = "confidence",
                            coefficients = mod$coefficients,
                            vcov = vcov(mod)))
#> ✔ Finished data preparation.
#>          p1    p2        p3        p4 .Sp .Proportion .Group    .Pred   .Lower
#> 1 0.7777778 0.000 0.1111111 0.1111111  p2       0.000      1 10.41016 9.705914
#> 2 0.7770000 0.001 0.1110000 0.1110000  p2       0.001      1 10.41029 9.706958
#> 3 0.7762222 0.002 0.1108889 0.1108889  p2       0.002      1 10.41041 9.708001
#> 4 0.7754444 0.003 0.1107778 0.1107778  p2       0.003      1 10.41054 9.709043
#> 5 0.7746667 0.004 0.1106667 0.1106667  p2       0.004      1 10.41066 9.710084
#> 6 0.7738889 0.005 0.1105556 0.1105556  p2       0.005      1 10.41079 9.711123
#>     .Upper .Marginal .Threshold .MarEffect  .Effect
#> 1 11.11441 0.1252706      0.067   Negative decrease
#> 2 11.11361 0.1252706      0.067   Negative decrease
#> 3 11.11282 0.1252706      0.067   Negative decrease
#> 4 11.11203 0.1252706      0.067   Negative decrease
#> 5 11.11124 0.1252706      0.067   Negative decrease
#> 6 11.11045 0.1252706      0.067   Negative decrease

## Can also create only the intermediary communities without predictions
## by specifying prediction = FALSE.
## Any additional columns can then be added and the `add_prediction` function
## can be manually called.
## Note: If calling the `add_prediction` function manually, the data would
## not contain information about the marginal effect of changing the species
## interest
effects_data <- visualise_effects_data(data = init_comms, prop = 1:4,
                                       var_interest = "p2",
                                       effect = "decrease",
                                       prediction = FALSE)
#> ✔ Finished data preparation.
head(effects_data)
#>          p1    p2        p3        p4 .Sp .Proportion .Group  .Effect
#> 1 0.7777778 0.000 0.1111111 0.1111111  p2       0.000      1 decrease
#> 2 0.7770000 0.001 0.1110000 0.1110000  p2       0.001      1 decrease
#> 3 0.7762222 0.002 0.1108889 0.1108889  p2       0.002      1 decrease
#> 4 0.7754444 0.003 0.1107778 0.1107778  p2       0.003      1 decrease
#> 5 0.7746667 0.004 0.1106667 0.1106667  p2       0.004      1 decrease
#> 6 0.7738889 0.005 0.1105556 0.1105556  p2       0.005      1 decrease
## Prediction using model object
head(add_prediction(data = effects_data, model = mod, interval = "prediction"))
#>          p1    p2        p3        p4 .Sp .Proportion .Group  .Effect    .Pred
#> 1 0.7777778 0.000 0.1111111 0.1111111  p2       0.000      1 decrease 10.41016
#> 2 0.7770000 0.001 0.1110000 0.1110000  p2       0.001      1 decrease 10.41029
#> 3 0.7762222 0.002 0.1108889 0.1108889  p2       0.002      1 decrease 10.41041
#> 4 0.7754444 0.003 0.1107778 0.1107778  p2       0.003      1 decrease 10.41054
#> 5 0.7746667 0.004 0.1106667 0.1106667  p2       0.004      1 decrease 10.41066
#> 6 0.7738889 0.005 0.1105556 0.1105556  p2       0.005      1 decrease 10.41079
#>     .Lower   .Upper
#> 1 7.738930 13.08139
#> 2 7.739308 13.08126
#> 3 7.739686 13.08114
#> 4 7.740062 13.08101
#> 5 7.740439 13.08088
#> 6 7.740814 13.08076
## Prediction using regression coefficients
head(add_prediction(data = effects_data, coefficients = mod$coefficients))
#>          p1    p2        p3        p4 .Sp .Proportion .Group  .Effect    .Pred
#> 1 0.7777778 0.000 0.1111111 0.1111111  p2       0.000      1 decrease 10.41016
#> 2 0.7770000 0.001 0.1110000 0.1110000  p2       0.001      1 decrease 10.41029
#> 3 0.7762222 0.002 0.1108889 0.1108889  p2       0.002      1 decrease 10.41041
#> 4 0.7754444 0.003 0.1107778 0.1107778  p2       0.003      1 decrease 10.41054
#> 5 0.7746667 0.004 0.1106667 0.1106667  p2       0.004      1 decrease 10.41066
#> 6 0.7738889 0.005 0.1105556 0.1105556  p2       0.005      1 decrease 10.41079
```
