# Add predictions and uncertainty interval to data using either a model object or model coefficients.

This function accepts a data.frame and either a model object or
coefficients and adds columns containing the predictions and associated
uncertainty to the data. When a model object is specified, the function
uses
[`get_predicted()`](https://easystats.github.io/insight/reference/get_predicted.html)
from the
[`insight`](https://easystats.github.io/insight/reference/insight-package.html)
package under the hood to generate the predictions.

## Usage

``` r
add_prediction(
  data,
  model = NULL,
  coefficients = NULL,
  coeff_cols = NULL,
  vcov = NULL,
  interval = c("none", "confidence", "prediction"),
  conf.level = 0.95
)
```

## Arguments

- data:

  A data-frame containing appropriate values for all the terms in the
  model.

- model:

  A regression model object which will be used to make predictions for
  the observations in \`data\`. Will override \`coefficients\` if
  specified.

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

- interval:

  Type of interval to calculate:

  "none" (default)

  :   No interval to be calculated.

  "confidence"

  :   Calculate a confidence interval.

  "prediction"

  :   Calculate a prediction interval.

- conf.level:

  The confidence level for calculating confidence/prediction intervals.
  Default is 0.95.

## Value

A data-frame with the following additional columns

- .Pred:

  The predicted response for each observation.

- .Lower:

  The lower limit of the confidence/prediction interval for each
  observation (will be same as ".Pred" if using \`coefficients\` and
  \`vcov\` is not specified).

- .Upper:

  The lower limit of the confidence/prediction interval for each
  observation (will be same as ".Pred" if using \`coefficients\` and
  \`vcov\` is not specified).

## Examples

``` r
library(DImodels)
data(sim1)

# Fit a model
mod <- lm(response ~ 0 + p1 + p2 + p3 + p4 + p1:p2 + p3:p4, data = sim1)

# Create new data for adding predictions
newdata <- head(sim1[sim1$block == 1,])
print(newdata)
#>    community block   p1   p2   p3   p4 response
#> 1          1     1 0.70 0.10 0.10 0.10   10.815
#> 5          2     1 0.10 0.70 0.10 0.10    6.724
#> 9          3     1 0.10 0.10 0.70 0.10   11.135
#> 13         4     1 0.10 0.10 0.10 0.70    9.916
#> 17         5     1 0.25 0.25 0.25 0.25   11.374
#> 21         6     1 0.40 0.40 0.10 0.10   12.492

# Add predictions to data
add_prediction(data = newdata, model = mod)
#>    community block   p1   p2   p3   p4 response     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 10.601405

# Adding predictions to data with confidence interval
add_prediction(data = newdata, model = mod, interval = "confidence")
#>    community block   p1   p2   p3   p4 response     .Pred   .Lower    .Upper
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 10.466907 9.798148 11.135666
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 10.253858 9.585100 10.922617
#> 9          3     1 0.10 0.10 0.70 0.10   11.135  9.615319 8.946561 10.284078
#> 13         4     1 0.10 0.10 0.10 0.70    9.916  8.584323 7.915564  9.253081
#> 17         5     1 0.25 0.25 0.25 0.25   11.374  9.795128 9.328378 10.261878
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 10.601405 9.531207 11.671603

# Calculate prediction intervals instead
add_prediction(data = newdata, model = mod, interval = "prediction")
#>    community block   p1   p2   p3   p4 response     .Pred   .Lower   .Upper
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 10.466907 7.769709 13.16410
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 10.253858 7.556661 12.95106
#> 9          3     1 0.10 0.10 0.70 0.10   11.135  9.615319 6.918122 12.31252
#> 13         4     1 0.10 0.10 0.10 0.70    9.916  8.584323 5.887125 11.28152
#> 17         5     1 0.25 0.25 0.25 0.25   11.374  9.795128 7.140794 12.44946
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 10.601405 7.777762 13.42505

# Default is a 95% interval, change to 99%
add_prediction(data = newdata, model = mod, interval = "prediction",
               conf.level = 0.99)
#>    community block   p1   p2   p3   p4 response     .Pred   .Lower   .Upper
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 10.466907 6.874932 14.05888
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 10.253858 6.661883 13.84583
#> 9          3     1 0.10 0.10 0.70 0.10   11.135  9.615319 6.023344 13.20729
#> 13         4     1 0.10 0.10 0.10 0.70    9.916  8.584323 4.992347 12.17630
#> 17         5     1 0.25 0.25 0.25 0.25   11.374  9.795128 6.260235 13.33002
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 10.601405 6.841037 14.36177

####################################################################
##### Use model coefficients for prediction
coeffs <- mod$coefficients

# Would now have to add columns corresponding to each coefficient in the
# data and ensure there is an appropriate mapping between data columns and
# the coefficients.
newdata$`p1:p2` = newdata$p1 * newdata$p2
newdata$`p3:p4` = newdata$p3 * newdata$p4

# If the coefficients are named then the function will try to
# perform matching between data columns and the coefficients
# Notice that confidence intervals are not produced if we don't
# specify a variance covariance matrix
add_prediction(data = newdata, coefficients = coeffs)
#>    community block   p1   p2   p3   p4 response  p1:p2  p3:p4     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 0.0700 0.0100 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 0.0700 0.0100 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135 0.0100 0.0700  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916 0.0100 0.0700  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374 0.0625 0.0625  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 0.1600 0.0100 10.601405

# However, if the coefficients are not named
# The user would have to manually specify the subset
# of data columns arranged according to the coefficients
coeffs <- unname(coeffs)

subset_data <- newdata[, c(3:6, 8,9)]
subset_data # Notice now we have the exact columns in data as in coefficients
#>      p1   p2   p3   p4  p1:p2  p3:p4
#> 1  0.70 0.10 0.10 0.10 0.0700 0.0100
#> 5  0.10 0.70 0.10 0.10 0.0700 0.0100
#> 9  0.10 0.10 0.70 0.10 0.0100 0.0700
#> 13 0.10 0.10 0.10 0.70 0.0100 0.0700
#> 17 0.25 0.25 0.25 0.25 0.0625 0.0625
#> 21 0.40 0.40 0.10 0.10 0.1600 0.0100
add_prediction(data = subset_data, coefficients = coeffs)
#>      p1   p2   p3   p4  p1:p2  p3:p4     .Pred
#> 1  0.70 0.10 0.10 0.10 0.0700 0.0100 10.466907
#> 5  0.10 0.70 0.10 0.10 0.0700 0.0100 10.253858
#> 9  0.10 0.10 0.70 0.10 0.0100 0.0700  9.615319
#> 13 0.10 0.10 0.10 0.70 0.0100 0.0700  8.584323
#> 17 0.25 0.25 0.25 0.25 0.0625 0.0625  9.795128
#> 21 0.40 0.40 0.10 0.10 0.1600 0.0100 10.601405

# Or specify a selection (either by name or index) in coeff_cols
add_prediction(data = newdata, coefficients = coeffs,
               coeff_cols = c("p1", "p2", "p3", "p4", "p1:p2", "p3:p4"))
#>    community block   p1   p2   p3   p4 response  p1:p2  p3:p4     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 0.0700 0.0100 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 0.0700 0.0100 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135 0.0100 0.0700  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916 0.0100 0.0700  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374 0.0625 0.0625  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 0.1600 0.0100 10.601405

add_prediction(data = newdata, coefficients = coeffs,
               coeff_cols = c(3, 4, 5, 6, 8, 9))
#>    community block   p1   p2   p3   p4 response  p1:p2  p3:p4     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 0.0700 0.0100 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 0.0700 0.0100 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135 0.0100 0.0700  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916 0.0100 0.0700  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374 0.0625 0.0625  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 0.1600 0.0100 10.601405

# Adding confidence intervals when using model coefficients
coeffs <- mod$coefficients
# We need to provide a variance-covariance matrix to calculate the CI
# when using `coefficients` argument. The following warning will be given
add_prediction(data = newdata, coefficients = coeffs,
               interval = "confidence")
#> Warning: `vcov` was not specified so uncertainty intervals cannot be calculated.
#> ℹ The ".Upper" and ".Lower" columns will contain the same value as the ".Pred"
#>   column.
#>    community block   p1   p2   p3   p4 response  p1:p2  p3:p4     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 0.0700 0.0100 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 0.0700 0.0100 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135 0.0100 0.0700  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916 0.0100 0.0700  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374 0.0625 0.0625  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 0.1600 0.0100 10.601405
#>       .Lower    .Upper
#> 1  10.466907 10.466907
#> 5  10.253858 10.253858
#> 9   9.615319  9.615319
#> 13  8.584323  8.584323
#> 17  9.795128  9.795128
#> 21 10.601405 10.601405

vcov_mat <- vcov(mod)
add_prediction(data = newdata, coefficients = coeffs,
               interval = "confidence", vcov = vcov_mat)
#>    community block   p1   p2   p3   p4 response  p1:p2  p3:p4     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 0.0700 0.0100 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 0.0700 0.0100 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135 0.0100 0.0700  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916 0.0100 0.0700  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374 0.0625 0.0625  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 0.1600 0.0100 10.601405
#>      .Lower    .Upper
#> 1  9.813131 11.120684
#> 5  9.600082 10.907635
#> 9  8.961543 10.269096
#> 13 7.930546  9.238099
#> 17 9.338835 10.251422
#> 21 9.555182 11.647627

# Currently both confidence and prediction intervals will be the same when
# using this method
add_prediction(data = newdata, coefficients = coeffs,
               interval = "prediction", vcov = vcov_mat)
#>    community block   p1   p2   p3   p4 response  p1:p2  p3:p4     .Pred
#> 1          1     1 0.70 0.10 0.10 0.10   10.815 0.0700 0.0100 10.466907
#> 5          2     1 0.10 0.70 0.10 0.10    6.724 0.0700 0.0100 10.253858
#> 9          3     1 0.10 0.10 0.70 0.10   11.135 0.0100 0.0700  9.615319
#> 13         4     1 0.10 0.10 0.10 0.70    9.916 0.0100 0.0700  8.584323
#> 17         5     1 0.25 0.25 0.25 0.25   11.374 0.0625 0.0625  9.795128
#> 21         6     1 0.40 0.40 0.10 0.10   12.492 0.1600 0.0100 10.601405
#>      .Lower    .Upper
#> 1  9.813131 11.120684
#> 5  9.600082 10.907635
#> 9  8.961543 10.269096
#> 13 7.930546  9.238099
#> 17 9.338835 10.251422
#> 21 9.555182 11.647627
```
