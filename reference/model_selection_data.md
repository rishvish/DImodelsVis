# Prepare data for visualising model selection

The data preparation function for visualising model selection. The
output of this function can be passed to the
[`model_selection_plot`](https://rishvish.github.io/DImodelsVis/reference/model_selection_plot.md)
function for showing a visual comparison between the information
criteria for different models. It is also possible to visualise a
breakup of the information criteria into deviance (goodness-of-fit) and
penalty terms for each model.

## Usage

``` r
model_selection_data(
  models,
  metric = c("AIC", "BIC", "AICc", "BICc", "deviance"),
  sort = FALSE,
  breakup = FALSE,
  model_names = names(models)
)
```

## Arguments

- models:

  List of statistical regression model objects.

- metric:

  Metric used for comparisons between models. Takes values from c("AIC",
  "BIC", "AICc", "BICc", "logLik"). Can choose a single or multiple
  metrics for comparing the different models.

- sort:

  A boolean value indicating whether to sort the model from highest to
  lowest value of chosen metric.

- breakup:

  A boolean value indicating whether to breakup the metric value into
  deviance (defined as -2\*loglikelihood) and penalty components. Will
  work only if a single metric out of "AIC", "AICc", "BIC", or "BICc" is
  chosen to plot.

- model_names:

  A character string describing the names to display on X-axis for each
  model in order they appear in the models parameter.

## Value

A data-frame with multiple columns containing values of several
information criteria for each model specified in \`models\`.

- model_name:

  An identifier name for each model object to be shown on X-axis.

- deviance:

  The deviance values for each model object.

- logLik:

  The -2\*Log-Likelihood values for each model object.

- AIC:

  The Akaike information criteria (AIC) values for each model object.

- BIC:

  The Bayesian information criteria (BIC) values for each model object.

- AICc:

  The corrected AIC (AICc) values for each model object.

- BICc:

  The corrected BIC (BICc) values for each model object.

- Component:

  The names of the components to be shown in the plot.

- Value:

  The values for the components to be shown in the plot.

## Examples

``` r
## Fit different candidate models
mod1 <- lm(mpg ~ disp, data = mtcars)
mod2 <- lm(mpg ~ disp + hp, data = mtcars)
mod3 <- lm(mpg ~ disp + hp + wt, data = mtcars)
mod4 <- lm(mpg ~ disp + hp + wt + carb, data = mtcars)

## Group models into list
models_list <- list("Model 1" = mod1, "Model 2" = mod2,
                    "Model 3" = mod3, "Model 4" = mod4)

## Prepare data for visualisation
## Specific metric
model_selection_data(models = models_list,
                     metric = c("AIC"))
#> # A tibble: 4 × 9
#>   model_name deviance logLik   AIC   BIC  AICc  BICc Component Value
#>   <fct>         <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     <dbl>
#> 1 Model 1        317.   164.  170.  175.  171.  176. AIC        170.
#> 2 Model 2        283.   161.  169.  174.  170.  177. AIC        169.
#> 3 Model 3        195.   149.  159.  166.  161.  170. AIC        159.
#> 4 Model 4        194.   148.  161.  169.  164.  175. AIC        161.

## Multiple metrics can be plotted together as well
model_selection_data(models = models_list,
                     metric = c("AIC", "BIC"))
#> # A tibble: 8 × 9
#>   model_name deviance logLik   AIC   BIC  AICc  BICc Component Value
#>   <fct>         <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     <dbl>
#> 1 Model 1        317.   164.  170.  175.  171.  176. AIC        170.
#> 2 Model 1        317.   164.  170.  175.  171.  176. BIC        175.
#> 3 Model 2        283.   161.  169.  174.  170.  177. AIC        169.
#> 4 Model 2        283.   161.  169.  174.  170.  177. BIC        174.
#> 5 Model 3        195.   149.  159.  166.  161.  170. AIC        159.
#> 6 Model 3        195.   149.  159.  166.  161.  170. BIC        166.
#> 7 Model 4        194.   148.  161.  169.  164.  175. AIC        161.
#> 8 Model 4        194.   148.  161.  169.  164.  175. BIC        169.

## If single metric is specified then breakup of metric
## between goodness of fit and penalty can also be visualised
model_selection_data(models = models_list,
                     metric = c("AICc"),
                     breakup = TRUE)
#> # A tibble: 8 × 9
#>   model_name deviance logLik   AIC   BIC  AICc  BICc Component        Value
#>   <fct>         <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <chr>            <dbl>
#> 1 Model 1        317.   164.  170.  175.  171.  176. Goodness of fit 164.  
#> 2 Model 1        317.   164.  170.  175.  171.  176. Penalty           6.87
#> 3 Model 2        283.   161.  169.  174.  170.  177. Goodness of fit 161.  
#> 4 Model 2        283.   161.  169.  174.  170.  177. Penalty           9.48
#> 5 Model 3        195.   149.  159.  166.  161.  170. Goodness of fit 149.  
#> 6 Model 3        195.   149.  159.  166.  161.  170. Penalty          12.3 
#> 7 Model 4        194.   148.  161.  169.  164.  175. Goodness of fit 148.  
#> 8 Model 4        194.   148.  161.  169.  164.  175. Penalty          15.4 

## Sort models
model_selection_data(models = models_list,
                     metric = c("AICc"),
                     breakup = TRUE, sort = TRUE)
#> # A tibble: 8 × 9
#>   model_name deviance logLik   AIC   BIC  AICc  BICc Component        Value
#>   <fct>         <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <chr>            <dbl>
#> 1 Model 3        195.   149.  159.  166.  161.  170. Goodness of fit 149.  
#> 2 Model 3        195.   149.  159.  166.  161.  170. Penalty          12.3 
#> 3 Model 4        194.   148.  161.  169.  164.  175. Goodness of fit 148.  
#> 4 Model 4        194.   148.  161.  169.  164.  175. Penalty          15.4 
#> 5 Model 2        283.   161.  169.  174.  170.  177. Goodness of fit 161.  
#> 6 Model 2        283.   161.  169.  174.  170.  177. Penalty           9.48
#> 7 Model 1        317.   164.  170.  175.  171.  176. Goodness of fit 164.  
#> 8 Model 1        317.   164.  170.  175.  171.  176. Penalty           6.87

## If multiple metrics are specified then sorting
## will be done on first metric specified in list (AIC in this case)
model_selection_data(models = models_list,
                     metric = c("AIC", "BIC", "AICc", "BICc"), sort = TRUE)
#> # A tibble: 16 × 9
#>    model_name deviance logLik   AIC   BIC  AICc  BICc Component Value
#>    <fct>         <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     <dbl>
#>  1 Model 3        195.   149.  159.  166.  161.  170. AIC        159.
#>  2 Model 3        195.   149.  159.  166.  161.  170. BIC        166.
#>  3 Model 3        195.   149.  159.  166.  161.  170. AICc       161.
#>  4 Model 3        195.   149.  159.  166.  161.  170. BICc       170.
#>  5 Model 4        194.   148.  161.  169.  164.  175. AIC        161.
#>  6 Model 4        194.   148.  161.  169.  164.  175. BIC        169.
#>  7 Model 4        194.   148.  161.  169.  164.  175. AICc       164.
#>  8 Model 4        194.   148.  161.  169.  164.  175. BICc       175.
#>  9 Model 2        283.   161.  169.  174.  170.  177. AIC        169.
#> 10 Model 2        283.   161.  169.  174.  170.  177. BIC        174.
#> 11 Model 2        283.   161.  169.  174.  170.  177. AICc       170.
#> 12 Model 2        283.   161.  169.  174.  170.  177. BICc       177.
#> 13 Model 1        317.   164.  170.  175.  171.  176. AIC        170.
#> 14 Model 1        317.   164.  170.  175.  171.  176. BIC        175.
#> 15 Model 1        317.   164.  170.  175.  171.  176. AICc       171.
#> 16 Model 1        317.   164.  170.  175.  171.  176. BICc       176.
```
