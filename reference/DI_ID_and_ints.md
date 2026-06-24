# Add identity and interaction terms used in a Diversity-Interactions (DI) model to new data

Utility functions that accept a fitted Diversity-Interactions (DI) model
object along with a data frame and adds the necessary identity and
interaction structures to the data for making predictions using the
model object specified in `` `model` ``.

## Usage

``` r
add_ID_terms(data, model)

add_interaction_terms(data, model)
```

## Arguments

- data:

  A data-frame with species proportions that sum to 1 to create the
  identity effect groupings.

- model:

  A Diversity Interactions model object fit using the
  [`DI()`](https://rdrr.io/pkg/DImodels/man/DI.html) or
  [`autoDI()`](https://rdrr.io/pkg/DImodels/man/autoDI.html) functions
  from the
  [`DImodels`](https://rdrr.io/pkg/DImodels/man/DImodels-package.html)
  or
  [`DImulti()`](https://di-laurabyrne.github.io/DImodelsMulti/reference/DImulti.html)
  from the
  [`DImodelsMulti`](https://di-laurabyrne.github.io/DImodelsMulti/reference/DImodelsMulti.html)
  R packages.

## Value

The original data-frame with additional columns appended at the end
describing the identity and interactions terms present in the model
object.

## Examples

``` r
library(DImodels)
data(sim1)

# Fit DI models with different ID effect groupings
mod1 <- DI(y = "response", prop = 3:6,
           data = sim1, DImodel = "AV") # No ID grouping
#> Fitted model: Average interactions 'AV' DImodel
mod2 <- DI(y = "response", prop = 3:6,
           data = sim1, DImodel = "AV",
           ID = c("ID1", "ID1", "ID2", "ID2"))
#> Fitted model: Average interactions 'AV' DImodel
mod3 <- DI(y = "response", prop = 3:6,
           data = sim1, DImodel = "AV",
           ID = c("ID1", "ID1", "ID1", "ID1"))
#> Fitted model: Average interactions 'AV' DImodel

# Create new data for adding interaction terms
newdata <- sim1[sim1$block == 1, 3:6]
print(head(newdata))
#>      p1   p2   p3   p4
#> 1  0.70 0.10 0.10 0.10
#> 5  0.10 0.70 0.10 0.10
#> 9  0.10 0.10 0.70 0.10
#> 13 0.10 0.10 0.10 0.70
#> 17 0.25 0.25 0.25 0.25
#> 21 0.40 0.40 0.10 0.10

add_ID_terms(data = newdata, model = mod1)
#>      p1   p2   p3   p4 p1_ID p2_ID p3_ID p4_ID
#> 1  0.70 0.10 0.10 0.10  0.70  0.10  0.10  0.10
#> 5  0.10 0.70 0.10 0.10  0.10  0.70  0.10  0.10
#> 9  0.10 0.10 0.70 0.10  0.10  0.10  0.70  0.10
#> 13 0.10 0.10 0.10 0.70  0.10  0.10  0.10  0.70
#> 17 0.25 0.25 0.25 0.25  0.25  0.25  0.25  0.25
#> 21 0.40 0.40 0.10 0.10  0.40  0.40  0.10  0.10
#> 25 0.40 0.10 0.40 0.10  0.40  0.10  0.40  0.10
#> 29 0.40 0.10 0.10 0.40  0.40  0.10  0.10  0.40
#> 33 0.10 0.40 0.40 0.10  0.10  0.40  0.40  0.10
#> 37 0.10 0.40 0.10 0.40  0.10  0.40  0.10  0.40
#> 41 0.10 0.10 0.40 0.40  0.10  0.10  0.40  0.40
#> 45 1.00 0.00 0.00 0.00  1.00  0.00  0.00  0.00
#> 49 0.00 1.00 0.00 0.00  0.00  1.00  0.00  0.00
#> 53 0.00 0.00 1.00 0.00  0.00  0.00  1.00  0.00
#> 57 0.00 0.00 0.00 1.00  0.00  0.00  0.00  1.00
add_ID_terms(data = newdata, model = mod2)
#>      p1   p2   p3   p4 ID1 ID2
#> 1  0.70 0.10 0.10 0.10 0.8 0.2
#> 5  0.10 0.70 0.10 0.10 0.8 0.2
#> 9  0.10 0.10 0.70 0.10 0.2 0.8
#> 13 0.10 0.10 0.10 0.70 0.2 0.8
#> 17 0.25 0.25 0.25 0.25 0.5 0.5
#> 21 0.40 0.40 0.10 0.10 0.8 0.2
#> 25 0.40 0.10 0.40 0.10 0.5 0.5
#> 29 0.40 0.10 0.10 0.40 0.5 0.5
#> 33 0.10 0.40 0.40 0.10 0.5 0.5
#> 37 0.10 0.40 0.10 0.40 0.5 0.5
#> 41 0.10 0.10 0.40 0.40 0.2 0.8
#> 45 1.00 0.00 0.00 0.00 1.0 0.0
#> 49 0.00 1.00 0.00 0.00 1.0 0.0
#> 53 0.00 0.00 1.00 0.00 0.0 1.0
#> 57 0.00 0.00 0.00 1.00 0.0 1.0
add_ID_terms(data = newdata, model = mod3)
#>      p1   p2   p3   p4 ID1
#> 1  0.70 0.10 0.10 0.10   1
#> 5  0.10 0.70 0.10 0.10   1
#> 9  0.10 0.10 0.70 0.10   1
#> 13 0.10 0.10 0.10 0.70   1
#> 17 0.25 0.25 0.25 0.25   1
#> 21 0.40 0.40 0.10 0.10   1
#> 25 0.40 0.10 0.40 0.10   1
#> 29 0.40 0.10 0.10 0.40   1
#> 33 0.10 0.40 0.40 0.10   1
#> 37 0.10 0.40 0.10 0.40   1
#> 41 0.10 0.10 0.40 0.40   1
#> 45 1.00 0.00 0.00 0.00   1
#> 49 0.00 1.00 0.00 0.00   1
#> 53 0.00 0.00 1.00 0.00   1
#> 57 0.00 0.00 0.00 1.00   1
library(DImodels)
data(sim1)

# Fit different DI models
mod1 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "AV")
#> Fitted model: Average interactions 'AV' DImodel
mod2 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "FULL")
#> Fitted model: Separate pairwise interactions 'FULL' DImodel
mod3 <- DI(y = "response", prop = 3:6, data = sim1, DImodel = "ADD")
#> Fitted model: Additive species contributions to interactions 'ADD' DImodel
mod4 <- DI(y = "response", prop = 3:6, data = sim1,
           FG = c("G", "G", "H", "H"), DImodel = "FG")
#> Fitted model: Functional group effects 'FG' DImodel

# Create new data for adding interaction terms
newdata <- sim1[sim1$block == 1, 3:6]
print(head(newdata))
#>      p1   p2   p3   p4
#> 1  0.70 0.10 0.10 0.10
#> 5  0.10 0.70 0.10 0.10
#> 9  0.10 0.10 0.70 0.10
#> 13 0.10 0.10 0.10 0.70
#> 17 0.25 0.25 0.25 0.25
#> 21 0.40 0.40 0.10 0.10

add_interaction_terms(data = newdata, model = mod1)
#>      p1   p2   p3   p4    AV
#> 1  0.70 0.10 0.10 0.10 0.240
#> 5  0.10 0.70 0.10 0.10 0.240
#> 9  0.10 0.10 0.70 0.10 0.240
#> 13 0.10 0.10 0.10 0.70 0.240
#> 17 0.25 0.25 0.25 0.25 0.375
#> 21 0.40 0.40 0.10 0.10 0.330
#> 25 0.40 0.10 0.40 0.10 0.330
#> 29 0.40 0.10 0.10 0.40 0.330
#> 33 0.10 0.40 0.40 0.10 0.330
#> 37 0.10 0.40 0.10 0.40 0.330
#> 41 0.10 0.10 0.40 0.40 0.330
#> 45 1.00 0.00 0.00 0.00 0.000
#> 49 0.00 1.00 0.00 0.00 0.000
#> 53 0.00 0.00 1.00 0.00 0.000
#> 57 0.00 0.00 0.00 1.00 0.000
add_interaction_terms(data = newdata, model = mod2)
#>      p1   p2   p3   p4  p1:p2  p1:p3  p1:p4  p2:p3  p2:p4  p3:p4
#> 1  0.70 0.10 0.10 0.10 0.0700 0.0700 0.0700 0.0100 0.0100 0.0100
#> 5  0.10 0.70 0.10 0.10 0.0700 0.0100 0.0100 0.0700 0.0700 0.0100
#> 9  0.10 0.10 0.70 0.10 0.0100 0.0700 0.0100 0.0700 0.0100 0.0700
#> 13 0.10 0.10 0.10 0.70 0.0100 0.0100 0.0700 0.0100 0.0700 0.0700
#> 17 0.25 0.25 0.25 0.25 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625
#> 21 0.40 0.40 0.10 0.10 0.1600 0.0400 0.0400 0.0400 0.0400 0.0100
#> 25 0.40 0.10 0.40 0.10 0.0400 0.1600 0.0400 0.0400 0.0100 0.0400
#> 29 0.40 0.10 0.10 0.40 0.0400 0.0400 0.1600 0.0100 0.0400 0.0400
#> 33 0.10 0.40 0.40 0.10 0.0400 0.0400 0.0100 0.1600 0.0400 0.0400
#> 37 0.10 0.40 0.10 0.40 0.0400 0.0100 0.0400 0.0400 0.1600 0.0400
#> 41 0.10 0.10 0.40 0.40 0.0100 0.0400 0.0400 0.0400 0.0400 0.1600
#> 45 1.00 0.00 0.00 0.00 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
#> 49 0.00 1.00 0.00 0.00 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
#> 53 0.00 0.00 1.00 0.00 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
#> 57 0.00 0.00 0.00 1.00 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
add_interaction_terms(data = newdata, model = mod3)
#>      p1   p2   p3   p4 p1_add p2_add p3_add p4_add
#> 1  0.70 0.10 0.10 0.10 0.2100 0.0900 0.0900 0.0900
#> 5  0.10 0.70 0.10 0.10 0.0900 0.2100 0.0900 0.0900
#> 9  0.10 0.10 0.70 0.10 0.0900 0.0900 0.2100 0.0900
#> 13 0.10 0.10 0.10 0.70 0.0900 0.0900 0.0900 0.2100
#> 17 0.25 0.25 0.25 0.25 0.1875 0.1875 0.1875 0.1875
#> 21 0.40 0.40 0.10 0.10 0.2400 0.2400 0.0900 0.0900
#> 25 0.40 0.10 0.40 0.10 0.2400 0.0900 0.2400 0.0900
#> 29 0.40 0.10 0.10 0.40 0.2400 0.0900 0.0900 0.2400
#> 33 0.10 0.40 0.40 0.10 0.0900 0.2400 0.2400 0.0900
#> 37 0.10 0.40 0.10 0.40 0.0900 0.2400 0.0900 0.2400
#> 41 0.10 0.10 0.40 0.40 0.0900 0.0900 0.2400 0.2400
#> 45 1.00 0.00 0.00 0.00 0.0000 0.0000 0.0000 0.0000
#> 49 0.00 1.00 0.00 0.00 0.0000 0.0000 0.0000 0.0000
#> 53 0.00 0.00 1.00 0.00 0.0000 0.0000 0.0000 0.0000
#> 57 0.00 0.00 0.00 1.00 0.0000 0.0000 0.0000 0.0000
add_interaction_terms(data = newdata, model = mod4)
#>      p1   p2   p3   p4 FG_.bfg_G_H FG_.wfg_G FG_.wfg_H
#> 1  0.70 0.10 0.10 0.10      0.1600    0.0700    0.0100
#> 5  0.10 0.70 0.10 0.10      0.1600    0.0700    0.0100
#> 9  0.10 0.10 0.70 0.10      0.1600    0.0100    0.0700
#> 13 0.10 0.10 0.10 0.70      0.1600    0.0100    0.0700
#> 17 0.25 0.25 0.25 0.25      0.2500    0.0625    0.0625
#> 21 0.40 0.40 0.10 0.10      0.1600    0.1600    0.0100
#> 25 0.40 0.10 0.40 0.10      0.2500    0.0400    0.0400
#> 29 0.40 0.10 0.10 0.40      0.2500    0.0400    0.0400
#> 33 0.10 0.40 0.40 0.10      0.2500    0.0400    0.0400
#> 37 0.10 0.40 0.10 0.40      0.2500    0.0400    0.0400
#> 41 0.10 0.10 0.40 0.40      0.1600    0.0100    0.1600
#> 45 1.00 0.00 0.00 0.00      0.0000    0.0000    0.0000
#> 49 0.00 1.00 0.00 0.00      0.0000    0.0000    0.0000
#> 53 0.00 0.00 1.00 0.00      0.0000    0.0000    0.0000
#> 57 0.00 0.00 0.00 1.00      0.0000    0.0000    0.0000
```
