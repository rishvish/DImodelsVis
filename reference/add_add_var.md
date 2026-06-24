# Add additional variables to data via cartesian product

Utility function for incorporating any additional variables into new
data via their cartesian product. Each row in the data will be
replicated and new columns will be added for each variable specified in
\`add_var\` with values corresponding to their cartesian product.

## Usage

``` r
add_add_var(data, add_var = NULL)
```

## Arguments

- data:

  A data frame containing the data in which to add the additional
  variables.

- add_var:

  A named list or data-frame specifying the names and corresponding
  values of each new variable to add to the data. If a list is
  specified, each row in the data would be replicated for each unique
  combination of values of the specified variables (i.e., their
  cartesian product) in \`add_var\`, while specifying a data-frame would
  replicate each row in the data for each row in add_var (i.e., merge
  the two data-frames).

## Value

A data-frame with all additional columns specified in \`add_var\` and
the following additional column.

- .add_str_ID:

  A unique identifier describing each element from the cartesian product
  of all variables specified in \`add_var\`.

## Examples

``` r
test_data <- data.frame(diag(1, 3))
print(test_data)
#>   X1 X2 X3
#> 1  1  0  0
#> 2  0  1  0
#> 3  0  0  1

## Adding a single variable
add_add_var(data = test_data,
            add_var = list("Var1" = c(10, 20)))
#>   X1 X2 X3 Var1 .add_str_ID
#> 1  1  0  0   10    Var1: 10
#> 2  0  1  0   10    Var1: 10
#> 3  0  0  1   10    Var1: 10
#> 4  1  0  0   20    Var1: 20
#> 5  0  1  0   20    Var1: 20
#> 6  0  0  1   20    Var1: 20

## Specifying multiple variables as a list will add values for
##  each unique combination
add_add_var(data = test_data,
            add_var = list("Var1" = c(10, 20),
                           "Var2" = c(30, 40)))
#>    X1 X2 X3 Var1 Var2          .add_str_ID
#> 1   1  0  0   10   30 Var1: 10; \tVar2: 30
#> 2   0  1  0   10   30 Var1: 10; \tVar2: 30
#> 3   0  0  1   10   30 Var1: 10; \tVar2: 30
#> 4   1  0  0   20   30 Var1: 20; \tVar2: 30
#> 5   0  1  0   20   30 Var1: 20; \tVar2: 30
#> 6   0  0  1   20   30 Var1: 20; \tVar2: 30
#> 7   1  0  0   10   40 Var1: 10; \tVar2: 40
#> 8   0  1  0   10   40 Var1: 10; \tVar2: 40
#> 9   0  0  1   10   40 Var1: 10; \tVar2: 40
#> 10  1  0  0   20   40 Var1: 20; \tVar2: 40
#> 11  0  1  0   20   40 Var1: 20; \tVar2: 40
#> 12  0  0  1   20   40 Var1: 20; \tVar2: 40

## Specifying add_var as a data.frame would simply merge the two data-frames
add_add_var(data = test_data,
            add_var = data.frame("Var1" = c(10, 20),
                                 "Var2" = c(30, 40)))
#>   X1 X2 X3 Var1 Var2          .add_str_ID
#> 1  1  0  0   10   30 Var1: 10; \tVar2: 30
#> 2  0  1  0   10   30 Var1: 10; \tVar2: 30
#> 3  0  0  1   10   30 Var1: 10; \tVar2: 30
#> 4  1  0  0   20   40 Var1: 20; \tVar2: 40
#> 5  0  1  0   20   40 Var1: 20; \tVar2: 40
#> 6  0  0  1   20   40 Var1: 20; \tVar2: 40

## If the list specified in `add_var` is not named, then the additional
## variables will be automatically named Var1, Var2, Var3, etc.
add_add_var(data = test_data,
            add_var = list(c(1, 2), c(3, 4)))
#>    X1 X2 X3 Var1 Var2        .add_str_ID
#> 1   1  0  0    1    3 Var1: 1; \tVar2: 3
#> 2   0  1  0    1    3 Var1: 1; \tVar2: 3
#> 3   0  0  1    1    3 Var1: 1; \tVar2: 3
#> 4   1  0  0    2    3 Var1: 2; \tVar2: 3
#> 5   0  1  0    2    3 Var1: 2; \tVar2: 3
#> 6   0  0  1    2    3 Var1: 2; \tVar2: 3
#> 7   1  0  0    1    4 Var1: 1; \tVar2: 4
#> 8   0  1  0    1    4 Var1: 1; \tVar2: 4
#> 9   0  0  1    1    4 Var1: 1; \tVar2: 4
#> 10  1  0  0    2    4 Var1: 2; \tVar2: 4
#> 11  0  1  0    2    4 Var1: 2; \tVar2: 4
#> 12  0  0  1    2    4 Var1: 2; \tVar2: 4
```
