# Copy attributes from one object to another

This function copies over any additional attributes from \`source\` into
\`target\`. Any attributes already present in \`target\` would be left
untouched. This function is useful after manipulating the data from the
`*_data` preparation functions to ensure any attributes necessary for
creating the plot aren't lost.

## Usage

``` r
copy_attributes(target, source)
```

## Arguments

- target:

  The object to which attributes should be added.

- source:

  The object whose attributes to copy.

## Value

The object specified in \`target\` with all additional attributes in
\`source\` object.

## Examples

``` r

## Simple example
a <- data.frame(Var1 = runif(1:10), Var2 = runif(1:10))
b <- data.frame(Var3 = runif(1:10), Var4 = runif(1:10))
attr(b, "attr1") <- "Lorem"
attr(b, "attr2") <- "ipsum"

print(attributes(a))
#> $names
#> [1] "Var1" "Var2"
#> 
#> $class
#> [1] "data.frame"
#> 
#> $row.names
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
print(attributes(b))
#> $names
#> [1] "Var3" "Var4"
#> 
#> $class
#> [1] "data.frame"
#> 
#> $row.names
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $attr1
#> [1] "Lorem"
#> 
#> $attr2
#> [1] "ipsum"
#> 

## Copy over attributes of `b` into `a`
print(copy_attributes(target = a, source = b))
#>         Var1       Var2
#> 1  0.9907123 0.71339728
#> 2  0.9705209 0.06521611
#> 3  0.3891828 0.35420680
#> 4  0.4611865 0.82519942
#> 5  0.3152418 0.27381825
#> 6  0.1746759 0.57004495
#> 7  0.5315735 0.33571908
#> 8  0.4936370 0.59626279
#> 9  0.7793086 0.19151803
#> 10 0.2041783 0.94776394
## Note the attributes already present in `a` are left untouched

## Can also be used in the dplyr pipeline
library(dplyr)

iris_sub <- iris[1:10, ]
attr(iris_sub, "attr1") <- "Lorem"
attr(iris_sub, "attr2") <- "ipsum"
attributes(iris_sub)
#> $names
#> [1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"     
#> 
#> $row.names
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $class
#> [1] "data.frame"
#> 
#> $attr1
#> [1] "Lorem"
#> 
#> $attr2
#> [1] "ipsum"
#> 

## Grouping can drop attributes we set
iris_sub %>%
   group_by(Species) %>%
   summarise(mean(Sepal.Length)) %>%
   attributes()
#> $names
#> [1] "Species"            "mean(Sepal.Length)"
#> 
#> $row.names
#> [1] 1
#> 
#> $class
#> [1] "tbl_df"     "tbl"        "data.frame"
#> 

## Use copy_attributes with `iris_sub` object as source
##  to add the attributes again
iris_sub %>%
   group_by(Species) %>%
   summarise(mean(Sepal.Length)) %>%
   copy_attributes(source = iris_sub) %>%
   attributes()
#> $names
#> [1] "Species"            "mean(Sepal.Length)"
#> 
#> $row.names
#> [1] 1
#> 
#> $class
#> [1] "tbl_df"     "tbl"        "data.frame"
#> 
#> $attr1
#> [1] "Lorem"
#> 
#> $attr2
#> [1] "ipsum"
#> 
```
