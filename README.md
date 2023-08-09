
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MultiGroupSequential

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/MultiGroupSequential)](https://CRAN.R-project.org/package=MultiGroupSequential)
<!-- badges: end -->

It is often challenging to strongly control the family-wise type-1 error
rate in the group-sequential trials with multiple endpoints
(hypotheses). The inflation of type-1 error rate comes from two sources
(1) repeated testing individual hypothesis and (2) simultaneous testing
multiple hypotheses. The **MultiGroupSequential** package is intended to
help researchers to tackle this challenge.

The procedures provided include the sequential procedures described in
[Luo and Quan (2023)](https://doi.org/10.1080/19466315.2023.2191989) and
the graphical procedure proposed by [Maurer and Bretz
(2013)](https://doi.org/10.1080/19466315.2013.807748). Luo and Quan
(2023) describes three procedures and functions to implement these
procedures

1.  `seqgspgx()` implements a sequential graphical procedure based on
    the group-sequential p-values.
2.  `seqgsphh()` implements a sequential Hochberg/Hommel procedure based
    on the group-sequential p-values.
3.  `seqqvalhh()` implements a sequential Hochberg/Hommel procedure
    based on the q-values.

In addition, `seqmbgx()` implements the sequential graphical procedure
described in Maurer and Bretz (2013).

## Installation

You can install the released version of MultiGroupSequential from CRAN:

``` r
install.packages("MultiGroupSequential")
```

## Usage

For example, to use the sequential graphical procedure based on group
sequential p-values.

- The input matrix `pm` has
  - Rows for different hypotheses, and
  - Columns for the group sequential p-values at different times.
- `alpha` is the overall family-wise type-1 error rate.
- `W` is the weights of the graph assigned to each hypothesis, whereas
  `G` holds the transition matrix of the graph.

The procedures implemented here will usually give output list with
elements:

- `rejected`: the index set of rejected hypotheses
- `decisionsm`: rejection decision for each hypothesis (row) at each
  time point (column)
- `cumdecisionsm`: cumulative rejection decision for each hypothesis
  (row) at each time point (column)

``` r
library(MultiGroupSequential)
seqgspgx(
  pm = rbind(c(0.02, 0.03, 0.01), c(0.03, 0.04, 0.01)),
  alpha = 0.025,
  W = c(0.6, 0.4),
  G = rbind(c(0, 1), c(1, 0))
)
#> $rejected
#> [1] 1 2
#> 
#> $decisionsm
#>      [,1] [,2] [,3]
#> [1,]    0    0    1
#> [2,]    0    0    1
#> 
#> $cumdecisionsm
#>      [,1] [,2] [,3]
#> [1,]    0    0    1
#> [2,]    0    0    1
```
