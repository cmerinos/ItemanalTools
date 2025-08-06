# Itemanalysis

[![R-CMD-check](https://github.com/tuusuario/Itemanalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tuusuario/Itemanalysis/actions)

## Overview

**Itemanalysis** is an R package for in-depth item-level psychometric analysis. It provides functions to assess:

-   **Effective Number of Nominal Options (AENO)**
-   **Item-level associations with external variables**
-   **Monotonic trends in response categories of Likert-type items**
-   **Between-group differences in both central tendency and dispersion**
-   **Deviation of observed correlations from a reference value**
-   **Associations between items and categorical grouping variables**
-   **Summary of item-level association coefficients**


The package is aimed at researchers and practitioners conducting scale development, validation, and item refinement using both classical and modern psychometric methods.

## Installation

To install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("tuusuario/Itemanalysis")
```

Replace `"tuusuario"` with your actual GitHub username.

## Example: FSI (Factorial Simplicity Index)

``` r
library(Itemanalysis)

# Simulated factor loadings (4 items, 2 factors)
loadings <- matrix(c(
  0.8, 0.2,
  0.75, 0.3,
  0.1, 0.9,
  0.2, 0.85
), nrow = 4, byrow = TRUE)

# Calculate FSI
fsi <- FSI(loadings)
print(fsi$FSI.i)
```

## Selected Functions

| Function | Description |
|------------------|------------------------------------------------------|
| `FSI()` | Calculates the Factorial Simplicity Index (FSI) |
| `HofmannFac()` | Computes Hofmann’s index of factor complexity |
| `RGWil()` | Computes rank-biserial correlation by groups (Wilcoxon-based) |
| `SVALsingle()` | Substantive validity (psa, svc) for one item |
| `SVALmult()` | Substantive validity for multiple items |
| `BSIbootpp()` | Bootstrap-based confidence intervals for pre-post item change |
| `AENOmulti()` | Calculates entropy for each ordinal item |
| `FWitems()` | Friedman test and Kendall’s W for repeated-measures item data |
| `cat.order()` | Monotonicity test using category means across response levels |

## References

-   Fleming, J. S., & Merino Soto, C. (2005). *Medidas de simplicidad y de ajuste factorial*. Revista de Psicología, 23(2), 250–266.
-   Fleming, J. S. (2003). *Computing measures of simplicity of fit*. Behavior Research Methods, 35(4), 520–524.
-   Botella, J., Blázquez, D., Suero, M., & Juola, J. F. (2018). *Item Bootstrapping*. Frontiers in Psychology, 9, 223.
-   Hofmann, R. J. (1977). *Indices of factorial complexity*. Behavior Research Methods and Instrumentation, 9, 538–550.

## Contributing

Feel free to report issues or contribute improvements via pull requests. You can also suggest new features by opening an issue.

## License

This package is licensed under the **GPL-3** license.
