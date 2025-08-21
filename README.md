# Itemanalysis

[![R-CMD-check](https://github.com/tuusuario/Itemanalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tuusuario/Itemanalysis/actions)

## Overview

**Itemanalysis** is an R package for in-depth item-level psychometric analysis. It provides functions to assess:

-   **Effective Number of Nominal Options (Samejima's AENO)**
-   **Item-level associations with external variables**
-   **Monotonic trends in response categories of Likert-type items**
-   **Between-group differences in both central tendency and dispersion**
-   **Deviation of observed correlations from a reference value**
-   **Associations between items and categorical grouping variables**
-   **Summary of item-level association coefficients**
-   **Concordance and differences tests for multiple items**
-   **Quantifying item validity**

The package is aimed at researchers and practitioners conducting scale development, validation, and item refinement using both classical and modern psychometric methods. It is especially recommended for obtaining evidence on the statistical behavior of items within a content validity framework, with a focus on quantitative indicators.

## Installation

To install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("cmerinos/Itemanalysis")
```

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
|----------------------|--------------------------------------------------|
| `AENO()` | Calculates the Factorial Simplicity Index (FSI) |
| `AENO()` | Computes Hofmann’s index of factor complexity |
| `cat.order()` | Monotonicity test using category means across response levels |
| `catPlot()` | Plot of response category trend |
| `RGitems()` | Nonparametric effect size Wilcoxon Rank-Biserial in two grups |
| `CucconiMult()` | Multigroup Cucconi Test for Location and Scale Differences |
| `DmIndex()` | Item-Level Deviation Index from a Theoretical Association Value |
| `epsilonItems()` | Nonparametric effect size epsilon squared across multiple groups |
| `FWitems()` | Friedman test and Kendall’s W for repeated-measures item data |
| `MurakamiPairs()` | All-pairs nonparametric comparisons between groups using the two-groups Cucconi test |
| `RotheryItems()` | Nonparametric concordance index for repeated ordinal measurements, based in Rothery (1979) |
| `Spearitems()` | Spearman's rho or Kendall's tau correlation between each ordinal item and a continuous external criterion |
| `summaryItemAssoc()` | Summary Statistics for Item-Level Association Coefficients |

## References

-   Merino-Soto, C., Juárez-García, A., Salinas-Escudero, G., & Toledano-Toledano, F. (2022). Item-level psychometric analysis of the psychosocial processes at Work Scale (PROPSIT) in workers. *International Journal of Environmental Research and Public Health*, 19(13), 7972. <https://doi.org/10.3390/ijerph19137972>
-   Rosario-Hernández, E., Rovira-Millán, L. V., Merino-Soto, C., & Angulo-Ramos, M. (2023). Review of the psychometric properties of the Patient Health Questionnaire-9 (PHQ-9) Spanish version in a sample of Puerto Rican workers. *Frontiers in psychiatry*, 14, 1024676. <https://doi.org/10.3389/fpsyt.2023.1024676>
-   Jin, Q., & Yu, K. (2023). Adaptation and validation of the university-to-work success scale among Chinese university graduates. *Frontiers in psychology*, 14, 1258746. <https://doi.org/10.3389/fpsyg.2023.1258746>
-   Cabedo-Peris, J., Merino-Soto, C., Chans, G. M., & Martí-Vilar, M. (2024). Exploring the Loss Aversion Scale's psychometric properties in Spain. *Scientific reports*, 14(1), 15756. <https://doi.org/10.1038/s41598-024-66695-6>
-   Merino-Soto, C., Angulo-Ramos, M., Rovira-Millán, L. V., & Rosario-Hernández, E. (2023). Psychometric properties of the generalized anxiety disorder-7 (GAD-7) in a sample of workers. *Frontiers in psychiatry*, 14, 999242. <https://doi.org/10.3389/fpsyt.2023.999242>

## Contributing

Feel free to report issues or contribute improvements via pull requests. You can also suggest new features by opening an issue.

## License

This package is licensed under the **GPL-3** license.
