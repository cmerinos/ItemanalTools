# ItemanalTools 1.0.0
 * 1sr summision: CRAN
 * Edition: format

# ItemanalTools 0.1.5
 * Edition: format of roxygen documentation for several version
 
## Internal changes
 * `itemThresholds()`: the standard errors of the differences between
  consecutive thresholds (`diff_tests$se.diff`) are now computed using the
  full covariance matrix of the thresholds, i.e.
  `SE(tau_{k+1} - tau_k) = sqrt(SE_k^2 + SE_{k+1}^2 - 2*Cov(tau_k, tau_{k+1}))`.
  Previously, the covariance term was ignored. As a result, the values of
  `z`, `p.value`, and any downstream statistics that depended on them have
  changed. The new computation is the statistically correct one
  
  * `itemThresholds()` gains a monotonicity diagnostic. For each adjacent pair
  of thresholds, the function now returns:
  - the raw difference `tau_{k+1} - tau_k`,
  - a z statistic based on the corrected standard error,
  - a two-tailed p-value (`p.value`),
  - a one-tailed p-value for the monotonicity test
    (H1: `tau_{k+1} - tau_k < 0`, `p.one.sided`),
  - a binary flag `disordered` (1 if the raw difference is negative).
  A new element `monotonicity` collects this diagnostic in a dedicated
  table, and a new element `any_disord` flags the item as a whole.

 * `itemThresholds()` now reports a standardized effect size (Cohen's d) for
  the difference between consecutive thresholds. The conversion depends on
  the `link` argument:
  - with `link = "probit"`, the raw differences are already on a scale with
    SD = 1, so `d = tau_{k+1} - tau_k`;
  - with `link = "logit"`, raw differences are divided by `pi/sqrt(3)`;
  - with other links, `d.cohen` is `NA`.
  A new column `redundant` flags pairs with `d < 0.20` as functionally
  redundant categories. See the `@details` section of `itemThresholds()`
  for interpretation benchmarks.

* `itemThresholdsMult()` gains a `scale_summary` component: a one-row
  data frame with aggregate, descriptive diagnostics of the analyzed items.
  It includes counts (`n_items`, `n_disord_items`, `n_redundant_pairs`),
  summary statistics of the Equidistance Index (`EI.mean`, `EI.median`,
  `EI.min`, `EI.max`, `EI.sd`), summary statistics of Cohen's d
  (`d.mean`, `d.median`, `d.min`, `d.max`), and a global flag
  `scale_ok`. The aggregate view is descriptive and does not assume a
  latent dimension underlying the items.

* `itemThresholdsMult()` per-item `summary` now includes two new columns:
  `any_disord` (from `itemThresholds()`) and `n_redundant`
  (number of adjacent pairs flagged as redundant in that item).

* `itemThresholdsMult()` formatted tables per item now include the columns
  `p.diff`, `p.disord`, `d.cohen`, and `disordered`, mirroring the new
  outputs of `itemThresholds()`.
  
 
# ItemanalTools 0.1.4
 - Remove: RotheryItems
 - Check equations: MurakiPairs
 
# ItemanaTools 0.1.3
 - Package, project and web: changing name
 
# Itemanalysis 0.1.2
 - Non ASCII characters: removing

# Itemanalysis 0.1.1
 - Changes pre CRAN
 - Creation of development brach
 
# Itemanalysis 0.1.0

## Initial steps: creation of functions

* Provides a collection of tools for item-level statistical analysis in psychometrics, including:
  - Effective Number of Nominal Options (AENO)
  - Item-level association coefficients with external variables
  - Equidistance in the thresholds of Likert-type item categories
  - Monotonicity evaluation for Likert-type item categories
  - Between-group comparisons of central tendency and dispersion
  - Deviation of observed correlations from reference values
  - Association tests between items and categorical grouping variables
  - Concordance and difference tests across multiple items

* Includes visualization tools for item-level indices.

* Designed for researchers and practitioners in psychological, educational, and allied field measurement, with a focus on scale development, and content validity.
