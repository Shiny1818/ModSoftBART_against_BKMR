---
output:
  pdf_document: default
  html_document: default
---

## ModSoftBart 1.0.1

- Added `omega_group` to `Hypers` to store the within-group component probability profile, and `alpha_comp` to control sparsity across components.

- Added `UpdateOmegaGroup()` and `get_var_counts_by_variable()` to support grouped predictors.

- Updated `SampleVar()`, `UpdateS()`, and `get_tree_counts()` to implement hierarchical sampling under grouped predictors using a Dirichlet prior. When `group` is not `NULL`, `get_tree_counts()` now returns an `Rcpp::List` instead of a matrix.

- Updated `UpdateAlpha()` so that, when `update_alpha = TRUE`, both group and components Dirichlet prior hyperparameters are updated using a scaled beta-prime prior.


