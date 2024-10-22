---
output:
  pdf_document: default
  html_document: default
---
# SoftBart

## SoftBart 1.0.1

- Added vignette to the package explaining generic usage. This can be accessed
  by running `browseVignettes('SoftBart')`.
  

## SoftBart 1.1.1

- Introduced `omega_group` in `Hypers` to store the within-group components probability profile, and `alpha_comp` to control the sparsity across components.

- Introduced new functions `UpdateOmegaGroup()`, `get_var_counts_by_variable()` to facilitate handling the scenario where grouping exists.
- Updated `SampleVar()`, `UpdateS()` and `get_tree_counts()` functions so as to implement the hierarchical sampling by use of Dirichlet prior when grouping exists. 
  Now the output type of `get_tree_counts` is a Rcpp::List rather than a matrix.
- Updated `UpdateAlpha()` function so that when `update_alpha = TRUE`, the hyperparameters in Dirichlet priors are updated with scaled beta prime prior. 

