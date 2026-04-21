# ModSoftBart

`ModSoftBart` extends the general soft BART regression framework to support explicit grouped predictor structure in models that combine linear covariate effects with a nonparametric tree ensemble function.

## Overview

The original `SoftBart` package includes separate main functions for regression, probit, and general regression models. In the original package, the general regression interface allows linear covariate effects to be modeled jointly with the nonparametric tree-based mean function.

`ModSoftBart` builds on that general modeling framework. Starting from the structure of `gsoftbart_regression()`, the package introduces explicit handling of grouped predictors through hierarchical prior updates and related hyperparameter modifications. It also provides a corresponding binary outcome model through `pb_gsbart()`.

## Main functions

The primary user-facing functions are:

-   `gsbart_regression()` for continuous outcomes
-   `pb_gsbart()` for binary outcomes

Post-fitting utilities include:

-   `calPIPs_gsbart()` for calculating group-level and component-level posterior inclusion probabilities (PIPs) from a fitted `gsbart` model

## Main features

-   SoftBART models with linear covariate effects and a nonparametric tree-based mean function
-   explicit handling of grouped predictors
-   hierarchical updates for grouped predictor selection
-   continuous and binary outcome modeling
-   utilities for grouped predictor summaries and updates

## Relationship to `SoftBart`

`ModSoftBart` is a modified and extended version of the general modeling framework in `SoftBart`. It is built from the structure of `gsoftbart_regression()`, with additional methodology for grouped predictor selection and related prior hyperparameters. The package also adds a binary outcome counterpart through `pb_gsbart()`.

The package builds on the original `SoftBart` code base developed by Antonio R. Linero.

## Methodological background

The package is based on the Soft Bayesian Additive Regression Trees framework described in:

-   Linero, A. R. and Yang, Y. (2018). *Bayesian tree ensembles that adapt to smoothness and sparsity*. Journal of the Royal Statistical Society, Series B.

## Installation

Install the development version from GitHub with:

``` r
# install.packages("remotes")
if (!requireNamespace("remotes", quietly = TRUE)) {

  install.packages("remotes")

}
remotes::install_github("Shiny1818/ModSoftBART_against_BKMR", subdir = "ModSoftBART-master")
```
