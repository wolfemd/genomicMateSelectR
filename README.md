# Genomic Mate Selection in R

<!-- badges: start -->

<!-- badges: end -->

Welcome to the R package that provides support for genomic-enabled prediction and mate selection for plant and animal breeding: `library(genomicMateSelectR)`.

The primary functions of `genomicMateSelectR` predict the means and variances in performance among progeny of crosses based on parent data in order to support selection of mates in a breeding program. Supports diploid organisms with phased, chromosome- or linkage-group ordered biallelic marker data, and a centimorgan-scale genetic map. Additional functions automate cross-validation estimation of prediction accuracy, and more.

The package includes what might be too many extra functions. Indeed, it spans the entire pipeline used for **genomic mate selection** for the [NextGen Cassava Breeding programs](https://www.nextgencassava.org/) using data downloads from [Cassavabase](https://www.cassavabase.org). Functions to: automate cross-validation procedures, compute accuracy on a selection index, make predictions of individual and cross performances on a multi-trait selection index. Prediction models including additive-effects only ("A"), additive-plus-dominance ("AD") and a directional dominance model ("DirDom"). Functions used in cleaning and curating breeding pipeline field data *plus* handling and imputing genomics data are also included, but most users are not likely to find these of interest. Not everything has been equally tested or documented at this stage.

## Installation

You can install [**genomicMateSelectR** package from my GitHub](https://www.github.com/wolfemd/genomicMateSelectR/) with:

``` {.r}
devtools::install_github("wolfemd/genomicMateSelectR", ref = 'master') 
```

## Get Started

**CHECK OUT THE NEW VIGNETTE! --> [Getting starting predicting crosses](articles/start_here.html)**

**More to come!**

## Feature highlights

-   Allows for parents to be of arbitrary heterozygosity/homozygosity (outbred or inbred)
-   Predicts the **additive** *and* **dominance** genetic variances in the $F_1$
-   Support for a directional dominance model (`modelType="DirDom"`) to incorporate genome-wide homozygosity-effects (inbreeding) into predictions.
-   Predicts genetic **variances** *and* ***co***-**variances**.
-   Support for multi-trait **selection index** via `selInd=TRUE` and `SIwts=` arguments.
-   Functions to implement parent-wise cross-validation as described in the paper, also for "standard" cross-validation
-   Handles simple (one trait, one cross) predictions, but built for complex (multi-trait, many crosses) scenarios.
-   Single estimate of marker effects from REML or MCMC (posterior mean effects --> predicts "variance of posterior means") supported
-   Posterior Mean Variance (PMV) *also* supported: the estimator of Lehermeier et al. 2017b which computes the predicted variance across a sample of marker effects, e.g. the thinned MCMC samples, usually stored on disk. For the multi-trait case, a multivariate Bayesian model is required as only a marker effects for each trait must be computed on the same Gibbs chain. *New version of cross variance predictions are not optimized for PMV so users wanted that are recommended to return to `predCrossVar`*.
-   The relatively stable functions used to implement imputation and data cleaning for the NextGen Cassava Breeding GS programs.
-   Various data helper functions

## Core function overview

| **Genomic mate selection functions** | Top level functions to predict the performance of individual genotypes and the usefulness of potential crosses.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|--------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `runGenomicPredictions()`            | Run GBLUP model using `sommer::mmer`, potentially on multiple traits. Returns genomic BLUPs (GEBV and GETGV). If requested, returns backsolved marker effects (equivalent to ridge regression / SNP-BLUP).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `predictCrosses()`                   | Predict potentially for multiple traits, the means, variances and trait-trait covariances in a set of user-requested crosses-to-evaluate. Output enables easy ranking of potential crosses. Potentially computes the usefulness criteria, $UC\_{parent}$ and $UC\_{variety}$. Provides users the option to predict cross usefulness (means and variances) on a linear multi-trait selection index, taking into account trait-trait covariances within each cross, using a set of user-supplied weights. Utilizes the functions `predCrossVars()` and `predCrossMeans()` under-the-hood. This function is designed to work with `runGenomicPredictions()` , taking SNP-effects matrices output from that function. |

| **Cross-validation functions** | Functions to automate cross-validation                                                                                                               |
|--------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| `runParentWiseCrossVal()`      | Assess the accuracy of predicted previously unobserved crosses.                                                                                      |
| `runCrossVal()`                | Run k-fold cross-validation and assess the accuracy of predicted previously unobserved genotypes (individuals) based on the available training data. |

| **Cross-prediction functions** | Functions that predict cross means and variances |
|--------------------------------|--------------------------------------------------|
| `predCrossVars()`              | Predict cross variances and covariances          |
| `predCrossMeans()`             | Predict cross means                              |

## NextGen Cassava GS pipeline functions

In addition, a host of functions developed to support processing both the field trial data stored on the [Cassavabase](https://www.cassavabase.org/) and the genotyping/genomics data (imputation, file conversions).

There are two "families" of functions distinguished in the [Reference](reference/index.html) as "cassavabase_pheno_pipeline" and "imputation_functions".

## Relationship to `predCrossVar`

`library(genomicMateSelectR)` descends from and extends the [**predCrossVar**](https://wolfemd.github.io/predCrossVar/) R package. `library(predCrossVar)` was build alongside an [an initial study](https://www.biorxiv.org/content/10.1101/2021.01.05.425443v1), which showed promising results regarding the prediction of genetic variance in cassava crosses. Subsequent to that initial study, the code was completed and improved. The considerable subsequent analyses leading to the functions collected in `library(genomicMateSelectR)` are completely documented [here](https://wolfemd.github.io/implementGMSinCassava/).
