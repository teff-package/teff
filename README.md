# teff (treatment effects classification)

This is an R Package for predicting which individuals would response best to treatment given their values in several feature data.

## What it does

<code>teff</code> is a software package to predict the effect of treating an individual given the individual's profile in some feature data. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of treatment effects is based on inferences using random causal forest as implemented in the package <core>grf</code> by Tibshirani et al. 

With the extracted profiles, new individuals with feature data can be targeted and classified. If treatment and effect data is available for these new individuals the package can test
whether the association between the treatment and the effect is indeed different across subpopulations. 

Applications include:

- Predicting treatment effect of brodalumab treatment of psoriasis from gene expression data at baseline. 

- Extracting gene expression profiles for which men and women have high differences on immune cell count in blood. In this application sex is considered as a treatment. 

## How it does it

The classification is based on the application of random causal forest to identify the individuals with significant treatment effects. Individuals with significant treatment effects are considered for those whose confidence intervals for the treatment estimate do not overlap 0. 

When possible, single consensus profiles of individuals with high, and low, treatment effects are obtained from majority votes of features adjusted for covariates.

## Additional functions

The package include functions to extract feature and treatment data from transcriptomic and methylomic studies. 

## Install with

<code>
library(devtools)
install_github("teff-package/teff")
</code>

and see vignettes [teff]((https://alejandro-isglobal.github.io/teff/teff.html) and [dimorphisms](https://alejandro-isglobal.github.io/teff/dimorphisms.html) 
