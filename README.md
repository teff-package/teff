# teff (treatment effects estimation)

This is an R Package for estimating which individuals would respond best to a treatment on an outcome, given their values in several feature data.

## What it does

<code>teff</code> is a software package to estimate the expected effect of treating a **single** individual given the individual's profile in some feature data. The individual treatment effect is the expected difference in outcome between two treatments (treatment/no-treatment, drug/placebo) when feature data are kept constant. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of individual treatment effects is based on causal inference using causal random forest as implemented in the package [grf](https://github.com/grf-labs/grf) by Tibshirani et al.

Feature profiles associated to positive and negative effects of treatment can then be extracted, based on the feature data of individuals with significant treatment effects. If such profiles can be identified then new individuals can be targeted and classified into groups of positive, negative or neutral treatment effects. When the outcome and treatment data are also available for these new individuals, the package can test whether the association between the treatment and outcome is indeed different across groups of associated positive and negative treatment effects. 

Applications include:

- Forecasting treatment effect of brodalumab treatment of psoriasis from gene expression data at baseline. See [application](https://alejandro-isglobal.github.io/teff/teff.html)

- Extracting gene expression profiles for which men and women have high differences on immune cell count in blood. In this application sex is considered as a treatment. See [application](https://alejandro-isglobal.github.io/teff/dimorphisms.html) .

## How it does it

The classification is based on the application of causal random forest to identify the individuals with significant treatment effects. Individuals with significant treatment effects are those whose confidence intervals for the treatment effect estimate do not overlap 0 (or a defined threshold). 

When possible, single consensus profiles of individuals with positive, and negative, treatment effects are obtained from majority votes of features adjusted for covariates.

## Additional functions

The package includes functions to extract feature and treatment data from transcriptomic and methylomic studies. 

## Install

Install with 


<code>
library(devtools)
</code>
</br><code>
install_github("teff-package/teff")
</code>


</br>See vignettes  [teff](https://alejandro-isglobal.github.io/teff/teff.html) and [dimorphisms](https://alejandro-isglobal.github.io/teff/dimorphisms.html) 
online
