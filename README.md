# teff (treatment effects prediction)

This is an R Package for predicting which individuals would respond best to a treatment on an outcome, given their values in several feature data.

## What it does

<code>teff</code> is a software package to predict the effect of treating a **single** individual given the individual's profile in some feature data. The effect of treatment is the estimated difference of an outcome between treating and and not treating when the feature data are kept constant. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of treatment effects is based on inferences using random causal forest as implemented in the package [grf](https://github.com/grf-labs/grf) by Tibshirani et al.

Feature profiles associated to positive and negative effects of treatment can then be extracted, based on the feature data of individuals with significant treatment effects. If such profiles can be identified then new individuals can be targeted and classified into groups of positive, negative or neutral treatment effects. When the outcome and treatment data are also available for these new individuals, the package can test whether the association between the treatment and outcome is indeed different across groups of associated positive and negative treatment effects. 

Applications include:

- Forcasting treatment effect of brodalumab treatment of psoriasis from gene expression data at baseline. See [application](https://alejandro-isglobal.github.io/teff/teff.html)

- Extracting gene expression profiles for which men and women have high differences on immune cell count in blood. In this application sex is considered as a treatment. See [application](https://alejandro-isglobal.github.io/teff/dimorphisms.html) .

## How it does it

The classification is based on the application of random causal forest to identify the individuals with significant treatment effects. Individuals with significant treatment effects are considered for those whose confidence intervals for the treatment estimate do not overlap 0. 

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
