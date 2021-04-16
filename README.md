# teff (treatment effects prediction)

This is an R Package for predicting which individuals would response best to treatment given their values in several feature data.

## What it does

<code>teff</code> is a software package to predict the effect of treating a single individual given the individual's profile in some feature data. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of treatment effects is based on inferences using random causal forest as implemented in the package <core>grf</code> by Tibshirani et al.

With the extracted profiles, new individuals according to their feature data can be targeted and classified into groups where the treatment is significantly positive or negative on the effect. If treatment and effect data is available for these new individuals the package can test
whether the association between the treatment and the effect is indeed different across groups of associated positive and negative treatment effects. 

Applications include:

- Predicting treatment effect of brodalumab treatment of psoriasis from gene expression data at baseline. 

- Extracting gene expression profiles for which men and women have high differences on immune cell count in blood. In this application sex is considered as a treatment. 

## How it does it

The classification is based on the application of random causal forest to identify the individuals with significant treatment effects. Individuals with significant treatment effects are considered for those whose confidence intervals for the treatment estimate do not overlap 0. 

When possible, single consensus profiles of individuals with positive, and negative, treatment effects are obtained from majority votes of features adjusted for covariates.

## Additional functions

The package include functions to extract feature and treatment data from transcriptomic and methylomic studies. 

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
