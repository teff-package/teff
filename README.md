# teff (treatment effects) for transcriptomic data

This is an R library for estimating which individuals would respond best to a treatment on an outcome, given their values in gene transcription data. 


## What it does

<code>teff</code> is a software package to estimate the expected effect of treating a **single** individual given the individual's profile in some high dimensional feature data. The individual treatment effect is the expected difference in outcome between two treatments (treatment/no-treatment, drug/placebo) when feature data are kept constant. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of individual treatment effects is based on causal inference using causal random forest (CRF) as implemented in the package [grf](https://github.com/grf-labs/grf) by Tibshirani et al.

The analysis strategy is based on the selection of features, subjects and extraction of a binary profile to positive, or even negative, effects of treatment. Selection of features are obtained by selecting transcription levels whose association with treatment response is significantly modulated by the treatment. Selection of individuals are those with significant treatment effects, as estimated by CRF. An average binarized profile of significant individuals is then constructed.  If such profile can be identified then new individuals can be targeted and classified into groups of positive (or negative) and neutral treatment effects. When treatment and treatment response data are also available for these new individuals, the package can test whether the association between the treatment and outcome is indeed different across groups of associated positive and negative treatment effects. 

Applications include:

- Forecasting treatment effect of brodalumab treatment of psoriasis from gene expression data at baseline. See [application](https://teff-package.github.io/teff/teff.html)

- Extracting gene expression profiles for which men and women have high differences on immune cell count in blood. In this application sex is considered as a treatment. See [application](https://teff-package.github.io/teff/dimorphisms.html) .

## How it does it

Feature selection is performed with transcriptome-wide differential gene expression for the interaction between treatment and response.  

Individual targeting into groups of different treatment effects is based on the application of CRF to identify the individuals with significant treatment effects. Individuals with significant treatment effects are those whose confidence intervals for the treatment effect estimate do not overlap 0 (or a defined threshold). 

When possible, single consensus profiles of individuals with positive, and negative, treatment effects are obtained from majority votes of features adjusted for covariates.

## Additional functions

The package includes functions to extract feature and treatment data from transcriptomic studies. 

## Install

Install with 


<code>
library(devtools)
</code>
</br><code>
install_github("teff-package/teff")
</code>



</br>See [application to Brodalumab treatment of psoriasis](https://teff-package.github.io/teff/SupplementaryCaceres2021.html)


Cite: Alejandro Caceres and Juan R Gonzalez. teff: estimation of Treatment EFFects on transcriptomic data with casual random forest, 2022; under submission. 
