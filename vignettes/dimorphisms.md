---
title: "Profiling sexual-dimorphisms with teff"
author: 
  - name: Alejandro Caceres
    affiliation: 
    - ISGlobal, Barcelona, Spain
    - Department of Mathematics, Escola d'Enginyeria de Barcelona Est (EEBE) Universitat Politècnica de Catalunya, Barcelona Spain.
    email: alejandro.cacere@isglobal.org 
  - name: Luis A. Perez-Jurado
    affiliation: 
    - Genetics Unit, Universitat Pompeu Fabra, Barcelona, Spain
  - name: Juan R. González
    affiliation: 
    - ISGlobal, Barcelona, Spain
package: teff
output: 
  BiocStyle::html_document:
    number_sections: true
    toc_float: yes
vignette: >
  %\VignetteIndexEntry{Profiling sexual-dimorphisms with teff}
  %\VignetteEngine{knitr::knitr}
  \usepackage[utf8]{inputenc}
---


# Introduction

<code>teff</code> is a software package to predict the effect of treating an individual given the individual's profile in some feature data. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of treatment effects is based on inferences using random causal forest as implemented in the package <core>grf</code> by Tibshirani et al. 

Here, we show how to use the package to estimate transcriptomic profiles with high immune sexual dimorphism. We use sex as a treatment and therefore search for the transcription profiles for which the effect of sex is high on the amount of immune cells in blood.

We have inferred the amount of T cells in blood in the GTEx individuals and have preselected the transcriptomic features with a differential expression analysis of the interaction between T cell count and sex (more details in Caceres et al. 2021, under preparation). 


# Data



```r
library(teff)
```

The <core>tcell</code> data set is a list that contains the  
<core>teffdata</core> with variables: <core>eff</core> for the effect, namely the infered T-cell count in blood, <core>t</core> for treatment; 1:male, 2:female, and covariates such as <core>age</core>, <core>bmi</core> and  <core>cov5</core>, which is a surrogate variable estimated for the transcription data, and associated to latent technical differences in their measurement.  


```r
names(tcell)
```

```
## [1] "teffdata" "features"
```

```r
head(tcell$teffdata)
```

```
##                               eff t age   bmi         cov5
## GTEX-N7MS-0007-SM-26GMV 23.111463 1  61 26.28  0.076541692
## GTEX-NFK9-0006-SM-3GACS 14.026118 1  49 27.12  0.011393626
## GTEX-NL3G-0007-SM-4SOIF 14.195506 2  67 40.35  0.008561203
## GTEX-NL4W-0006-SM-2I3GH  1.531071 1  53 29.29  0.103587529
## GTEX-NPJ7-0006-SM-3GACR 15.961988 2  68 29.23 -0.033921612
## GTEX-N7MT-0007-SM-3GACQ 25.003704 2  68 24.41 -0.044006008
```
The transcriptomic feature data contains the expression levels of 14, found significant in the differential expression analysis of the interaction between T cell count and sex. 


```r
dim(tcell$features)
```

```
## [1] 428  14
```

```r
head(tcell$features)
```

```
##                             DDX3Y    DDX3X     KDM5D    KDM5C      PRKY
## GTEX-N7MS-0007-SM-26GMV 7.1325489 8.877976 6.7455258 6.745526 5.8524410
## GTEX-NFK9-0006-SM-3GACS 4.1589332 6.821898 4.4484399 5.076471 3.3109363
## GTEX-NL3G-0007-SM-4SOIF 0.8516047 7.928420 0.8516047 6.684495 0.8516047
## GTEX-NL4W-0006-SM-2I3GH 3.1063862 4.691349 3.1063862 4.691349 1.5214237
## GTEX-NPJ7-0006-SM-3GACR 1.0960979 6.824018 1.0960979 6.710808 1.0960979
## GTEX-N7MT-0007-SM-3GACQ 1.2111039 8.031283 1.2111039 7.360851 1.2111039
##                             PRKX    RPS4Y1     RPS4X    TXLNGY    TXLNG
## GTEX-N7MS-0007-SM-26GMV 5.852441 6.2150111  8.536939 5.8524410 4.630049
## GTEX-NFK9-0006-SM-3GACS 3.310936 6.1984616  8.283629 3.3109363 2.573971
## GTEX-NL3G-0007-SM-4SOIF 4.758495 0.8516047  9.195901 0.8516047 3.173533
## GTEX-NL4W-0006-SM-2I3GH 1.521424 4.6913487  7.076013 1.5214237 1.521424
## GTEX-NPJ7-0006-SM-3GACR 4.266023 1.0960979 10.745354 1.0960979 2.681060
## GTEX-N7MT-0007-SM-3GACQ 5.298567 1.2111039  9.616245 1.2111039 4.018459
##                             USP9Y    USP9X      XIST      TSIX
## GTEX-N7MS-0007-SM-26GMV 6.2150111 8.330488 3.0450861 3.0450861
## GTEX-NFK9-0006-SM-3GACS 2.5739707 5.381326 0.9890082 0.9890082
## GTEX-NL3G-0007-SM-4SOIF 0.8516047 5.375167 8.4882293 8.3028158
## GTEX-NL4W-0006-SM-2I3GH 1.5214237 3.843352 3.1063862 3.1063862
## GTEX-NPJ7-0006-SM-3GACR 1.0960979 5.739954 6.5879510 6.5223627
## GTEX-N7MT-0007-SM-3GACQ 1.2111039 6.702957 8.3091359 8.1536184
```

# Analysis

We then use <code>predictteff</code> to estimate the effect of sex on T cell count in a subsample of test individuals. The function randomly selects 80\% of individuals to grow the forest using the feature data, adjusted by covariates. The predictor is applied on the \20% of left-out individuals, who are used to estimate the effect of sex on each of them, given their gene expression levels across the genes. For the immune sexual dimorphism, the genes within the features belong to sex chromosomes. For those that are coupled with homologous genes between the sexual chromosomes, we estimate the average levels between homologs and use those as the features in the random forest. The matrix <code>homologous</code> couples the homologs between rows.  


```r
homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)

pred <- predicteff(tcell, featuresinf=homologous, profile=TRUE)

pred
```

```
## object of class: pteff 
## Estimated treatment effects in $predictions: 
##  [1]   6.6485833  -5.7028466  -3.8289295  -3.8192694  -3.7057803   4.4653225
##  [7]   6.6794771   2.5649402 -12.8958318 -11.9609187   6.0041126  -1.8374791
## [13]  -6.3450333 -12.4676735 -10.0255919   2.4823997  -5.5275620  -6.9856559
## [19]  -5.5720955   4.0512810  14.0750727  10.3782952   3.0360949   6.9510375
## [25]  -5.1137459   6.0374294   2.1693907   9.7728061  13.3242438  -2.4570938
## [31]  -5.9521672 -10.6135121   1.1786464  -6.5074171  -6.7046309  -6.3080904
## [37]   4.8647779   2.8521892 -11.8356216 -12.4006121  -0.7404155  -6.7135344
## [43]  -5.8881988  -4.3798962  -5.7743037  -1.0602535   3.8517047  10.4898261
## [49]  -2.2121442   4.5478292  -5.5397219   8.7262339  -4.0027492  14.0994602
## [55]  12.2207131  -3.9611612  -4.0048134   1.6783265  -3.0231956  -3.5509816
## [61]  -5.9151559  -6.5992608   9.7602988 -12.6630520   1.1142965  -0.7789992
## [67]  -6.9613671  -5.1730401 -10.6536716  -0.6666106   0.7942378  14.1245064
## [73]  -0.3412159   4.4328506 -10.1408222  -9.0526475  10.5106305  -2.7657815
## [79] -12.4594603   2.5919673  -9.6278145   7.4425594   7.4480063  -2.7459211
## [85]   0.9192135
```
We can plot the prediction with confidence intervals. Treated refers to females, and not treated to males. 

```r
plot(pred)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

For each individual <code>predicteff</code> predicts a sex effect with 95\% confidence intervals. Therefore, individuals that have a significantly low effect of sex are those for which their specific transcription levels would result in lower T cell count when observed in females than in males. Those with significantly high effects of sex are those for which their transcription profiles would result in higher T cell count when observed in females than males.  

<code>predicteff</code> allows for a parameter <code>profile</code> that when it is set to <code>TRUE</code> it attempts to create a common profile of gene expression levels across all individuals with significantly low sex effect. It also created a profile for the group of individuals with significantly high sex effects.  


```r
pred$profile
```

```
## $profhigh
##      PRKY-PRKX TSIX-TSIX KDM5D-KDM5C USP9Y-USP9X RPS4Y1-RPS4X XIST-XIST
## [1,]     FALSE      TRUE       FALSE       FALSE        FALSE      TRUE
##      TXLNGY-TXLNG DDX3Y-DDX3X
## [1,]        FALSE       FALSE
## 
## $proflow
##      PRKY-PRKX TSIX-TSIX KDM5D-KDM5C USP9Y-USP9X RPS4Y1-RPS4X XIST-XIST
## [1,]      TRUE     FALSE        TRUE        TRUE         TRUE     FALSE
##      TXLNGY-TXLNG DDX3Y-DDX3X
## [1,]         TRUE        TRUE
```

The profile for high dimorphism means, for instance, that the group of individuals who at the same time up-regulate XIST and TSIX and down-regulate PRKY-PRKX, TSIX-TSIX, KDM5D-KDM5C, USP9Y-USP9X, RPS4Y1-RPS4X, and DDX3Y-DDX3X presents more T cell counts in females than in men. 


The sexual dimorphic profiles of immune cell count can be used to target individuals from other studies and test whether the profiles present different sex-associated risks or survival of different diseases. The following code shows how to set up data for the transcriptomic study in arthritis GSE17755 from GEO.  


```r
library(GEOquery)
gsm <- getGEO("GSE17755")
gsm <- gsm[[1]]

data4teff <- feateff(gsm, tname="gender:ch1",
                     reft=c("male", "female"),
                     effname="disease:ch1",
                     refeff=c("healthy","arthritis"),
                     covnames="age:ch1", covtype="n",
                     sva=TRUE, UsegeneSymbol=TRUE)
```

The result has been saved in the data structure <code>data4teff</code> within the <code>teff</code> package, with only the genes within the sexual dimorphism profiles. We see how the data structure is a list with fields <core>features</core> and  <core>teffdata</core>


```r
data(data4teff)
names(data4teff)
```

```
## [1] "features" "teffdata"
```

```r
names(data4teff$teffdata)
```

```
##  [1] "eff"     "t"       "age.ch1" "cov1"    "cov2"    "cov3"    "cov4"   
##  [8] "cov5"    "cov6"    "cov7"    "cov8"    "cov9"    "cov10"   "cov11"  
## [15] "cov12"   "cov13"   "cov14"   "cov15"   "cov16"   "cov17"   "cov18"  
## [22] "cov19"   "cov20"   "cov21"   "cov22"   "cov23"   "cov24"   "cov25"
```

with variables: <core>eff</core> for the effect, namely artrhitis diagnose, <core>t</core> for treatment; 1:male, 2:female, and covariates such as <core>age</core>, and surrogate variables.


The <core>target</core> function classifies individuals into their sexually dimorphic groups and produces a plot on the targeting, based on gene expression data available. It also tests the significance of the association of the effect (arthritis disease) with the interaction between sex and the groups of sexual dimorphism. 


```r
res <- target(data4teff, pred, plot=TRUE, effect="highandlow", featuresinf=homologous, nmcov="age.ch1", model="binomial")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

```r
res
```

```
## object of class: tarteff 
## 
## classification into 
##   low treatment effect: -1
##   neutral: 0
##   high treatment: 1
##  
##  -1   0   1 
##  61 104  57 
## 
## interaction fitted model: binomial
##   Estimate Std. Error    z value   Pr(>|z|) 
##  0.5460645  0.4813214  1.1345110  0.2565803
```

When the outcome is continuous a plot of the interaction can also be obtained. Here, we show the targeting of the in the <code>tcell</code> data set, where the effect is the T Cell count. Clearly, for this case, the association of the effect with the interaction between sex and the groups of sexual dimorphism is high because this data and this effect were used to infer the groups. 



```r
data(tcell)
homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
pf <- predicteff(tcell, featuresinf=homologous, profile=TRUE)
res <- target(tcell, pf, effect="highandlow", featuresinf=homologous, nmcov="age", model="log2")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
res
```

```
## object of class: tarteff 
## 
## classification into 
##   low treatment effect: -1
##   neutral: 0
##   high treatment: 1
##  
##  -1   0   1 
##  79 288  61 
## 
## interaction fitted model: log2
##     Estimate   Std. Error      t value     Pr(>|t|) 
## 1.630906e+01 2.178184e+00 7.487459e+00 4.121995e-13
```

the highly significant interaction is illustrated in the plot

```r
plot(res)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

