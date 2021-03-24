---
title: "Predicting treatment effects with teff"
author: 
  - name: Alejandro Caceres
    affiliation: 
    - ISGlobal, Barcelona, Spain
    - Department of Mathematics, Escola d'Enginyeria de Barcelona Est (EEBE) Universitat Politècnica de Catalunya, Barcelona Spain.
    email: alejandro.cacere@isglobal.org 
  - name: Juan R. González
    affiliation: 
    - ISGlobal, Barcelona, Spain
package: teff
output: 
  BiocStyle::html_document:
    number_sections: true
    toc_float: yes
vignette: >
  %\VignetteIndexEntry{Predicting treatment effects with teff}
  %\VignetteEngine{knitr::knitr}
  \usepackage[utf8]{inputenc}
---


# Introduction

<code>teff</code> is a software package to predict the effect of treating an individual given the individual's profile in some feature data. The package focuses on transcriptomic features for which surrogate covariates need to be estimated. The estimation of treatment effects is based on inferences using random causal forest as implemented in the package <core>grf</code> by Tibshirani et al. 

Here, we show how to use the package to estimate transcriptomic profiles of psoriasis patients that strongly respond to brodalumab treatment (more details in Caceres et al. 2021, under preparation). Publicly available data has been downloaded from GEO for accession number GSE117468. This is data from a clinical trial of brodalumab treatment of psoriasis patients.


# Analysis



```r
library(teff)
```

We ﬁrst obtained clinical data relating to age, BMI, psoriasis area-and-severity-index (PASI) at baseline and at week 12 after treatment, and brodalumab or placebo treatment. 

The <core>psoriasis</code> data set is a list that contains the  
<core>teffdata</core> with variables: <core>eff</core> for the effect, namely the improvement of PASI score at week 12 with 1: improved PASI after treatment, 0: did not improve PASI after treatment; <core>t</core> for treatment with 0:placebo, 1:brodalumab; and covariates such as <core>age</core>, <core>bmi</core> and  <core>cov5</core>..., which are surrogate variables estimated for the transcription data, and associated to latent technical differences in their measurement.  


```r
names(psoriasis)
```

```
## [1] "features" "teffdata"
```

```r
head(psoriasis$teffdata[,1:10])
```

```
##            t eff age    bmi        cov5        cov6         cov7        cov8
## GSM3300910 1   1  53 20.750 -0.03890177  0.03467942  0.043746932  0.17594309
## GSM3300916 0   1  51 35.235  0.14126324  0.10711059 -0.001265844 -0.06115758
## GSM3300920 0   0  47 35.471 -0.15755427  0.20506609 -0.079110601  0.07190374
## GSM3300928 1   1  38 33.272  0.04814467 -0.08004999 -0.165959371  0.03569621
## GSM3300932 0   0  47 36.553 -0.07157232 -0.11915946 -0.038394625 -0.04629123
## GSM3300936 0   0  64 32.189 -0.30802680 -0.09920047  0.043504869  0.04527702
##                   cov9        cov10
## GSM3300910 0.135287841 -0.177583007
## GSM3300916 0.122662569  0.030069280
## GSM3300920 0.078727884  0.105476906
## GSM3300928 0.009297803 -0.026598687
## GSM3300932 0.044703103 -0.009509721
## GSM3300936 0.061549094  0.094933813
```
The transcriptomic feature data contains the expression levels of 87 genes in none lesional skin at baseline that were found significant in the differential expression analysis of the interaction between PASI improvement and treatment (more details in Caceres et al. 2021, under preparation). 



```r
dim(psoriasis$features)
```

```
## [1] 96 87
```

```r
head(psoriasis$features[,1:10])
```

```
##            204622_x_at 211868_x_at 215565_at 210090_at 215036_at 234884_x_at
## GSM3300910    5.701753    6.212269  4.118702  3.707908  4.646778    6.879799
## GSM3300916    6.147227    4.517039  2.648588  3.408026  3.747902    4.097929
## GSM3300920    7.408559    6.305645  3.864454  3.387308  4.632612    6.389898
## GSM3300928    6.450375    5.298880  2.722913  3.428508  4.526163    4.319476
## GSM3300932    6.163380    5.337662  2.561234  4.372990  4.285132    4.196663
## GSM3300936    6.811362    5.596383  2.915567  4.455007  4.085423    4.628715
##            217281_x_at 207768_at 214973_x_at 217179_x_at
## GSM3300910    5.774919  3.429258    6.339571    4.746741
## GSM3300916    4.549441  3.463358    3.714813    3.619801
## GSM3300920    5.276760  3.234945    4.289022    5.813204
## GSM3300928    4.668434  3.408525    4.088177    3.754877
## GSM3300932    4.275372  3.614708    3.585289    3.405797
## GSM3300936    4.719125  3.678755    3.793925    4.038573
```

We then use <code>predictteff</code> to estimate the effect of treatment on PASI improvement on a subsample of test individuals. The function randomly selects 80\% of individuals to grow the forest using the feature data, adjusted by covariates. The predictor is applied on the \20% of left-out individuals, who are used to estimate the effect of treatment on each of them, given their gene expression levels across all the genes. 



```r
pso <- predicteff(psoriasis, dup=TRUE)

pso
```

```
## object of class: pteff 
## Estimated treatment effects in $predictions: 
##  [1] 0.3216790 0.3731556 0.3370073 0.3153833 0.3923582 0.3457552 0.3543223
##  [8] 0.3098761 0.3436701 0.2782222 0.4316714 0.3584010 0.3480320 0.3873781
## [15] 0.3374915 0.2904020 0.3832157 0.3148482 0.2839521 0.3216790 0.3731556
## [22] 0.3370073 0.3153833 0.3923582 0.3457552 0.3543223 0.3098761 0.3436701
## [29] 0.2782222 0.4316714 0.3584010 0.3480320 0.3873781 0.3374915 0.2904020
## [36] 0.3832157 0.3148482 0.2839521
```

We can plot the prediction of treatment effect with confidence intervals. We have implemented the <code>dup</code> option for handling data sets with few individuals. 


```r
plot(pso)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)


This is a prediction of the probability of treatment for each individual at baseline on nonlesional skin. Not treated refers to patients who ended up receiving placebo, and treated to those under brodalumab. We can compare the prediction with the observed improvement after week 12. 



```r
plot(pso, rk =pasi12w, xlab="Observed PASI week 12")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

We see that while individuals receiving placebo did not improve PASI, some would have strongly benefited from treatment. For individuals who ended up un brodalumab treatment, we can see a strong relationship between predicted response given, by the associated treatment effect at baseline, and the finally observed PASI improvement at week 12.  
