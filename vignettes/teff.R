## ---- message=FALSE, warning=FALSE--------------------------------------------
library(teff)

## -----------------------------------------------------------------------------
names(psoriasis)
head(psoriasis$teffdata[,1:10])

## -----------------------------------------------------------------------------
dim(psoriasis$features)
head(psoriasis$features[,1:10])

## -----------------------------------------------------------------------------
pso <- predicteff(psoriasis, dup=TRUE)

pso

## -----------------------------------------------------------------------------
plot(pso)

## -----------------------------------------------------------------------------
plot(pso, rk =pasiw12, xlab="Observed PASI week 12")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(drc)

response <- pasiw12[pso$subsids]
predictions <- pso$predictions
treatment <- factor(pso$treatment, labels = c("placebo", "brodalumab"))

mod <- drm(response*100~predictions, treatment, fct=LL.3())

plot(mod, pch=16, col=c("orange", "blue"), 
     legendPos=c(0.36,-0.25), 
     ylab="PASI improvement week12", 
     xlab="Predicted treatment effect")

