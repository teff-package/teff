## ---- message=FALSE, warning=FALSE--------------------------------------------
library(teff)

## -----------------------------------------------------------------------------
names(tcell)
head(tcell$teffdata)

## -----------------------------------------------------------------------------
dim(tcell$features)
head(tcell$features)

## -----------------------------------------------------------------------------
homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)

pred <- predicteff(tcell, featuresinf=homologous, profile=TRUE)

pred


## -----------------------------------------------------------------------------
plot(pred)

## -----------------------------------------------------------------------------
pred$profile

## ---- eval=FALSE--------------------------------------------------------------
#  library(GEOquery)
#  gsm <- getGEO("GSE17755")
#  gsm <- gsm[[1]]
#  
#  data4teff <- feateff(gsm, tname="gender:ch1",
#                       reft=c("male", "female"),
#                       effname="disease:ch1",
#                       refeff=c("healthy","arthritis"),
#                       covnames="age:ch1", covtype="n",
#                       sva=TRUE, UsegeneSymbol=TRUE)
#  
#  

## -----------------------------------------------------------------------------
data(data4teff)
names(data4teff)
names(data4teff$teffdata)

## -----------------------------------------------------------------------------
res <- target(data4teff, pred, plot=TRUE, effect="highandlow", featuresinf=homologous, nmcov="age.ch1", model="binomial")

res

## -----------------------------------------------------------------------------
data(tcell)
homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
pf <- predicteff(tcell, featuresinf=homologous, profile=TRUE)
res <- target(tcell, pf, effect="highandlow", featuresinf=homologous, nmcov="age", model="log2")
res

## -----------------------------------------------------------------------------
plot(res)

