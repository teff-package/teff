## ---- message=FALSE, warning=FALSE, echo=FALSE--------------------------------
library(teff)

## ---- echo=FALSE--------------------------------------------------------------
Sys.setenv(VROOM_CONNECTION_SIZE=5000072)

## -----------------------------------------------------------------------------
names(tcell)
head(tcell$teffdata)

## -----------------------------------------------------------------------------
dim(tcell$features)
head(tcell$features)

## -----------------------------------------------------------------------------
homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C",
                      "PRKY","PRKX","RPS4Y1","RPS4X",
                      "TXLNGY", "TXLNG", "USP9Y", "USP9X",
                      "XIST", "XIST", "TSIX", "TSIX"), nrow=2)

pred <- predicteff(tcell, featuresinf=homologous, profile=TRUE)

pred

## -----------------------------------------------------------------------------
plotPredict(pred, 
            ctrl.plot=list(lb=c("Male", "Female"),
                           wht="topleft", whs = "bottomright"))

## -----------------------------------------------------------------------------
pred$profile

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  library(GEOquery)
#  gsm <- getGEO("GSE17755", AnnotGPL = TRUE)
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
res <- target(data4teff, pred, plot=TRUE, effect="positive", featuresinf=homologous, nmcov="age.ch1", model="binomial")

res

## -----------------------------------------------------------------------------
data(tcell)
homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
pf <- predicteff(tcell, featuresinf=homologous, profile=TRUE)
res <- target(tcell, pf, effect="positiveandnegative", featuresinf=homologous, nmcov="age", model="log2")
res

## -----------------------------------------------------------------------------
plotTarget(res,  labs=c("Sexual dimorphism", "T cell count", "Condition", "Male", "Female"))

## -----------------------------------------------------------------------------
boxPlot(res, 
        labs=c("Sexual dimorphism", "T cell count", "Condition", "Male", "Female"), 
        lg="bottomleft")

