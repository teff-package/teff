#' Predicts treatmet effects of individuals in a randon subsample
#'
#' @details This function sets up feature and treatment-effects data, fits
#' random causal forest and identify the individuals with significant
#' treatment effects. Individuals with significant treatment effects
#' are considered for those whose confidence intervals for the treatment
#' estimate do not overlap 0. Single consensus profiles of individuals with
#' high, and low, treatment effects are obtained from majority votes of
#' adjusted features, binarized over the population means.
#'
#' The result is two profiles, associated with high and low
#' treatment effects, given by logical vectors across the
#' features. The logical value of a given profile at feature indicates whether the
#' adjusted feature of a new individual should be higher than the feature population
#' mean if the individual is successfully targeted by the profile. See \link[teff]{taget}.
#'
#' @export
#' @param x a \code{list} with a fields \code{teffdata} and \code{features}.
#' \code{teffdata} is a \code{data.frame} (or \code{matrix}) with the treatment
#' \code{$t} and effect \code{$eff} variables, and covariates across
#' subjects. \code{features} is a \code{matrix} which the profiling
#' is of subjects if performed. The features are adjusted for the covariates
#' on the \code{teffdata} before fitting the causal random forest forest.
#'
#' @param featuresinf a \code{vector} of characters with the names of the
#' features to be used; a \code{matrix} of characters whose columns are
#' names of features whose values will be averaged (Default:NULL).
#' @param cores an integer with the number of cores for parallel
#' computation (Default:1)
#' @param seed and integer with the random seed for splitting data into
#' train (80%) and test (20%) sets.
#' @param plot.overlap a logical. If \code{TRUE} then it plots the overlap of
#' adjusted feature data across treatments. Parameter available for less
#' than 20 features (Default: \code{FALSE}).
#' @param quant a number from 0 to 1 with the quantile of features to be selected
#' with top information score from the causal forest.
#' By default it selects all the features (Default: Inf).
#' @param profile a logical. If \code{TRUE} then it estimates a profile of binarized
#' feature data for the individuals with significantly high and low treatment effects, respectively.
#' @param dup a logical that indicates whether the feature and teff data should
#' be duplicated in case of small datasets.
#' @return  a \code{list} of class \code{pteff} with fields:
#' \describe{
#' \item{predictions:}{a \code{vector} with the estimated treatment effect
#' of the individuals in the test set.}
#' \item{featurenames:}{a \code{vector} with the names of the features used.}
#' \item{cl:}{a \code{vector} with the lower limit of the 95% confidence
#' intervals for the estimated treatment effect.}
#' \item{cu:}{a \code{vector} with the upper limit of the 95% confidence
#' intervals for the estimated treatment effect.}
#' \item{subsids:}{a \code{vector} with ids of subjects in the test set.}
#' \item{treatment:}{a \code{vector} with treatment effect in the test set.}
#' \item{profile:}{a \code{list} with fields \code{profhigh} and \code{proflow}
#' that are matrices with binarized feature data for the individuals with
#' significantly high and low treatment effects, respectively.}
#' }
#' @export
#'
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' predicteff(tcell, featuresinf=homologous, profile=TRUE)
#'
predicteff <- function(x,
                    featuresinf=NULL,
                    cores=1,
                    seed=1234,
                    plot.overlap=FALSE,
                    quant=Inf,
                    dup=FALSE,
                    profile=FALSE){

  ##Data set up
  ########################

  teffdata <- x$teffdata
  features <- t(x$features)
  subsids <- colnames(features)

  #redefine features on X by taking the average across featuresinf rows
  if(is.matrix(featuresinf)){

    X <- lapply((1:ncol(featuresinf)), function(i)
      colMeans(features[featuresinf[,i], ]))

    X <- do.call(cbind,X)
    colnames(X) <-  sapply(1:ncol(featuresinf), function(i) paste0(featuresinf[,i], collapse="-"))

    if(plot.overlap){
      wh <- which(colSums(sapply(data.frame(featuresinf), duplicated))==0)
      colnames(X) <-  paste(featuresinf[1,],c(featuresinf[2,wh]," "," "), sep="-")
    }

    nms <- colnames(X)
    featuresinf <- nms

  }else{
    if(is.null(featuresinf)) featuresinf <- rownames(features)
    selfeatures <- which(rownames(features)%in%featuresinf)
    X <- t(features[as.numeric(selfeatures), ])
    nms <- colnames(X)
  }

  #obtain residuals on XX of features data after adjusting for covariates in teffdata
  XX <- parallel::mclapply(1:ncol(X), function(i){
    covvars <- colnames(teffdata)
    xi <- X[,i]
    fla <- paste("xi ~", paste(covvars, collapse="+"))

    lm(as.formula(fla) , data=data.frame(teffdata))$residuals
  }, mc.cores = cores)


  #redefine feature data
  X <- do.call(cbind, XX)
  colnames(X) <- nms


  #Formating data for RCF
  ########################

  #effect
  Y <- teffdata[,"eff"]

  #treatment
  W <- as.numeric(as.factor(teffdata[,"t"]))==2

  #randomly selection of a test comprising 20% of subjects
  set.seed(seed)
  sm <- sample(1:nrow(X),floor(nrow(X)*0.2))
  smt <- (1:nrow(X))[-sm]

  if(dup){
    sm <- c(sm,sm)
    smt <- c(smt,smt)
  }

  X.test <- X[sm,]
  Y.test <- Y[sm]
  W.test <- W[sm]

  #randomly select a training
  X.train <- X[smt,]
  Y.train <- Y[smt]
  W.train <- W[smt]

  #plots covariate overlap
  ########################

  if(plot.overlap){
    if(ncol(X)>20) stop("plot.overlap not allowed for number of features > 20")

    print("... plots of t:eff interactions across tretments in interactions.pdf \n")

    pdf("./interactions.pdf")
      cc <- W.train
      cc[W.train] <- "red"
      cc[!W.train] <- "black"

      for(ii in 1:ncol(X.train)){

        plot(X.train[,ii], log2(Y.train+1), col=cc, main=colnames(X.train)[ii],pch=20, xlab="Mean expression residual", ylab="log2 cell abundancy")
        abline(lm(log2(Y.train+1)[W.train] ~ X.train[W.train,ii]),lwd=1.5, lty=2, col="red")
        abline(lm(log2(Y.train+1)[!W.train] ~ X.train[!W.train,ii]),lwd=1.5, lty=2)
        legend("topleft",legend=c("Not treated", "Treated"), col=c("black", "red"), pch=20, cex=0.7)
      }


    dev.off()


    pdf("./overlap.pdf")
      print("... plots of covariate overlap across tretments in overlap.pdf \n")
      for(ii in 1:ncol(X.train)){

        int <- seq(min(X.train)-0.5,max(X.train)+0.5, 0.5)

        hist(X.train[W.train,ii], freq=FALSE, br=int,
             border="red", main=gsub("\n ","-",colnames(X.train)[ii]),
             xlab="Mean expression residual", ylab="Density", ylim=c(0,0.9))

        hist(X.train[!W.train,ii], add=TRUE, freq=FALSE, br=int)

        legend("topright", legend=c("Not treated", "Treated"), col=c("black", "red"),
               lty=1,cex=0.7, bty = "n" )
      }

    dev.off()
  }

  #fit RCF
  ########################

  tau.forest <- grf::causal_forest(X.train, Y.train, W.train, num.trees = 4000)
  tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
  names(tau.hat$predictions) <- subsids[sm]

  sigma.hat <- sqrt(tau.hat$variance.estimates)

  v1 <- grf::variable_importance(tau.forest)
  o.imp <- order(-v1)
  imp <- o.imp[1]
  featurecommon <- data.frame(featurenames=colnames(X), imp=v1)[o.imp,]

  featureimp <- lapply(featuresinf, function(ll){
    #select features results in each study for common feature symbols
    out <- featurecommon[which(as.character(featurecommon$featurenames)%in%ll), ]
    #select probe with highest score
    if(length(nrow(out))!=0)
      out <- out[order(-out[,2])[1], ]
    out
  })

  featureimp <- do.call(rbind,featureimp)
  featureimp <- featureimp[order(featureimp$imp, decreasing = TRUE),]

  #extract confidence intervals
  cl <- tau.hat$predictions - 1.96*sigma.hat
  cu <- tau.hat$predictions + 1.96*sigma.hat
  #output
  res <- list(predictions=tau.hat$predictions, featurenames=featureimp,  cl=cl, cu=cu, subsids=subsids[sm], treatment = W.test*1)


  #add profiling to output
  if (profile){

    sighethigh <- (cl>0) + 1
    sighetlow <-  (cu<0) + 1


    #get profile for top featurenames define by quant
    cutoff <- max(which(cumsum(featureimp$imp) < quant))
    colimp <- featureimp$imp[1:cutoff]

    #select trasncription profiles (residuals) for featurenames passing the cutoff
    whimp <- as.numeric(rownames(featureimp))
    Xscale <- scale(X.test)[, whimp[1:cutoff]]

    #Extract a binarized profiles across important featurenames, where featurenames up or
    #downregulated in their study
    profall <- lapply(1:ncol(Xscale), function(i) (Xscale[,i] > 0))
    profall <- do.call(cbind, profall)

    #Profiles of individuals with significant heterogeneous treatment-effects
    profhigh <- profall[sighethigh==2,]
    proflow <- profall[sighetlow==2,]


    if(is.matrix(profhigh)==FALSE)
      profhigh <- matrix(profhigh, ncol=ncol(profall))

    if(is.matrix(proflow)==FALSE)
      proflow <- matrix(proflow, ncol=ncol(profall))

    #summarize subject profiles into a single one representing individuals with
    #high treatment effects
    profhigh <- matrix(colMeans(profhigh, na.rm=TRUE)>0.5, nrow=1)

    #and low treatment effects
    proflow <- matrix(colMeans(proflow, na.rm=TRUE)>0.5, nrow=1)

    colnames(profhigh) <- colnames(Xscale)
    colnames(proflow) <- colnames(Xscale)

    ##########
    #gather output
    res <- list(predictions=tau.hat$predictions, featurenames=featureimp,  cl=cl, cu=cu, subsids=subsids[sm], treatment = W.test*1, profile=list(profhigh=profhigh,proflow=proflow))
  }


  attr(res, "class") <- "pteff"

  return(res)
}
