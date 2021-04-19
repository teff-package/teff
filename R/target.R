#' Targets individuals into groups of positive, negative or neutral effect of treatment on an outcome.
#'
#' @details This function uses feature data of individuals to classify them
#' into subpopulations associated with the positive and negative effect of treating them
#' according to a given outcome. The
#' classification is performed by targeting the feature data (adjusted by
#' covariates and binarized) to the profiles computed in \link[teff]{profile}.
#'
#' The function tests whether the classification of the subjects into groups of positive
#' and negative treatment effects modulates
#' the association of the treatment with the outcome, fitting a model for the outcome as function of
#' the interaction between the classification and the treatment.
#'
#' Models on the outcome to test the interaction include general lineal
#' models, beta regression and proportional hazards models.
#'
#' The function can be used to target new individuals not used in the profiling
#' and/or on the effects of other types of outcomes also expected from the treatment.
#'
#' @export
#' @param  x a \code{list} with a fields \code{teffdata} and \code{features}.
#' \code{teffdata} the data for treatment, outcome on which the effect is measured and covariates across
#' subjects. \code{features} data is the data on which the profiling
#' is done. Adjusted features for the variables on the \code{teffdata} are
#' used to fit the forest and extract the profiles of individuals with
#' significant  treatments.
#' @param pteffObject object of class \code{pteff}
#' @param effect \code{character} with one of the values: "positive", "negative", "positiveandnegative", indicating
#' which subpopulations to be targeted: subpopulation with positive treatment effects only, negative
#' treatment effects only or both, respectively (Default: positive).
#' @param featuresinf a \code{vector} of characters with the names of the
#' features to be used; a \code{matrix} of characters whose columns are
#' names of features whose values will be averaged (Default:NULL).
#' @param plot \code{logical} indicating a whether a plot of the targeting will be produced
#' @param lb \code{character} indicating labels for features plot.
#' @param model \code{character} with one of the values corresponding to the family argument of a
#' general linear model for the outcome (\link[stats]{glm}), such as: "gaussian", "binomial", etc. It also allows "beta" for beta
#' regression, "hazard" for proportional hazards model \link[survival]{coxph}, and "log2" for log2
#' transformation of the effects (Default: "gaussian").
#' @param match \code{numeric} value between 0 and 1 indicating the level of the match for the targeing
#' across all the binarized features in the profile (positive or negative). A value of 1 means that all the binarized
#' features must take identical values in the profile (i,e "positive") for an individual to be clsssified
#' in its subpopulation (of positive treatment effects).
#' @param nmcov vector of \code{character} with the names of the covariates to be included in the models.
#' @param cores an integer with the number of cores for parallel
#' computation (Default:1).
#' @param ... additional parameters for the \code{image} function
#' @return   a \code{list} of class \code{tarteff} with fields:
#' \describe{
#' \item{classification:}{classification of individuals into subpopulations of
#' expected positive and/or negative treatment effects}
#' \item{summary.model:}{a \code{summary} of the model used to test the association of the outcome with the
#' interaction of between
#' the profile ("positiveandnegative": positive=1, neutral=0, negative=-1; "positive": positive=1, neutral=0,
#' "negative": negative=1, neutral=0) with the effect.}
#' }
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous, profile=TRUE)
#' res <- target(tcell, pf, effect="positiveandnegative", featuresinf=homologous, nmcov="age", model="log2")
#' res
#'
target <- function(x,
                   pteffObject,
                   effect="positive",
                   featuresinf=NULL,
                   plot=TRUE,
                   lb=NULL,
                   model="gaussian",
                   match=0.8, nmcov=NULL, cores=1, ...){
  teffdata <- x$teffdata
  features <- t(x$features)

  if(effect=="positive"){
    nmfeatures <- unlist(strsplit(colnames(pteffObject$profile$profpositive), "-"))
  }else{
    nmfeatures <- unlist(strsplit(colnames(pteffObject$profile$profnegative), "-"))
  }
  nmfeatures <- unique(nmfeatures)

  selfeatures <- which(rownames(features)%in%nmfeatures)
  #use featuresinf
  if (length(featuresinf) != 0)
    selfeatures <- which(rownames(features)%in%as.vector(featuresinf))

  X <- t(features[as.numeric(selfeatures), ])
  nms <-  colnames(t(features[selfeatures,]))

  #compute mean across multiple variables of the same feature
  ################

  XX <- lapply(unique(nms), function(x)
    if(is.matrix(ncol(X[,x])) == TRUE)
      rowMeans(X[,x])
    else X[,x])
  XX <- do.call(cbind,XX)
  colnames(XX) <- unique(nms)
  X <- XX

  #mean over features in featuresinf
  ############

  if (is.matrix(featuresinf)){

    selhomfeatures <- featuresinf[1, ] %in% colnames(X)
    selfeaturesinf <- featuresinf[, selhomfeatures]
    Xhom <- lapply(1:ncol(selfeaturesinf),
                   function(i){
                     if (selfeaturesinf[2,i] %in% colnames(X)){
                       out <- rowMeans(X[,selfeaturesinf[, i]], na.rm = TRUE)
                     }else{
                       out <- X[, selfeaturesinf[1,i]]
                     }
                     out
                   })
    Xhom <- do.call(cbind, Xhom)
    colnames(Xhom) <- paste(selfeaturesinf[1,],selfeaturesinf[2,], sep="-")

    #redefine expression data
    X <- Xhom
    nms <- colnames(X)
  }



  #obtain residuals on XX of features data after adjusting for covariates in teffdata
  ############

  XX <- parallel::mclapply(1:ncol(X), function(i){
    covvars <- colnames(teffdata)
    xi <- X[,i]
    fla <- paste("xi ~", paste(covvars, collapse="+"))

    lm(as.formula(fla) , data=data.frame(teffdata))$residuals
  }, mc.cores = cores)


  #redefine feature data
  X <- do.call(cbind, XX)
  colnames(X) <- nms


  #binarize feature data (profiles to compare)
  Xscale <- scale(X)
  profcomp <- lapply(1:ncol(Xscale), function(i) (Xscale[,i]>0))
  profcomp <- do.call(cbind, profcomp)
  colnames(profcomp) <- colnames(X)

  #compare the binarized feature data with the reference profiles (prof)
  #and compute the highest level of matching (80%) useing targetprofile()

  if(effect=="positive" | effect=="positiveandnegative"){
    pref <- pteffObject$profile$profpositive
    pref <- matrix(pref[,colnames(profcomp)], nrow=nrow(pref))
    pf1 <- targetprofile(binfeatures=profcomp, profileRef=pref, match=match)
    names(pf1) <- colnames(features)
  }

  if(effect=="negative" | effect=="positiveandnegative"){
    pref <- pteffObject$profile$profnegative
    pref <- matrix(pref[,colnames(profcomp)], nrow=nrow(pref))
    pf2 <- targetprofile(binfeatures=profcomp, profileRef= pref,  match=match)
    names(pf2) <- colnames(features)
  }

  if(effect=="positiveandnegative"){
    pf3 <- pf1
    pf3[pf2==1]<- -1
    names(pf3) <- colnames(features)
  }

  if(effect=="positive") pf <-as.numeric(pf1)
  if(effect=="negative") pf <-as.numeric(pf2)
  if(effect=="positiveandnegative") pf <-as.numeric(pf3)

  if(plot==TRUE){
    #Produce plots of targeted individuals
    if(effect=="positive"){
      im1 <- rbind(profcomp[pf1,], profcomp[!pf1,])

      raster::plot(raster::raster(im1), axes = FALSE, box=FALSE , ylab = "",xlab="", col = c( "green4", "red4"), legend = FALSE)

      whchprof <- sum(!pf1)/length(pf1)
      lines(c(0,1), c(whchprof,whchprof), col="black", lwd=2)

      if(is.null(lb))
        lb <- colnames(Xscale)

      ll <- ncol(im1)
      axis(1, at = seq(1/ll/2,1-1/ll/2,length=ll), labels = lb[1:ncol(Xscale)], cex.axis = 0.7, las = 2)
      axis(2, at =  c((1+whchprof)/2, (whchprof)/2)  , labels = c("Positive", "Neutral"), pos=0)    }

    if(effect=="negative"){
      im1 <- rbind(profcomp[pf2,], profcomp[!pf2,])

      plot(raster::raster(im1), axes = FALSE, box=FALSE , ylab = "",xlab="", col = c("red4", "green4"), legend = FALSE)

      whchprof <- sum(!pf2)/length(pf2)
      lines(c(0,1), c(whchprof,whchprof), col="black", lwd=2)

      if(is.null(lb))
        lb <- colnames(Xscale)

      ll <- ncol(im1)
      axis(1, at = seq(1/ll/2,1-1/ll/2,length=ll), labels = lb[1:ncol(Xscale)], cex.axis = 0.7, las = 2)
      axis(2, at =  c((1+whchprof)/2, (whchprof)/2)  , labels = c("Negative", "Neutral"), pos=0)
    }

    if(effect=="positiveandnegative"){
      im1 <- rbind(profcomp[pf3==-1,],profcomp[pf3==0,],profcomp[pf3==1,] )

      plot(raster::raster(im1), axes = FALSE, box=FALSE , ylab = "",xlab="", col = c("red4", "green4"), legend = FALSE)

      pff <- pf3==-1
      whchprof1 <- sum(!pff)/length(pff)
      lines(c(0,1),c(whchprof1,whchprof1),col="black", lwd=2)

      pff <- pf3==1
      whchprof2 <- 1-sum(!pff)/length(pff)
      lines(c(0,1),c(whchprof2,whchprof2),col="black", lwd=2)


      if(is.null(lb))
       lb <- colnames(Xscale)

      ll <- ncol(im1)
      axis(1, at = seq(1/ll/2,1-1/ll/2,length=ll) , labels = lb[1:ncol(Xscale)], cex.axis=0.6, las=2)
      axis(2, at = c((whchprof1+1)/2,whchprof2/2), labels = c("Positive", "Negative"), pos=0)
    }

  }

  names(pf) <- colnames(features)
  res <- list(classification=pf, summary.model=NULL, effect=effect, model=NULL, teffdata=NULL)


  ##test the profile-treatment interaction on the effect using the specified model
  if(length(model)!=0)
  {
    Y <- teffdata[,"eff"]
    W <- as.numeric(as.factor(teffdata[,"t"]))==2
    dat <- data.frame(Y=Y, W=W, pf=pf)
    fla <- as.formula("Y ~ W*pf")

    #add covariates specified by nmcov
    if(length(nmcov)>0){
      dat <- data.frame(dat, data.frame(teffdata)[nmcov])
      fla <- as.formula(paste("Y ~ W*pf", paste(nmcov, collapse="+"), sep="+"))
    }

    if(model=="beta"){
      pdf(plot.name)
      Y[Y==0] <- 1e-5
      Y[Y==1] <- 1-1e-5
      m1 <- summary(betareg(fla, data=dat))
    }else{
      if(model=="hazard")
      {
        time <-  x$teffdata[,"time"]
        event <- x$teffdata[,"event"]

        dat <- data.frame(dat, time=time, event=event)

        fla <- as.formula("Surv(eff, event) ~ W*pf")
        if(length(nmcov)>0)
          fla <- as.formula(paste("Surv(time, event) ~ W*pf", paste(nmcov, collapse = "+"), sep="+"))

        m1 <- coxph(fla, data = dat)

      }else{
        if(model=="log2"){
          Y <- log2(Y+1)
          m1 <- summary(glm(fla, data=dat))
        }else{
          m1 <- summary(glm(fla, data=dat, family=model))
        }
      }
    }
    names(pf) <- colnames(features)
    res <- list(classification=pf, summary.model=m1, effect=effect, model=model, teffdata=data.frame(t=W, eff=Y))
  }

  #return profiles and results of fitted model
  attr(res, "class") <- "tarteff"

  return(res)
}

#' Match subjects to profiles
#'
#' @keywords internal
#' @param binfeatures binarized feature data with respect to the mean value over
#' the subjects in the feature data to be targeted
#' @param  profileRef reference profile of positive or negative treatment effects as obtained by
#' \link[teff]{profile}
#' @param match \code{numeric} value between 0 and 1 indicating the level of the match for the targeing
#' across all the binarized features in the profile. A value of 1 means that all the binarized
#' features must take identical values in the profile (i,e "positive") for an individual to be clsssified
#' in its subpopulation (of positive treatment effects).
#' @return a logical vector indicating the features that matched the profile

targetprofile <- function(binfeatures,
                          profileRef,
                          match = 0.8){

  res <- sapply(1:nrow(binfeatures),
                function(j){
                  map <- sapply(1:nrow(profileRef), function(i)
                  mean(binfeatures[j,] == profileRef[i,]))
                  wm <- which(map == max(map))[1]
                  out <- map[wm]
                  names(out) <- wm
                  out
                })
  #determine is the matching is greater than mtch
  res <- res > match
  return(res)
}
