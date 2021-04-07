#' Targets individuals into groups of high, low or neutral treatment effects
#'
#' @details This function uses feature data of individuals to classify them
#' into subpopulations associated with high and low treatment effects. The
#' classification is performed by targeting the feature data (adjusted by
#' covariates and binarized) to the profiles computed in \link[teff]{profile}.
#'
#' The function tests whether the classification of the subjects stratifies
#' the association of the treatment with the effect, testing the
#' interaction between subpopulation classification and the treatment.
#'
#' Models on the effects to test the interaction include the general lineal
#' models, beta regression and proportional hazards models.
#'
#' The function can be used to target new individuals not used in the profiling
#' and/or other types of effects expected from the treatment.
#'
#' @export
#' @param  x a \code{list} with a fields \code{teffdata} and \code{features}.
#' \code{teffdata} the data for treatment, effect and covariates across
#' subjects. \code{features} data is the data on which the profiling
#' is done. Adjusted features for the variables on the \code{teffdata} are
#' used to fit the forest and extract the profiles of individuals with
#' significant  treatments.
#' @param pteffObject object of class \code{pteff}
#' @param effect \code{character} with one of the values: "high", "low", "highandlow", indicating
#' which subpopulation to be targeted: subpopulation with high treatment effect only, low
#' treatment effect only or both, respectively (Default: high).
#' @param featuresinf a \code{vector} of characters with the names of the
#' features to be used; a \code{matrix} of characters whose columns are
#' names of features whose values will be averaged (Default:NULL).
#' @param plot \code{logical} indicating a whether a plot of the targeting will be produced
#' @param lb \code{character} indicating labels for features plot.
#' @param model \code{character} with one of the values corresponding to the family argument of a
#' general linear model (\link[stats]{glm}), such as: "gaussian", "binomial", etc. It also allows "beta" for beta
#' regression, "hazard" for proportional hazards model \link[survival]{coxph}, and "log2" for log2
#' transformation of the effects (Default: "gaussian").
#' @param match \code{numeric} value between 0 and 1 indicating the level of the match for the targeing
#' across all the binarized features in the profile. A value of 1 means that all the binarized
#' features must take identical values in the profile (i,e "high") for an individual to be clsssified
#' in its subpopulation (of high treatment effects).
#' @param nmcov vector of \code{character} with the names of the covariates to be included in the models.
#' @param cores an integer with the number of cores for parallel
#' computation (Default:1).
#' @param ... additional parameters for the \code{image} function
#' @return   a \code{list} of class \code{tarteff} with fields:
#' \describe{
#' \item{classification:}{classification of individuals into subpopulations of
#' expected high and/or low treatment effects}
#' \item{summary.model:}{a \code{summary} of the model used to test the interaction between
#' the profile ("highandlow": high=1, neutral=0, low=-1; "high": high=1, neutral=0,
#' "low": low=1, neutral=0) with the effect.}
#' }
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous, profile=TRUE)
#' res <- target(tcell, pf, effect="highandlow", featuresinf=homologous, nmcov="age", model="log2")
#' res
#'
target <- function(x,
                   pteffObject,
                   effect="high",
                   featuresinf=NULL,
                   plot=TRUE,
                   lb=NULL,
                   model="gaussian",
                   match=0.8, nmcov=NULL, cores=1, ...){
  teffdata <- x$teffdata
  features <- t(x$features)

  if(effect=="high"){
    nmfeatures <- unlist(strsplit(colnames(pteffObject$profile$profhigh), "-"))
  }else{
    nmfeatures <- unlist(strsplit(colnames(pteffObject$profile$proflow), "-"))
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

  if(effect=="high" | effect=="highandlow"){
    pref <- pteffObject$profile$profhigh
    pref <- matrix(pref[,colnames(profcomp)], nrow=nrow(pref))
    pf1 <- targetprofile(binfeatures=profcomp, profileRef=pref, match=match)
    names(pf1) <- colnames(features)
  }

  if(effect=="low" | effect=="highandlow"){
    pref <- pteffObject$profile$proflow
    pref <- matrix(pref[,colnames(profcomp)], nrow=nrow(pref))
    pf2 <- targetprofile(binfeatures=profcomp, profileRef= pref,  match=match)
    names(pf2) <- colnames(features)
  }

  if(effect=="highandlow"){
    pf3 <- pf1
    pf3[pf2==1]<- -1
    names(pf3) <- colnames(features)
  }

  if(effect=="high") pf <-as.numeric(pf1)
  if(effect=="low") pf <-as.numeric(pf2)
  if(effect=="highandlow") pf <-as.numeric(pf3)

  if(plot==TRUE){
    #Produce plots of targeted individuals
    if(effect=="high"){
      im1 <- t(rbind(profcomp[!pf1,],profcomp[pf1,]))
      image(im1, axes = FALSE, ylab = "Subjects",
          xlab="", col = gray.colors(2))
      whchprof <- sum(!pf1)/length(pf1)
      lines(c(-1,100), c(whchprof,whchprof), col="red", lwd=2)

      if(is.null(lb))
        lb <- colnames(Xscale)

      axis(1, at = seq(0,1,length=nrow(im1)) , labels = lb[1:ncol(Xscale)], cex.axis = 0.7, las = 2)
      axis(2, at =  (1+whchprof)/2, labels = c("Positive"))
    }

    if(effect=="low"){
      im1 <- t(rbind(profcomp[!pf2,], profcomp[pf2,]))
      image(im1, axes = FALSE, ylab = "Subjects",
          xlab="", col = gray.colors(2))
      whchprof <- sum(!pf2)/length(pf2)
      lines(c(-1,100), c(whchprof,whchprof), col="red", lwd=2)

      if(is.null(lb))
        lb <- colnames(Xscale)

      axis(1, at = seq(0,1,length=nrow(im1)), labels = lb[1:ncol(Xscale)], cex.axis = 0.7, las = 2)
      axis(2, at =  (1+whchprof)/2  , labels = c("Negative"))
    }

    if(effect=="highandlow"){
      im1 <- t(rbind(profcomp[pf3==-1,],profcomp[pf3==0,],profcomp[pf3==1,] ))
      image(im1, axes = FALSE, ylab = "Subjects",
          xlab="", col = gray.colors(2), ...)
      pff <- pf3==1
      whchprof1 <- sum(!pff)/length(pff)
      lines(c(-1,100),c(whchprof1,whchprof1),col="red", lwd=2)

      pff <- pf3==-1
      whchprof2 <- 1-sum(!pff)/length(pff)
      lines(c(-1,100),c(whchprof2,whchprof2),col="green", lwd=2)


      if(is.null(lb))
       lb <- colnames(Xscale)

      axis(1, at = seq(0,1,length=nrow(im1)) , labels = lb[1:ncol(Xscale)], cex.axis=0.7, las=2)
      axis(2, at = c((whchprof1+1)/2,whchprof2/2), labels = c("Positive", "Negative"))
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
#' @param  profileRef reference profile of high or low treatment effects as obtained by
#' \link[teff]{profile}
#' @param match \code{numeric} value between 0 and 1 indicating the level of the match for the targeing
#' across all the binarized features in the profile. A value of 1 means that all the binarized
#' features must take identical values in the profile (i,e "high") for an individual to be clsssified
#' in its subpopulation (of high treatment effects).
#' @return a logical vector indivcating the featrues that matched the profile

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
