#' Plots predicted treatment effects of individuals from
#' \link[teff]{predicteff}
#'
#' @export
#' @param x object of class \code{pteff}
#' @param rk object of class \code{vector}, if null, treatment effect are
#' plotted against their ranking, if not then they ara plotted against rk values.
#' @param lb label of the y axis for treatment effect.
#' @param xlab label of the x axis.
#' @param ctrl.plot controls plot legends (when NULL then internally ctrl.plot <- list(lb=c("Not treated", "Treated"),
#' wht="bottomleft", whs = "topright"))
#' @return A plot on the current graphics device.
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous)
#' plotPredict(pf)

plotPredict <- function(x, ..., rk=NULL, lb="Associated treatment effect", xlab = "Subject Ranking", ctrl.plot = NULL){
  if(class(x)!="pteff"){
    stop("x should be of class pteff")
  }

  yrange <- c(x$cl, x$cu)

  colsigpositive <-  (x$cl>x$resplevel)
  colsignegative <- (x$cu<x$resplevel)

  colsighet <- colsigpositive+1
  colsighet[colsigpositive==1] <- 3
  colsighet[colsignegative==1] <- 3

  coltreatment <- rep("orange", length(x$treatment))
  coltreatment[x$treatment == 1] <- "blue"

  ranktau <- rank(x$predictions)

  if(!is.null(rk))
    ranktau <- rk

  graphics::plot(ranktau,  as.vector(x$predictions),
       ylim = c(min(yrange), max(yrange)), type = "p",
       pch = 16, ylab="", xlab=xlab,
       col = coltreatment, ...)


  title(ylab=lb, line=2)


  for(i in 1:length(x$predictions))
    graphics::lines(c(ranktau[i], ranktau[i]), c(x$cl[i],x$cu[i]), col = colsighet[i])

  graphics::points(ranktau, x$predictions, pch = 16,col = coltreatment)

  graphics::lines(c(-10,500), c(x$resplevel,x$resplevel), lwd=1.5, lty=2, col="red")

  if(is.null(ctrl.plot)){
    ctrl.plot <- list(lb=c("Not treated", "Treated"),
                      wht="bottomleft", whs = "topright")
  }

  graphics::legend(ctrl.plot$wht, legend=ctrl.plot$lb, pch=16, col=c("orange","blue"), bty="n" )
  graphics::legend(ctrl.plot$whs, legend=c("significant"), lty=1, col=3, bty="n" )

}

#' Prints pteff object
#'
#' @param x object of class \code{pteff}
#' @return object of class \code{pteff}
#' @export

print.pteff <- function(x){

  p1 <- x$predictions
  cat("object of class: pteff \n")
  cat("Estimated treatment effects in $predictions: \n")
  print(p1)
}

#' Plots estimated treatment effects of individuals from
#' \link[teff]{target}
#'
#' @export
#' @param x object of class \code{tarteff}
#' @param labs string of characters for the labels of the plot, it refers in order to labels to use for:
#'  Treatment effect group, Outcome, Treatment, and levels of the treatment like: Treated and not Treated.
#' @param labeff string of characters for the labels of the treatment effect default NULL, 0: neutral and 1: positive or negative; or 0: neutral, -1: negative and 1: positive.
#' @return A plot on the current graphics device
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous)
#' res <- target(tcell, pf, effect="positiveandnegative", featuresinf=homologous, nmcov="age", model="log2")
#' plotTarget(res)

plotTarget <- function(x, ..., labs=c("Treatment effect", "Outcome", "Treatment", "Not treated", "Treated"), labeff=NULL){

  if(length(x$model)==0){
    stop("not available plot: no iteraction model was fitted")
    return(NULL)
  }

  t <- factor(x$teffdata[,"t"], labels = labs[4:5])

  dd <- data.frame(eff=x$teffdata[,"eff"], t=t, pf=x$classification)
  dd <- dd[complete.cases(dd),]

  names(dd)[1] <- labs[2]
  names(dd)[2] <- labs[3]
  names(dd)[3] <- labs[1]

  if(!is.null(labeff))
    dd[[3]] <- factor(dd[[3]], labels = labeff)

  ggpubr::ggline(dd, x = labs[1], y = labs[2],
           add = "mean_ci", color = labs[3], palette = c("orange", "blue"),
           xlab=labs[1], main="", ylab=labs[2])
}

#' Box plots for inferred treatment effects of individuals from
#' \link[teff]{target}
#'
#' @export
#' @param x object of class \code{tarteff}
#' @param labs string of characters for the labels of the plot, it refers in order to labels to use for:
#' Treatment effect, Outcome, Treatment, and levels of the treatment like: Treated and not Treated.
#' @param lg position of the legend for treatment levels such as "topright"
#' @param labeff string of characters for the labels of the treatment effect default NULL, 0: neutral and 1: positive or negative; or 0: neutral, -1: negative and 1: positive.
#' @return A plot on the current graphics device#'
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous)
#' res <- target(tcell, pf, effect="positiveandnegative", featuresinf=homologous, nmcov="age", model="log2")
#' boxplot(res, lg="topright")

boxPlot <- function(x, labs=c("Treatment effect", "Outcome", "Treatment", "Not treated", "Treated"), lg=NULL, labeff=NULL){

    if(length(x$model)==0){
      stop("not available plot: no iteraction model was fitted")
      return(NULL)
    }

    t <- factor(x$teffdata[,"t"])

    dd <- data.frame(eff=x$teffdata[,"eff"], t=t, pf=x$classification)
    dd <- dd[complete.cases(dd),]

    names(dd)[1] <- labs[2]
    names(dd)[2] <- labs[3]
    names(dd)[3] <- labs[1]

    if(!is.null(labeff))
      dd[[3]] <- factor(dd[[3]], labels = labeff)

    fc <- factor(paste(dd[,2], as.factor(dd[,3]),  sep="-"))
    boxplot(dd[,1] ~ fc, col=rep(c("orange", "blue"), each=length(levels(factor(dd[,3])))), ylab=labs[2], xlab=labs[1], xaxt="n", yaxt="n", main="", cex.lab=1.4, cex.axis=1.3, cex.main=1.4)

    axis(1,at=1:(2*(length(levels(factor(dd[,3]))))),labels=rep(levels(factor(dd[,3])), 2), cex.axis=1.2, cex.lab=1.4)
    axis(2,cex.axis=1.4)

    if(length(lg)!=0)
      legend(lg, legend=labs[4:5], col=c("orange", "blue"), pch=16)
}


#' Prints tarteff object
#'
#' @param x object of class \code{tarteff}
#' @return object of class \code{tarteff}
#' @export

print.tarteff <- function(x){
  mt <- match(x$effect, c("positive", "positiveandnegative", "negative"))

  if(mt==1) lb <- "  positive treatment effect: 1\n"
  if(mt==2) lb <- "  negative treatment effect: -1\n  neutral: 0\n  positive treatment: 1\n"
  if(mt==3) lb <- "  negative treatment effect: 1\n  neutral: 0\n"

  cat("object of class: tarteff \n")
  cat("\n")
  cat("classification into \n")
  cat(lb)
  cat(" ")
  print(table(x$classification))

  if(length(x$model)!=0)
  {
    cat("\n")
    cat("interaction fitted model: ")
    cat(x$model)
    cat("\n")
    print(x$summary.model$coeff["WTRUE:pf",])
  }
}


#' returns the \code{mean} and the error limits defined by the
#'   \code{confidence interval}.
#'
#' @param x vector
#' @param ci confidence limit
#' @param character
#' @return vector
#' @export

mean_ci <- function(x, ci = 0.95, error.limit = "both"){
  length <- base::sum(!is.na(x))
  sd = stats::sd(x, na.rm=TRUE)
  se <- sd / sqrt(length)
  .mean <- base::mean(x, na.rm = TRUE)
  ci <- stats::qt(ci/2 + .5, length-1)*se
  data.frame(
    y =  .mean,
    ymin = .mean - ci,
    ymax = .mean + ci
  )
}

