#' Plots inferred treatment effects of individuals from
#' \link[teff]{predicteff}
#'
#' @param x object of class \code{pteff}
#' @param rk object of class \code{vector}, if null treatments are
#' plotted against their ranking, if not then they ara plotted against rk values.
#' @param lb label of the y axis for treatment effect.
#' @param xlab label of the x axis.
#' @return A plot on the current graphics device.
#' @export
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous)
#' plot(pf)

plot.pteff <- function(x, rk=NULL, lb="Associated treatment effect", xlab = "Subject Ranking", ...){
  if(class(x)!="pteff"){
    stop("x should be of class pteff")
  }

  yrange <- c(x$cl, x$cu)

  colsighigh <-  (x$cl>0)
  colsiglow <- (x$cu<0)

  colsighet <- colsighigh+1
  colsighet[colsighigh==1] <- 3
  colsighet[colsiglow==1] <- 3

  coltreatment <- rep("orange", length(x$treatment))
  coltreatment[x$treatment == 1] <- "blue"

  ranktau <- rank(x$predictions)

  if(!is.null(rk))
    ranktau <- rk


  plot(ranktau, x$predictions,
       ylim = c(min(yrange), max(yrange)), type = "p",
       pch = 16, ylab="", xlab=xlab,
       col = coltreatment, ...)

  title(ylab=lb, line=2)


  for(i in 1:length(x$predictions))
    lines(c(ranktau[i], ranktau[i]), c(x$cl[i],x$cu[i]), col = colsighet[i])

  points(ranktau, x$predictions, pch = 16,col = coltreatment)

  lines(c(-10,500), c(0,0), lwd=1.5, lty=2, col="red")

  legend("topright", c("Not treated", "Treated"), pch=16, col=c("orange","blue"), bty="n" )
  legend("bottomleft", c("significant"), lty=1, col=3, bty="n" )

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

#' Plots inferred treatment effects of individuals from
#' \link[teff]{target}
#'
#' @param x object of class \code{tarteff}
#' @return A plot on the current graphics device.
#' @export
#' @examples
#' data(tcell)
#' homologous<- matrix(c("DDX3Y","DDX3X","KDM5D","KDM5C","PRKY","PRKX","RPS4Y1","RPS4X","TXLNGY", "TXLNG", "USP9Y", "USP9X", "XIST", "XIST", "TSIX", "TSIX"), nrow=2)
#' pf <- predicteff(tcell, featuresinf=homologous)
#' res <- target(tcell, pf, effect="highandlow", featuresinf=homologous, nmcov="age", model="log2")
#' plot(res, lg="topright")

plot.tarteff <- function(x, ylab="Effect", lg=NULL, ...){

  if(length(x$model)==0){
    stop("not available plot: no iteraction model was fitted")
    return(NULL)
  }


  dd <- data.frame(eff=x$teffdata[,"eff"], t=x$teffdata[,"t"], pf=x$classification)
  dd <- dd[complete.cases(dd),]

  fc <- factor(paste(dd$t, dd$pf,  sep="-"))

  boxplot(dd$eff ~ fc, col=c("grey", "grey","grey", "red","red", "red"), ylab=ylab, xlab="", xaxt="n", yaxt="n", main="", cex.lab=1.4, cex.axis=1.3, cex.main=1.4)

  axis(1,at=1:6,labels=c("Neg.", "Neutr.", "Pos.", "Neg.", "Neutr.", "Pos"), cex.axis=1.2, cex.lab=1.4)
  axis(2,cex.axis=1.4)

  if(length(lg)!=0)
   legend(lg, legend=c("Not treated", "Treated"), col=c("grey", "red"), pch=16)
}


#' Prints tarteff object
#'
#' @param x object of class \code{tarteff}
#' @return object of class \code{tarteff}
#' @export

print.tarteff <- function(x){
  mt <- match(x$effect, c("high", "highandlow", "low"))

  if(mt==1) lb <- "  high treatment effect: 1\n"
  if(mt==2) lb <- "  low treatment effect: -1\n  neutral: 0\n  high treatment: 1\n"
  if(mt==3) lb <- "  low treatment effect: 1\n  neutral: 0\n"

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
