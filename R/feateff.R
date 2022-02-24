#' Extracts treatment, the outcome on which the effect is assessed, and feature data from transcriptomic and methylomic studies
#'
#'
#' @details This function extracts feature and treatment-effects data,
#' from  \code{eSet} or \code{SummarizedExperiment} objects
#' for profiling with \link[teff]{predicteff}. The function includes the option of adding
#' surrogate variables as additional covariates to the treatment-effect data.
#' @export
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or
#' \code{SummarizedExperiment}
#' @param tname \code{character} with name of treatment variable
#' @param effname \code{character} with name of outcome variable on which the effect is measured (time to event in the case of survival)
#' @param reft \code{character} vector that indicates no-treatment and treatment levels (Default: NULL)
#' @param refeff \code{character} vector that indicates no-outcome and outcome levels when the outcome on which the effect is measured is categorical (Default: NULL)
#' @param event \code{character} with name of event variable in the case of a survival outcome (Default: NULL)
#' @param covnames \code{character} vector with names of covariates (Default: NULL)
#' @param covtype \code{character} vector with character "n" indicates which covariates are numerical e.g. c("n", "n", "c") (Default: NULL)
#' @param sva \code{logical} indicates whether surrogate variable should be added as covariates (Default: NULL)
#' @param betas  \code{logical} indicates whether beta values be
#' used if \code{set} is a \code{GenomicRatioSet} (Default: TRUE)
#' @param UsegeneSymbol \code{logical} indicates whether genenames should be used as feature names (Default: FALSE)
#' @param rnaseq \code{logical} indicates if expression data is RNA-seq (TRUE) or microarray (FALSE, default)
#'
#' @examples
#'
#'\donttest{
#'#library(GEOquery)
#'#gsm <- getGEO("GSE17755")
#'#gsm <- gsm[[1]]
#'
#'#in this example we use sex as the treatment variable do detect groups of
#'# high sexual dimorphism in arthritis disease.
#'
#'#data4teff <- feateff(gsm, tname="gender:ch1", reft=c("male", "female"),
#'#                      effname="disease:ch1", refeff=c("healthy","arthritis"),
#'#                      covnames="age:ch1", covtype="n",
#'#                      sva=TRUE, UsegeneSymbol=TRUE)
#'
#'}


feateff <- function(set,
                    tname,
                    effname,
                    reft=NULL,
                    refeff=NULL,
                    event=NULL,
                    covnames = NULL,
                    covtype = NULL,
                    sva = FALSE,
                    betas = TRUE,
                    UsegeneSymbol=FALSE,
                    rnaseq=FALSE){


  ##get teffdata
  if (is(set, "eSet")){
    pFun <- Biobase::pData
  }
  else if (is(set, "SummarizedExperiment")){
    pFun <- SummarizedExperiment::colData
  } else{
    stop("set must be an eSet or a SummarizedExperiment derived object")
  }


  phenos <- pFun(set)

  message("effect and treatment variables converted to numeric")
  eff <- as.numeric(phenos[,effname])

  if(!is(refeff, "NULL")){
    eff[grep(refeff[1], phenos[,effname])] <- 0
    eff[grep(refeff[2], phenos[,effname])] <- 1
    eff <- as.numeric(eff)
  }

  t <- as.numeric(phenos[,tname])

  if(!is(refeff, "NULL")){
    t[grep(reft[1], phenos[,tname])] <- 1
    t[grep(reft[2], phenos[,tname])] <- 2
    t <- as.numeric(t)
  }


  wnumcov <- covnames[covtype=="n"]

  numcov <- lapply(phenos[wnumcov],
                 function(i) as.numeric(sapply(strsplit(i, ":"), function(j) j[[length(j)]]))
                 )

  numcov <- do.call(cbind,numcov)

  wnonumcov <- covnames[covtype!="n"]
  cov <- cbind(numcov,phenos[wnonumcov])

  teffdata <- data.frame(eff, t, cov)


  #get features data
  if (is(set, "ExpressionSet")){
    features <- Biobase::exprs(set)
    if(UsegeneSymbol==TRUE){
      genesIDs <- fData(set)
      rownames(features) <- genesIDs$"Gene symbol"
    }
  } else if (is(set, "SummarizedExperiment")){
    features <- SummarizedExperiment::assay(set)
  } else {
    stop("set must be an ExpressionSet, MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }


  #get features and teff data for complete.cases only
  ss <- sapply(1:nrow(features), function(x) sum(is.na(features[x,])))
  features <- features[ss==0, ]

  selcom <- complete.cases(teffdata) & complete.cases(t(features))
  teffdata <- teffdata[selcom, ]
  features <- features[, selcom]


  #complete covariates with SVAs
  if(sva==TRUE)
  {
    fla <- formula(paste0("~", paste0(names(teffdata), collapse = "+"), collapse=""))
    mod0 <- model.matrix(fla, data = teffdata)

    fla <- formula(paste0("~ eff:t+", paste0(names(teffdata), collapse = "+"), collapse=""))
    mod <- model.matrix(fla, data = teffdata)

    ns <- sva::num.sv(features, mod, method="be")

    if(rnaseq){
      ss <- sva::svaseq(features, mod, mod0, n.sv=ns)$sv
      v <- limma::voom(features, design = mod)
      features <- v$E

    }else{
      ss <- sva::sva(features, mod, mod0, n.sv=ns)$sv
    }

    colnames(ss) <- paste("cov", 1:ncol(ss), sep="")

    teffdata <- cbind(teffdata, ss)
  }

  res<-list(features=t(features), teffdata=teffdata)
  return(res)

}


