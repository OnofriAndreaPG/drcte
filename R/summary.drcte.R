summary.drcte <- function(object, od = FALSE, pool = TRUE,
                          units = NULL, ...)
{
  if(object$fit$method == "KDE"){
    parVec <- as.vector(coef(object))
    varMat <- NULL
    resVar <- NULL
    if (!is.null(resVar))
    {
        varMat.us <- varMat / (2*resVar)
    } else {
        varMat.us <- NULL
    }

    parNames <- object$"parNames"[[1]]
    resultMat <- matrix(NA, length(parVec), 1,
    dimnames = list(parNames, c("Estimate")))

    resultMat[, 1] <- parVec

    fctName <- deparse(object$call$fct)
    sumObj <- list(NA, varMat, resultMat, object$"boxcox", fctName, object$"robust", NULL, object$"type",
    df.residual(object), varMat.us, object$"fct"$"text", object$"fct"$"noParm", NULL)
    names(sumObj) <- c("resVar", "varMat", "coefficients", "boxcox", "fctName", "robust", "varParm", "type",
    "df.residual", "cov.unscaled", "text", "noParm", "rseMat")
    sumObj$robust <- "no"

  }  else if(object$fit$method == "NPMLE")
    {
    parVec <- NULL
    varMat <- NULL
    resVar <- NULL
    varMat.us <- NULL

    resultMat <- do.call(rbind, lapply(object$fit$icfitObj, function(x) x[,-c(1,2)]))
    fctName <- deparse(object$call$fct)
    # Naive SEs
    resultMat$Naive.SE <- sqrt(resultMat$cdf * (1 - resultMat$cdf))/sqrt(sum(resultMat$count))
    sumObj <- list(NA, varMat, resultMat, object$"boxcox", fctName, object$"robust", NULL, object$"type",
    df.residual(object), varMat.us, object$"fct"$"text", object$"fct"$"noParm", NULL)
    names(sumObj) <- c("resVar", "varMat", "coefficients", "boxcox", "fctName", "robust", "varParm", "type",
    "df.residual", "cov.unscaled", "text", "noParm", "rseMat")
    sumObj$robust <- "no"
    } else {

    if(!is.null(units)){
      vcovNew <- vcovCL(object, cluster = units)
      retMat <- coeftest(object, vcov. = vcovNew)
      class(object) <- "drc"
      sumObj <- summary(object, od = od, pool = pool)
      sumObj$coefficients <- retMat[,]
      sumObj$varMat <- vcovNew
      sumObj$robust <- "Cluster robust sandwich SEs"
    } else {
      class(object) <- "drc"
      sumObj <- summary(object, od = od, pool = pool)
      sumObj$robust <- "no"
    }
    }
  sumObj$resVar <- NULL

  class(sumObj) <- c("summary.drcte", "summary.drc")
  return(sumObj)
}

print.summary.drcte <- function(x, ...)
{
    object <- x
    cat("\n")
    cat(paste("Model fitted: ", object$"text", "\n", sep = ""))
    cat("\n")
    cat("Robust estimation:", object$"robust", "\n")
    cat("\n")
    if(object$fctName == "KDE()"){
      cat("Bandwidth estimates:\n\n")
    } else if (object$fctName == "NPMLE()"){
      cat("Turnbull's intervals and masses:\n\n")
    } else {
      cat("Parameter estimates:\n\n")
    }
    printCoefmat(object$"coefficients")
    invisible(object)
}
