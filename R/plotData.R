"plotData" <- function(x, xlim, confidence.level = 0.95, gridsize = 100,
                     type = c("average", "all", "bars", "none", "obs", "confidence"),
                     npmle.type = c("interpolation", "midpoint", "right", "left", "none")
                     ){
  # This function is used to explore 'drc' and 'drcte' objects
  # and extract the data to be used for ggplots
  # Date of editing: 19/01/2022

  object <- x
  type <- match.arg(type)
  npmle.type <- match.arg(npmle.type)
  method <- object$fit$method
  dataList <- object[["dataList"]]
  dlNames <- dataList[["names"]]
  doseName <- dlNames[["dName"]]
  respName <- dlNames[["orName"]]
  curveidName <- dlNames[["cNames"]]

  ## Determining logarithmic scales
  logX <- FALSE

  ## Constructing the observed data
  if(method == "NPMLE"){
    # To display the NPMLE (pkg. interval)
    # To be completed with several methods
    dose <- object$ICfit$npmle$time
    resp <- object$ICfit$npmle$cdf
    curveid <- object$ICfit$npmle[,1]
    plotid <- as.numeric(object$ICfit$npmle[,1])
  } else {
    # To display naive end-point estimator (upd. 21/12/21)
    dose <- dataList[["dose"]]
    resp <- dataList[["origResp"]]
    curveid <- dlNames$rNames[dataList$curveid]
    plotid <- dataList[["plotid"]]
  }

  # Capire se voglio plottare tutti i dati...
  plotPoints <- data.frame(dose, resp, curveid, plotid)
  # print(plotPoints)
  # print(c(doseName, "CDF", curveidName))
  if(is.null(curveidName)) curveidName <- "curveid"
  names(plotPoints)[1:(length(doseName) + 2)] <- c(doseName, "CDF", curveidName)

  if (!is.null(plotid))
  {  # used for event-time data
      assayNoOld <- as.vector(plotid)
  } else {
      assayNoOld <- as.vector(curveid)
  }
  uniAss <- unique(assayNoOld)
  numAss <- length(uniAss)

  # Get the prediction function
  plotFct <- (object$"curve")[[1]]
  if(method != "NPMLE" && method != "KDE"){
    logDose <- (object$"curve")[[2]]
    naPlot <- ifelse(is.null(object$"curve"$"naPlot"), FALSE, TRUE)
  }

  # Handling multiple covariates
  doseDim <- 1
  if(!is.vector(dose)){
    doseDim <- length(dose[1,])
    doseOld <- dose
    addVars <- dose[,-1]
    dose <- dose[,1]
  }

  ## Determining range of dose values
  if (missing(xlim))
  {
    if(method != "NPMLE" && method != "KDE")  xLimits <- c(min(dose), max(dose)) else xLimits <- c(0, max(dose))
  } else {
    xLimits <- xlim  # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
  }



  ## Handling small dose values
  conLevel <- round(min(dose[is.finite(dose)])) - 1

  if ((xLimits[1] < conLevel) && (logX || (!is.null(logDose))))
  {
      xLimits[1] <- conLevel
      smallDoses <- (dose < conLevel)
      dose[smallDoses] <- conLevel
      if (is.null(conName))
      {
          if (is.null(logDose)) {conName <- expression(0)} else {conName <- expression(-infinity)}
      }
  } else {
      conName <- NULL
  }
  if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}

  ## Constructing dose values for plotting
   if (doseDim == 1)
     {
     dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
     } else {
     dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
     addVars <- as.numeric(apply(as.data.frame(addVars), 2, function(col) unique(as.character(col))))
     dosePts <- expand.grid(dosePts, addVars)
     }
  if(method == "Parametric"){
    plotMat <- plotFct(dosePts)
    } else if(method == "KDE"){
    plotMat <- lapply(plotFct, function(x) x(dosePts))
    } else {
    plotMat <- NULL
  }

  plotMat <- as.data.frame(plotMat)
   if(method == "NPMLE"){
     retData <- NULL
   } else {
     doseName <- dlNames[["dName"]]
     if(length(dlNames[["rNames"]]) > 1){
       respName <- as.character(dlNames[["rNames"]])
     } else {
       respName <- "CDF"
     }
     retData <- data.frame(dosePts, as.data.frame(plotMat))
     colnames(retData) <- c(doseName, respName)
   }
   # print(retData)
   # Melt plot data
  if(method != "NPMLE" && numAss > 1){
    retData <- tidyr::pivot_longer(retData, names_to = dlNames$cNames,
                          values_to = "CDF",
                          cols = c(2:length(retData[1,])))
    retData <- as.data.frame(retData)
  }

  returnList <- list(plotPoints = plotPoints, plotFits = retData)
  return(returnList)
  # retData <- data.frame(dosePts, as.data.frame(plotMat))
  # colnames(retData) <- c(doseName, dlNames[["cNames"]])
  # returnList <- list(plotPoints = plotPoints, plotFits = retData)
  # return(invisible(returnList))

}




