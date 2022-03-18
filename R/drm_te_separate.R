drmte_sep_p <- function(formula, curveid, data, subset, fct,
                      start, na.action, control, lowerl,
                      upperl) {
  # We have two possibilities for separate fitting: one is for
  # models with a 'cured' fraction and the other one for all other models
  # Updated: 17/03/22
  fctName <- deparse(substitute(fct))
  grepl(fctName, "LL.3(", fixed=TRUE)
  if(grepl(fctName, "LL.3(", fixed=TRUE) |
     grepl(fctName, "L.3(", fixed=TRUE) |
     grepl(fctName, "LN.3(", fixed=TRUE) |
     grepl(fctName, "W1.3(", fixed=TRUE) |
     grepl(fctName, "W2.3(", fixed=TRUE) |
     grepl(fctName, "G.3(", fixed=TRUE) ){
    result <- drmte_sep1(formula = formula,
                               data = data, subset = subset,
                               fct = fct, start = start, na.action = na.action,
                               control = control,
                               lowerl = lowerl, upperl = upperl)
    return(result)
  } else {
    result <- drmte_sep2(formula = formula,
                               data = data, subset = subset,
                               fct = fct, start = start, na.action = na.action,
                               control = control,
                               lowerl = lowerl, upperl = upperl)
    return(result)
  }
  }

drmte_sep <- function(formula, data, subset, fct,
                      start, na.action, control, lowerl,
                      upperl) {
  # This function fits a time-to-event model. If the attempt fails,
  # a simpler model is fitted, where only the fraction of individuals
  # with event is estimated. assuming no time course of events
  # Updated: 12/01/22
  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")

  # Find total number of individuals per curve
  nTot <- sum(nSeeds)
  nMax <- sum(nSeeds[timeAf!=Inf])
  pMax0 <- nMax/nTot
  nFirst <- sum( nSeeds[timeBef == min(timeBef)] )
  pFirst <- nFirst/nTot
  tFirst <- timeAf[timeBef==min(timeBef)]
  tLast <- timeBef[!is.finite(timeAf)]

  pCum <- cumsum(tapply(nSeeds[is.finite(timeAf)], timeAf[is.finite(timeAf)], sum))/nTot
  tpCum <- as.numeric(names(pCum))

  # Fitting models
  cureMod <- try( drmte(nSeeds ~ timeBef + timeAf, fct = fct), silent = T)
  fct$fixed <- c(Inf, NA, 1)
  fct$noParm <- 1
  cureMod2 <- drm(pCum ~ tpCum, fct = EV.fr(), start = mean(pCum))

  # Deciding which model is to be used
  if( all(class(cureMod) != "try-error")) {
    p <- try( summary(cureMod), silent=T)
    if(any(class(p) == "try-error")) {
      class(cureMod) <- "try-error"
    } else if (coef(cureMod)[2] > 1 | coef(cureMod)[2] < 0) {
      cureMod <- drmte(nSeeds ~ timeBef + timeAf, fct = fct, upperl = c(NA, 1, NA))
    }
    }

  # Preparing and returning the results
  if(all(class(cureMod) != "try-error")){
    result <- cureMod
  } else {
    result <- cureMod2
    result$data <- data.frame(timeBef, timeAf, nSeeds,
                              curveid = 1, orig.curveid = 1, weights = 1)
    if(result$fit$hessian == 0) result$fit$hessian <- NaN
  }
  result
}

drmte_sep2 <- function(formula, data, subset, fct,
                      start, na.action, control, lowerl,
                      upperl) {
  # This function fits a time-to-event model separately for each curveid
  # level. It is similar to drmte_sep, but no attempt to fit a simpler model
  # is made. Simply, the model is fit the way it is and, if no convergence is
  # obtained, an error message is returned
  # Updated: 17/03/22
  callDetail <- match.call()

  fitMod <-drmte(formula = formula, fct = fct,
                 data = data, subset = subset) #,
                               # start = start, na.action = na.action,
                               # control = control,
                               # lowerl = lowerl, upperl = upperl)
  print(fitMod)
  # Preparing and returning the results
  if(fitMod != "try-error"){
    result <- fitMod
  } else {
    result <- "Model could not be fitted for this levels"
  }
}

"EV.fr" <- function(names = c("d"))
{
  ## This is an helper function, that specifies a model where only
  ## the fraction of individuals with event is estimated,
  ## while no time-course is estimated
  numParm <- 1
  if (!is.character(names) | !(length(names) == numParm))
  {stop("Not correct 'names' argument")}
  parmVec <- rep(0, numParm)

  ## Defining the function
  fct <- function(x, parm)
  {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)

    d <- parmMat[, 1]
    d
  }

  ## Defining self starter function
  ssfct <- function(dataf)
  {
    x <- dataf[, 1]
    y <- dataf[, 2]
    d <- mean(y)
    return(d)
  }

  ## Defining names
  pnames <- names

  ## Defining derivatives
  deriv1 <- function(x, parm){
    d1 <- 0
    cbind(d1)
  }
  ## Defining the ED function

  ## Defining the inverse function

  ## Defining descriptive text
  text <- "Horizontal line (mean model)"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, deriv1 = deriv1, names = pnames, text = text)

  class(returnList) <- "drcMean"
  invisible(returnList)
}


"sepFit2obj" <- function(fitList){
  # This is a helper function, that takes a list of drcte objects
  # and prepares a unique 'drcte' object to be returned.
  lenList <- length(fitList)
  oneFunction <- lenList==1

  uniCur <- names(fitList)
  numCur <- lenList
  numPar <- as.numeric(sapply(fitList, function(el) length(coef(el))))

  retList <- list()
  # ParameterMatrix parmMat (to be done)
  # parmMat <- lapply(fitList, function(el) el$"parmMat") # Da verificare
  # do.call(cbind, parmMat)

  # data matrix
  dataMat <- lapply(fitList, function(el) el$"data") # Da verificare
  dataMat <- do.call(rbind, dataMat)
  tmp <- mapply(function(i) rep(names(fitList)[i], length(fitList[[i]]$data$curveid)), 1:numCur, SIMPLIFY = F)
  dataMat$orig.curveid <- do.call(c, tmp)
  tmp <- mapply(function(i) rep(i, length(fitList[[i]]$data$curveid)), 1:numCur, SIMPLIFY = F)
  dataMat$curveid <- do.call(c, tmp)
  retList$"data" <- dataMat

  # Parameter names
  pnList <- lapply(fitList, function(el) el$"parNames")
  npVec <- as.vector(unlist(lapply(pnList, function(x){x[[2]]})))
  idVec <- rep(uniCur, numPar)
  aVec <- paste(npVec, idVec, sep = ":")
  retList$"parNames" <- list(aVec, npVec, idVec)

  # dataList
  tdataList <- lapply(fitList, function(el) el$"dataList"$dose)
  retList$"dataList"$dose <- as.numeric(do.call(c, tdataList))
  tdataList <- lapply(fitList, function(el) el$"dataList"$origResp)
  retList$"dataList"$origResp <- as.numeric(do.call(c, tdataList))
  retList$"dataList"$weights <- NULL
  tdataList <- lapply(fitList, function(el) el$"dataList"$curveid)
  tmp <- mapply(function(i) rep(i, length(tdataList[[i]])), 1:numCur, SIMPLIFY = F)
  retList$"dataList"$curveid <- as.numeric(do.call(c, tmp))
  retList$"dataList"$plotid <- as.numeric(do.call(c, tmp))
  retList$"dataList"$resp <- retList$"dataList"$origResp
  retList$"dataList"$names <- list()
  retList$"dataList"$names$dName <- fitList[[1]]$dataList$names$dName
  retList$"dataList"$names$orName <- fitList[[1]]$dataList$names$orName
  retList$"dataList"$names$wName <- fitList[[1]]$dataList$names$wName
  # print(uniCur)
  retList$"dataList"$names$rNames <- uniCur

  # To be completed
  retList$"objVal" <- fitList
  retList$"fit"$method <- "Parametric"
  retList$"parmMat" <- NULL

  plotFct <- function(x) {matrix(unlist(lapply(fitList, function(y)y$"curve"[[1]](x))), ncol = numCur)}
  retList$"curve" <- list(plotFct, fitList[[1]]$"curve"[[2]])

  # retList$"indexMat" <- matrix(c(1:(numCur * numPar)), numPar, numCur)

  coefVec <- as.vector(unlist(lapply(fitList, function(x){x$"fit"$"par"})))
  names(coefVec) <- aVec
  retList$"coefficients" <- coefVec

  # retList$"df.residual" <- sum(unlist(lapply(fitList, function(x){x$"df.residual"})))
  # retList$"minval" <- sum(unlist(lapply(fitList, function(x){x$"fit"$"value"})))

  class(retList) <- c("drcteList", "drcte", "drc")
  return(retList)
}
