"predict.drcte" <- function(object, newdata, se.fit = FALSE,
                            interval = FALSE, # = c("none", "confidence", "boot"),
                            level = 0.95,
                            na.action = na.pass, # od = FALSE, vcov. = vcov,
                            NPMLE.method = "interpolation",
                            units = NULL,
                            B = 200, ...)
{

  ## Checking arguments
  # object <- mod
  # newdata <- 180
  # cluster <- NULL
  cluster <- units
  # interval <- match.arg(interval) # modified as TRUE/FALSE
  respType <- object[["type"]]
  dataList <- object[["dataList"]]
  dataSet <- object[["data"]]
  ncol <- length(dataSet[1,])

  cName <- dataList[["names"]][["cNames"]]
  cLevs <- dataList[["names"]][["rNames"]]

  if (missing(newdata))
  {
    # Case 1: no newdata
    dataSet <- dataSet[is.finite(dataSet[,2]) == T, ]
    doseVec <- dataSet[,-c(1,(ncol - 3):ncol)]
    groupLevels <- dataSet[,(ncol - 2)] # as.character(dataList[["plotid"]])
  } else {
    # Case 2: newdata is given.
    if(is.vector(newdata)) newdata <- data.frame(newdata)
    if(is.vector(dataList[["dose"]]) == T){
      # Only one predictor in the model
      dName <- dataList[["names"]][["dName"]]
      if (any(names(newdata) %in% dName))
      {
        doseVec <- newdata[, dName]
      } else {
        doseVec <- newdata[, 1]
        #         warning("Dose variable not in 'newdata'")
      }
    } else {
      # More than one predictor in the model (eg: HT models)
      doseDim <- ncol(dataList[["dose"]])

      if(length(newdata[1,]) <= doseDim) stop("The number of dependent variables in newdata is not equal to the number of dependent variables in the model")
      doseVec <- newdata
    }

    # With new data, if more than one time-to-event curve is defined,
    # predictions are made for all curves
    nameLevels <- levels(factor(dataSet[,(ncol - 1)]))
    if(is.vector(doseVec)){
      groupLevels <- rep(nameLevels, each = length(doseVec))
      doseVec <- rep(doseVec, length(nameLevels))
    } else {
      groupLevels <- rep(nameLevels, each = length(doseVec[,1]))
      doseVec <- doseVec[rep(seq_len(nrow(doseVec)), length(nameLevels)),]
      # doseVec <- rep(doseVec, length(nameLevels))
    }

  }

  # counts the number of predictions to be derived
  noNewData <- length(groupLevels)

  ## Retrieving matrix of parameter estimates
  parmMat <- object[["parmMat"]]
  if(!is.null(parmMat)) pm <- t(parmMat[, groupLevels, drop = FALSE])

  ## Retrieving variance-covariance matrix for parametric fits
  if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){
    sumObj <- summary(object, od = od) #summary(object, od = od)

    if(is.null(cluster)){
      vcovMat <- vcov(object)
    } else {
      vcovMat <- vcovCL(object, cluster = cluster)
    }

  }

  ## Defining index matrix for parameter estimates
  indexMat <- object[["indexMat"]]

  ## Calculating predicted values
  retMat <- matrix(0, noNewData, 4)
  colnames(retMat) <- c("Prediction", "SE", "Lower", "Upper")
  objFct <- object[["fct"]]

  if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){
    retMat[, 1] <- objFct$"fct"(doseVec, pm)
  } else if (object$fit$method == "KDE"){
    # KDE
    listaFun <- function(x, y) object$curve[[1]][[x]](y)
    retMat[, 1] <- as.numeric(mapply(listaFun, groupLevels, doseVec))
  } else {
    # NPMLE
    # NPMLE.method = "interpolation"
    listaFun <- function(x, y) object$curve[[1]][[x]](y, NPMLE.method)
    retMat[, 1] <- as.numeric(mapply(listaFun, groupLevels, doseVec))
    }

  ## Checking if derivatives are available
  deriv1 <- objFct$"deriv1"
  if (is.null(deriv1) & object$fit$method != "KDE" & object$fit$method != "NPMLE")
  {
    return(retMat[, 1])
  }

  ## Calculating the quantile to be used in the confidence intervals
  # if (!identical(interval, "none")) {
  if (interval == TRUE) {
    tquan <- qnorm(1 - (1 - level)/2)
  }

  ## Calculating standard errors and/or confidence intervals
  if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){
    # if (se.fit || (!identical(interval, "none")))
    if (se.fit || interval == TRUE)
    {
      sumObjRV <- 0
      piMat <- indexMat[, groupLevels, drop = FALSE]
      for (rowIndex in 1:noNewData)
      {
        parmInd <- piMat[, rowIndex]
        varCov <- vcovMat[parmInd, parmInd]

        if(is.vector(doseVec)){
          dfEval <- deriv1(doseVec[rowIndex], pm[rowIndex, , drop = FALSE])
        } else{
          dfEval <- deriv1(doseVec[rowIndex, ], pm[rowIndex, , drop = FALSE])
        }

        varVal <- as.vector(dfEval) %*% varCov %*% as.vector(dfEval)
        retMat[rowIndex, 2] <- sqrt(varVal)

        if (!se.fit)
        {
          retMat[rowIndex, 3:4] <- retMat[rowIndex, 1] + (tquan * sqrt(varVal + sumObjRV)) * c(-1, 1)
        }
        # if (identical(interval, "confidence"))
        if (interval == TRUE)
        {
          retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal)
          retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal)
        }
        }
    }

  } else if(object$fit$method == "NPMLE"){
    # if(!se.fit & identical(interval, "none") | missing(newdata))
    if(se.fit == FALSE & interval == FALSE) # | missing(newdata))
    {
      return(retMat[, 1])
      end()
    } else if(se.fit == T){

      if(is.null(cluster)) df <- object$data else df <- data.frame(object$data, group = cluster)
      splitData <- by(df, object$data[,5], function(x) x)

      # if(!missing(newdata)){
      #   if(!any(names(newdata) %in% cName)){
      #     misLev <- c()
      #     cont <- 1
      #     for(i in 1:length(names(splitData))){
      #       if(any(unique(groupLevels) == names(splitData)[i]) == F)  {
      #         misLev[cont] <- i
      #         cont <- cont + 1}
      #     }
      #
      #     # if(!is.null(misLev)) splitData <- splitData[-misLev]
      #     splitData <- list(splitData[[1]])
      #   }
      # }

      confCal <- function(x, times.pred, B, groups = NULL) {
        L <- rep(x[,1], x[,3])
        R <- rep(x[,2], x[,3])
        if(!is.null(groups)) gr <- rep(x[,length(x[1,])], x[,3]) else gr <- NULL
        cis <- confint.predict(L, R, times.pred, B = B, groups = gr)
        cis
      }
      # splitDoseVec <- by(doseVec, dataSet[,5], function(x) x)
      splitDoseVec <- by(doseVec, groupLevels, function(x) x)

      cis <- list()
      for(i in 1:length(splitData)){
        # i <- 1
        tmp <- confCal(splitData[[i]], splitDoseVec[[i]], B = B, groups = cluster)
        cis[[i]] <- as.data.frame(tmp)
        # cis <- append(cis, as.data.frame(tmp))
      }

      # cis
      # cis <- lapply(splitData, function(x) confCal(x, unique(doseVec), B = B,
      #                                               groups = cluster))
      # cis <- lapply(splitData, function(x) confCal(x, doseVec, B = B,
      #                                              groups = cluster))
      # retMat[, 1] <- as.numeric(mapply(listaFun, groupLevels, doseVec))
      tmp <- data.frame(sapply(cis, function(x) x$se))
      tmp2 <- data.frame(sapply(cis, function(x) x$lower))
      tmp3 <- data.frame(sapply(cis, function(x) x$upper))

      if(length(tmp[1,]) == 1) {
        retMat[, 2] <- as.numeric(unlist(tmp))
        indCol <- 3
        if (interval == TRUE){
          retMat[, indCol] <- as.numeric(unlist(tmp2))
          retMat[, indCol + 1] <- as.numeric(unlist(tmp3))
          colnames(retMat)[indCol] <- "Lower"
          colnames(retMat)[indCol + 1] <- "Upper"
          indCol <- indCol + 2
        }
        retMat <- data.frame(retMat)
        retMat[, indCol] <- object$dataList$names$rNames
      } else {
        colnames(tmp) <- object$dataList$names$rNames
        tmp <- stack(tmp)
        tmp2 <- stack(tmp2)
        tmp3 <- stack(tmp3)

        retMat[, 2] <- tmp[,1]
        retMat <- data.frame(retMat)
        indCol <- 3

        if (interval == TRUE){
          retMat[, indCol] <- tmp2[, 1]
          retMat[, indCol + 1] <- tmp3[, 1]
          colnames(retMat)[indCol] <- "Lower"
          colnames(retMat)[indCol + 1] <- "Upper"
          indCol <- indCol + 2
        }
        retMat[, indCol] <- tmp[,2]
      }

      retMat <- data.frame(retMat)
      retMat[, indCol + 1] <- doseVec
      colnames(retMat)[indCol] <- cName
      colnames(retMat)[indCol + 1] <- "Time"

      return(retMat)
      end()
    }

  } else if(object$fit$method == "KDE"){
    # To be done
  }


  ## Keeping relevant indices
  keepInd <- 1
  if (se.fit) {keepInd <- c(keepInd, 2)}
  # if (!identical(interval, "none")) {keepInd <- c(keepInd, 3, 4)}
  if (interval == T) {keepInd <- c(keepInd, 3, 4)}

  return(retMat[, keepInd])  # , drop = FALSE])
}


