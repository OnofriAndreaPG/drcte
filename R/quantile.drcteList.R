##########################################################################
quantile.drcteList <- function(x, probs, restricted = FALSE,
                           interval = c("none", "delta", "boot"),
                           clevel = NULL, level = ifelse(!(interval == "none"), 0.95, NULL),
                           bound = TRUE, od = FALSE, vcov. = vcov, # robust = false,
                           display = FALSE, rate = F, B = 999, ...){


  quantileList <- function(el, probs, rate, display){

  EDlist <- try(quantile(el, probs = probs,
                         rate = rate, display = F),
                silent = T)
  if(any(class(EDlist) == "try-error") == T) {
    estim <- ifelse(rate == T, 0, NA)
    estim <- rep(estim, length(probs))
    SE <- rep(NA, length(probs))

    retDF <- data.frame(estim,
                        rLab = SE)
    colnames(retDF) <- c("Estimate", "SE")
    row.names(retDF) <- paste(":", probs*100, "%", sep = "")
    EDlist <- retDF
    }
  EDlist
  }

  ret <- lapply(x$separateFit, quantileList, probs = probs,
              rate = rate, display = F)
  retFin <- do.call(rbind, ret)
  row.names(retFin) <- sub(".1", ".", row.names(retFin), fixed = T)
  row.names(retFin) <- sub(".:", ":", row.names(retFin), fixed = T)
  return(retFin)
}
