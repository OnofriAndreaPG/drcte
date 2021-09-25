# Functions for bootstrap resampling (23/7/2021)

resample.cens <- function(dat, cluster = NULL, replace) {
  # This function resamples the data with replacement
  # dat is a data.frame, containing the survival data, in the form left/right
  # One row per individual
  # exit early for trivial data
  # if(length(dat[,1]) == 1 || all(replace == FALSE)) {
    # return(dat)
  # }
  if(is.null(cluster)){
    sel <- sample(1:length(dat[,1]), replace = replace[[1]])
    ret <- dat[sel, ]
  } else {
    cls <- sample(unique(cluster), replace = replace[[2]])
    # subset on the sampled clustering factors
    sub <- lapply(cls, function(b) subset(dat, cluster == b))
    sub
    # sample lower levels of hierarchy (if any)
    # if(length(cluster) > 1)
    sub <- lapply(sub, resample.cens, replace = replace[[1]])
    # join and return samples
    ret <- do.call(rbind, sub)
  }
  if(is.vector(ret)){
    ret <- data.frame(ret)
    names(ret) <- names(dat)
    row.names(ret) <- 1:length(dat[,1])
  } else {
    # print(ret)
    row.names(ret) <- 1:length(dat[,1])
  }
  ret
}

simulateTE <- function(start, end, count, B = 1, groups = NULL){
  # This function produces a new time-to-event dataset, by bootstrap
  # resampling the given dataset. It returns a list of datasets
  # set.seed(seed)
  # Superseded !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  res <- vector("list", length = B)
  for(i in 1:B){
    # bootSam <- sample(1:sum(count), replace = T)
    # startSam <- rep(start, count)[bootSam]
    # endSam <- rep(end, count)[bootSam]
    # bootSam <- sample(1:sum(count), replace = T)
    L <- rep(start, count)
    R <- rep(end, count)
    # print(R); print(L); print(count); print(groups)
    # stop()
    if(!is.null(groups)) groups <- rep(groups, count)
    df <- data.frame(L, R)

    if(!is.null(groups)){
      cat("\r Cluster Resampling:", i)
      newSam <- resample.cens(df, cluster = groups, replace = list(TRUE, TRUE))
      # print(newSam)
      stop()
    } else {
      cat("\r Resampling:", i)
      newSam <- resample.cens(df, replace = T)
    }

    startSam <- newSam$L; endSam <- newSam$R
    obj <- drcte:::getNPMLE(survival::Surv(startSam, endSam, type = "interval2") ~ 1)
    newStart <- obj$intmap[1,]
    newEnd <- obj$intmap[2,]
    n <- sum(count)
    newCount <- obj$pf*n
    df <- data.frame(startTime = newStart, endTime = newEnd,
                     count = newCount)
    res[[i]] <- df
  }
  return(res)
}

