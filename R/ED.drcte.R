ED.drcte <- function(object, respLev, interval = c("none", "delta", "boot"),
         clevel = NULL, level = ifelse(!(interval == "none"), 0.95, NULL),
         reference = c("control", "upper"), type = c("relative", "absolute"),
         lref, uref, bound = TRUE, od = FALSE, vcov. = vcov, # robust = false,
         display = TRUE, pool = TRUE, logBase = NULL, multcomp = FALSE,
         intType = "confidence", rate = F, B = 200, seed = 1234,
         units = NULL, ...)
{
    # Get information from the call
    interval <- match.arg(interval)
    reference <- match.arg(reference)
    type <- match.arg(type)

    # print(is.function(vcov.))
    # stop()

    ## Checking 'respLev' vector ... should be numbers between 0 and 100
    if ( (type == "relative") && (bound) )
    {
        if (any(respLev <= 0 | respLev >= 100))
        {
            stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
        }
    }

    if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){

        choice <- object$fct$name

        class(object) <- "drc"
        if(choice  == "LL.4" |  choice == "LL.3" | choice == "LL.2" |
           choice  == "LN.4" |  choice == "LN.3" | choice == "LN.2" |
           choice  == "W1.4" |  choice == "W1.3" | choice == "W1.2" |
           choice  == "W2.4" |  choice == "W2.3" | choice == "W2.2") {
            probLev <- respLev
            } else {
                if(type == "absolute") probLev <- respLev else probLev <- respLev/100
            }

        # if(type == "absolute") respLev2 <- respLev/100 else respLev2 <- respLev/100
        if(!is.null(units)){
            vcov. <- function(x) sandwich::vcovCL(x, cluster = units)
        }
        EDmat <- suppressWarnings(drc::ED(object = object,
                 respLev = probLev, #ifelse(type == "absolute", respLev/100, respLev),
                 interval = interval,
         clevel = clevel, level = level,
         reference = reference, type = type, lref = lref,  uref = uref,
         bound = bound, od = od, vcov. = vcov., # robust = false,
         display = F, pool = pool, logBase = logBase, multcomp = multcomp,
         ...))

        if(rate == T){
            GRval <- 1/EDmat[,1]
            GRes <- (EDmat[,2] * 1/EDmat[,1]^2)
            GRval[is.nan(GRval)==T] <- 0
            GRes[is.nan(GRes) == T] <- 0
            EDmat <- data.frame(Estimate = GRval, SE = GRes)
            if(interval == "delta"){
                EDmat$Lower <- EDmat$Estimate - EDmat$SE * qnorm(1 - (1 - level)/2)
                EDmat$Upper <- EDmat$Estimate + EDmat$SE * qnorm(1 - (1 - level)/2)
            } else if(interval == "inv"){
                # Da fare
            }
        }
        rn <- row.names(EDmat)
        if(type == "absolute") rowProb <- respLev*100 else rowProb <- respLev
        rn3 <- paste(sub(":.*", x = gsub("e:", "", rn, fixed=T), replacement = ""), rowProb, sep = ":")
        rn3 <- paste(rn3, "%", sep = "")
        row.names(EDmat) <- rn3


    } else if(object$fit$method == "NPMLE"){

        ## NPMLE

        ## Retrieving relevant quantities
        if(type == "absolute") probLev <- respLev else probLev <- respLev/100

        # obj <- object$fit$icfitObj
        if(is.null(units)) df <- object$data else df <- data.frame(object$data, group = units)

        obj <- by(df, object$data[,5], function(x) x)

        if(interval != "boot"){
            EDlist <- lapply(obj, function(x) quantileNPMLE(x[,1], x[,2], x[,3],
                                   probs = probLev, type = type, rate = rate))
            EDmat <- data.frame("Estimate" = unlist(EDlist))
        } else {
           if(!is.null(units)) {
               EDlist <- lapply(obj, function(x) quantileNPMLE.boot(x[,1], x[,2], x[,3],
                                   probs = probLev, type = type, rate = rate, B = B,
                                   cluster = x[,7]))
           } else {
               EDlist <- lapply(obj, function(x) quantileNPMLE.boot(x[,1], x[,2], x[,3],
                                   probs = probLev, type = type, rate = rate, B = B))
           }

           EDret <- function(x) {
               nam <- names(x)
               Xmat <- data.frame(x)
               Xmat
            }
           EDmat <- do.call(rbind, lapply(EDlist, EDret))
        }
    } else if(object$fit$method == "KDE") {
        if(type == "absolute") probLev <- respLev else probLev <- respLev/100
        obj <- object$curve[[1]]
        obj2 <- by(object$ICfit$icfit, factor(object$ICfit$icfit[,1]), function(x) x)
        if(interval != "boot"){
            EDlist <- lapply(obj, quantileKDE,
                             probs = probLev, type = type, rate = rate)
            EDmat <- data.frame("Estimate" = unlist(EDlist))
        } else {
           EDlist <- lapply(obj2, function(x) quantileKDE.boot(x$startTime,
                                                                 x$endTime, x$count,
                                   probs = probLev, type = type, rate = rate,
                                   B = B, seed = seed))
           EDret <- function(x) {
                Xmat <- data.frame(x)
                Xmat
                }
           EDmat <- do.call(rbind, lapply(EDlist, EDret))
        }
    }
    if (display)
    {
        cat("\n")
        cat(paste("Estimated times-to-event", "\n", sep = ""))
        if (identical(interval, "boot"))
         {
            if(is.null(units)) {
               intervalText <- paste("(bootstrap-based inference)\n", sep = "")
            } else {
               intervalText <- paste("(cluster robust bootstrap-based inference)\n", sep = "")
            }

            cat(intervalText)
         }
        cat("\n")
        print(EDmat, digits = 5)
    }
    invisible(EDmat)
}
