\name{quantile.drcte}
\alias{quantile.drcte}

\title{Estimating quantiles for time-to-event models}

\description{
This function estimates the quantiles for time-to-events data, either for the whole sample (restricted = F) or only for the fraction of individuals with events (restricted = T). The quantiles for rates (the inverse of times) are also available, by using the 'rate = T' option}

\usage{
  \method{quantile}{drcte}(x, probs, restricted = FALSE,
            interval = c("none", "delta", "boot"),
            clevel = NULL, level = ifelse(!(interval == "none"), 0.95, NULL),
            bound = TRUE, od = FALSE, vcov. = vcov, # robust = false,
            display = TRUE, rate = F, B = 1000, seed = 1234,  ...)
}

\arguments{
  \item{x}{an object of class 'drcte'.}
  \item{probs}{a vector with the sought probabilities; the values must be included between 0 and 1. }
  \item{restricted}{a logical value; if vector with the sought probabilities; the values must be included between 0 and 1. If FALSE, the whole sample is considered (the default), if TRUE, only the fraction of individuals with events are considered}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". See Details below for more explanation.}
  \item{clevel}{character string specifying the curve id in case on estimates for a specific curve is requested. By default estimates are shown for all curves.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default, but \code{sandwich} is also an option (for obtaining robust standard errors).}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{rate}{Logical: if FALSE quantiles for time-to-events are reported (default), if TRUE, rates (inverse of times) are reported.}
  \item{B}{Number of bootstrap resamples for bootstrap based inference}
  \item{seed}{seed for sampling methods}
  \item{...}{see the details section below.}
}

\details{
There are several options for calculating confidence intervals through the argument \code{interval}. The option "none" (default) does not provide confidence intervals. The option "delta" results in asymptotical Wald-type confidence intervals (using the delta method and the normal distribution) and it is only available for parametric time-to-event models. The option "boot" provides bootraspped confidence intervals and it is only available for NPMLE or KDE curves; it may be rather slow. For (\code{\link{KDE}} models, the additional arguments \code{lower} and \code{upper} may be supplied. These arguments specify the lower and upper limits of the bisection method used to find the ED values.  The lower and upper limits need to be smaller/larger than the EDx level to be calculated. The default limits are 0 and 1}
\value{
  An invisible matrix containing the shown matrix with two or more columns, containing the estimates and the corresponding estimated standard errors and possibly lower and upper confidence limits. Or, alternatively, a list with elements that may be plugged directly into \code{parm} in the package \emph{multcomp} (in case the argument \code{multcomp} is TRUE).
}

\seealso{
The 'ED()' function.
}

\author{Andrea Onofri}

\examples{
rm(list = ls())
start <- c(0, 6, 7, 10, 13, 16, 22, 23)
end <- c(6, 7, 10, 13, 16, 22, 23, Inf)
count <- c(3, 9, 3, 8, 1, 1, 5, 5)

# Fitting a non-parametric time-to-event model
mod <- drmte(count ~ start + end, fct = NPMLE())

# Getting the times required to reach three different fractions
# of individuals with events (i.e.: 0.1, 0.3 and 0.5), either including
# the whole lot of individuals, or only the individuals with events
quantile(mod, probs = c(0.10, 0.30, 0.50))
quantile(mod, probs = c(0.10, 0.30, 0.50), restricted = TRUE)
}

\keyword{models}
\keyword{time-to-event}
