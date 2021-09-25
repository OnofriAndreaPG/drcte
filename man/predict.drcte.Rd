\name{predict.drct}
\alias{predict.drcte}

\title{Prediction}

\description{
  Predicting probabilities for models of class 'drcte'.
}

\usage{

  \method{predict}{drcte}(object, newdata, se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95, na.action = na.pass, od = FALSE, vcov. = vcov,
  NPMLE.method = "interpolation", ...)
}

\arguments{
  \item{object}{an object of class 'drcte'.}
  \item{newdata}{An optional data frame in which to look for variables with which to predict.
  If omitted, the fitted values are used.}
  \item{se.fit}{logical. If TRUE standard errors are required.}
  \item{interval}{character string. Type of interval calculation: "none", "confidence".}
  \item{level}{Tolerance/confidence level.}
  \item{na.action}{function determining what should be done with missing values in 'newdata'.
  The default is to predict 'NA'.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{vcov.}{function providing the variance-covariance matrix. \code{\link{vcov}} is the default,
  but \code{sandwich} is also an option (for obtaining robust standard errors).}
  \item{NPMLE.method}{character string relating to the type of prediction for NPMLE. For interval censored observations there is not a unique MLE for a specific time, and, therefore, predictions can be obtained by three methods: "interpolation" (default; predictions are sought on the line connecting the two points bounding the non-unique MLE interval), "left" (taking the left side of the non-unique MLE interval) or "right" (taking the right side of the non-unique MLE interval. The argument is neglected for parametric or kernel based time-to-event models}
  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{
  A matrix with as many rows as there are dose values provided in 'newdata' or in the original dataset
  (in case 'newdata' is not specified) and, at most, 4 columns containing fitted, standard errors, lower and
  upper limits of confidence/prediction intervals.
}

%\references{}

\author{Andrea Onofri, borrowing code from Christian Ritz}

%\note{}

\seealso{For details are found in the help page for \code{\link{predict.lm}}.}

\examples{

data(chickweed)
modb <- drmte(count ~ start + end,
               data=chickweed,
               fct = LL.3())
predict(modb, se = TRUE, interval = "confidence")
}

\keyword{seed germination}
\keyword{time-to-event model}
\concept{prediction}
