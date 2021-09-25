\name{summary.drcte}
\alias{summary.drcte}

\title{Summarising time-to-event fits}

\description{
  'summary' compiles a comprehensive summary for objects of class 'drc'.
}

\usage{

  \method{summary}{drcte}(object, od = FALSE, pool = TRUE, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.}
  \item{...}{additional arguments.}
}

\value{
  A list of summary statistics that includes parameter estimates and estimated standard errors.
}

\author{Andrea Onofri, Christian Ritz}



%\examples{}

\keyword{models}
\keyword{time-to-event}
