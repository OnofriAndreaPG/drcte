\name{lognormal}
\alias{lognormal}
\alias{lognormal.2}
\title{
Log-normal distribution of times to an event
}
\description{
This function provides the truncated log-normal cumulative distribution of event times, to be used to fit time-to-event models.
}
\usage{
lognormal(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
lognormal.2(fixed = c(NA, NA, NA), names = c("log-b", "logit-d", "log-e"))
}

\arguments{
\item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed. It should only be used to fix the 'd' parameter to 1, which provides the usual log-logistic distribution, with no fraction of individuals that are never expected to get the event .
}
\item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.
}
}
\details{
The log-normal distribution of germination times is parameterised as:

\deqn{ P(t) = d \cdot \textrm{pnorm}((log(t + 0.000001) - e)/b)}

where 't' is the time and 'P(t)' is the proportion of individual with the event at time = t. The parameter 'e' is the median event time with a log-link, 'd' is the fraction of individuals that are susceptible to get the event (on a logit link) and 'b' is the standard deviation on a log-scale and a log-link. When d = 1, this parameterisation gives same results as 'survereg' in the 'survival package', with 'dist = "lognormal". The function 'lognormal.2()' has a similar parameterisation:

\deqn{ P(t) = d \cdot pnorm(b \cdot (log(t + 0.000001) - log(e)))}

which is close to the 'lnormal()' or 'LN.3()' function in 'drc' and parameter values are not given any link function.
}
\value{
The value returned is a list containing the nonlinear function, the self starter function and the parameter names.
}
\references{
Kleinbaum, D.G., Klein, M., 2005. Survival analysis. Springer Science, USA.
}
\author{
Andrea Onofri
}
\examples{
library(drcte)
data(chickweed)
modb <- drmte(count ~ start + end,
               data = chickweed,
               fct = lognormal())
}
