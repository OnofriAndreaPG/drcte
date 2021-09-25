\name{plot.drcte}
\alias{plot.drcte}

\title{Plotting fitted time-to-event models}

\description{
  \code{plot} displays fitted curves and observations in the same plot window,
  distinguishing between curves by different plot symbols and line types.
}

\usage{

  \method{plot}{drcte}(x, ..., add = FALSE, level = NULL, shading = TRUE,
  type = c("average", "all", "bars", "none", "obs", "confidence"),
  broken = FALSE, bp, bcontrol = NULL, conName = NULL, axes = TRUE,
  gridsize = 100, log = "", xtsty, xttrim = TRUE,
  xt = NULL, xtlab = NULL, xlab, xlim,
  yt = NULL, ytlab = NULL, ylab, ylim,
  cex, cex.axis = 1, col = FALSE, lty, pch,
  legend, legendText, legendPos, cex.legend = 1,
  normal = FALSE, normRef = 1, confidence.level = 0.95)
}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional graphical arguments. For instance, use \code{lwd=2} or \code{lwd=3} to increase the width of plot symbols.}
  \item{add}{logical. If TRUE then add to already existing plot.}
  \item{level}{vector of character strings. To plot only the curves specified by their identification number.}
  \item{shading}{For NPMLE, it prints shading on the graphs for the points where the likelihood value is not unique. Defaults to TRUE}
  \item{type}{a character string specifying how to plot the data. There are currently
  5 options:  "average" (averages and fitted curve(s); default),
  "none" (only the fitted curve(s)), "obs" (only the data points),
  "all" (all data points and fitted curve(s)),
  "bars" (averages and fitted curve(s) with model-based standard errors (see Details)), and "confidence" (confidence bands for fitted curve(s)). For nonparametric curves, }
  \item{broken}{logical. If TRUE the x axis is broken provided this axis is logarithmic
  (using functionality in the CRAN package 'plotrix').}
  \item{bp}{numeric value specifying the break point below which the dose is zero (the amount of stretching on
  the dose axis above zero in order to create the visual illusion of a logarithmic scale \emph{including} 0).
  The default is the base-10 value corresponding to the rounded value of the minimum of the log10 values of
  all positive dose values. This argument is only working for logarithmic dose axes.}
  \item{bcontrol}{a list with components \code{factor}, \code{style} and \code{width}.
  Controlling the appearance of the break (in case \code{broken} is \code{TRUE}).
  The component \code{factor} is the distance from the control to the break as a
  multiple of the value of \code{bp} (default is 2).
  The component \code{style} can take the values: \code{gap}, \code{slash} or \code{zigzag}.
  The component \code{width} is the width of the break symbol (default is 0.02).}
  \item{conName}{character string. Name on x axis for dose zero. Default is '"0"'.}
  \item{axes}{logical indicating whether both axes should be drawn on the plot.}
  \item{gridsize}{numeric. Number of points in the grid used for plotting the fitted curves.}
  \item{log}{a character string which contains '"x"' if the x axis is to be logarithmic, '"y"' if the y axis is to be logarithmic and '"xy"' or
  '"yx"' if both axes are to be logarithmic. The default is "x". The empty string "" yields the original axes.}
  \item{xtsty}{a character string specifying the dose axis style for arrangement of tick marks. By default ("base10")
  For a logarithmic axis by default only base 10 tick marks are shown ("base10"). Otherwise sensible
  equidistantly located tick marks are shown ("standard"), relying on \code{\link{axTicks}}.}
  \item{xttrim}{logical specifying if the number of tick marks should be trimmed in case too many tick marks
  are initially determined.}
  \item{xt}{a numeric vector containing the positions of the tick marks on the x axis.}
  \item{xtlab}{a vector containing the tick marks on the x axis.}
  \item{xlab}{an optional label for the x axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{yt}{a numeric vector, containing the positions of the tick marks on the y axis.}
  \item{ytlab}{a vector containing the tick marks on the y axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}
  \item{cex}{numeric or numeric vector specifying the size of plotting symbols and text
  (see \code{\link{par}} for details).}
  \item{cex.axis}{numeric value specifying the magnification to be used for axis annotation
  relative to the current setting of cex.}
  \item{col}{either logical or a vector of colours. If TRUE default colours are used.
  If FALSE (default) no colours are used.}
  \item{legend}{logical. If TRUE a legend is displayed.}
  \item{legendText}{a character string or vector of character strings specifying the legend text
  (the position of the upper right corner of the legend box).}
  \item{legendPos}{numeric vector of length 2 giving the position of the legend.}
  \item{cex.legend}{numeric specifying the legend text size.}
  \item{lty}{a numeric vector specifying the line types.}
  \item{pch}{a vector of plotting characters or symbols (see \code{\link{points}}).}
  \item{normal}{logical. If TRUE the plot of the normalized data and fitted curves are shown
  (for details see Weimer et al. (2012) for details).}
  \item{normRef}{numeric specifying the reference for the normalization (default is 1).}
  \item{confidence.level}{confidence level for error bars. Defaults to 0.95.}
}

\details{
  The use of \code{xlim} allows changing the range of the x axis, extrapolating the fitted dose-response curves.
  Note that changing the range on the x axis may also entail a change of the range on the y axis. Sometimes
  it may be useful to extend the upper limit on the y axis (using \code{ylim}) in order to fit a legend into
  the plot.

  See \code{\link{colors}} for the available colours.

  Suitable labels are automatically provided.

  The arguments \code{broken} and \code{bcontrol} rely on the function
  \code{link{axis.break}} with arguments
  \code{style} and \code{brw} in the package \code{plotrix}.

  The model-based standard errors used for the error bars are calculated
  as the fitted value plus/minus the estimated error times the
  1-(alpha/2) quantile in the t distribution with degrees of freedom
  equal to the residual degrees of freedom for the model (or using a
  standard normal distribution in case of binomial and poisson data),
  where alpha=1-confidence.level. The standard errors are obtained using
  the predict method with the arguments interval = "confidence"
  and level=confidence.level.
}

\value{
  An invisible data frame with the values used for plotting the fitted curves. The first column contains the dose values,
and the following columns (one for each curve) contain the fitted response values.
}

\author{
Christian Ritz and Jens C. Streibig. Contributions from Xiaoyan Wang and Greg Warnes.
}

\references{
  Weimer, M., Jiang, X., Ponta, O., Stanzel, S., Freyberger, A., Kopp-Schneider, A. (2012)
  The impact of data transformations on concentration-response modeling.
  \emph{Toxicology Letters}, \bold{213}, 292--298.
}

%\note{ ~~further notes~~ }

%\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{
library(drcte)
data(verbascum)
mod <- drmte(nSeeds ~ timeBef + timeAf, fct = NPMLE(),
                 data = verbascum, curveid = Species)
plot(mod)
}
\keyword{plot}

