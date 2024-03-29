% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_hist.R
\name{plot_hist}
\alias{plot_hist}
\title{Plotting Histograms of a Stressed Model}
\usage{
plot_hist(
  object,
  xCol = 1,
  wCol = "all",
  base = FALSE,
  x_limits,
  displ = TRUE,
  binwidth,
  displLines = FALSE
)
}
\arguments{
\item{object}{A \code{SWIM} or \code{SWIMw} object.}

\item{xCol}{Numeric or character, (name of) the column of the underlying data
of the \code{object} (\code{default = 1}).}

\item{wCol}{Vector, the columns of the scenario weights 
of the \code{object} corresponding to different 
stresses (\code{default = "all"}).}

\item{base}{Logical, if \code{TRUE}, statistics under the baseline 
are also returned (\code{default = "FALSE"}).}

\item{x_limits}{Vector, the limits of the x-axis of the plot, the
value for \code{xlim} in the \code{coord_cartesian}
function in \code{ggplot}.}

\item{displ}{Logical, if \code{TRUE} the plot is displayed,
otherwise the data.frame for customised plotting with
\code{ggplot} is returned (\code{default = TRUE}).}

\item{binwidth}{Numeric, the width of the bins used to plot
the histogram, the \code{binwidth} in the
\code{geom_freqpoly} function  in \code{ggplot}
(default corresponds to 30 bins).}

\item{displLines}{Logical, if \code{TRUE} lines are displayed
instead of bins (\code{default = FALSE}).}
}
\value{
If \code{displ = TRUE}, a histogram of the stochastic model
    under the scenario weights.

    If \code{displ = FALSE}, a data.frame for customised plotting with
    \code{ggplot}. The data.frame contains the columns: the column,
    \code{xCol}, of the data of the stressed model, \code{stress}
    (the stresses) and \code{value} (the values). \cr
    Denote by \code{res} the return of the function call, then
    \code{ggplot} can be called via:
    \deqn{ggplot(res, aes(x = res[ ,1], y = ..density.., weight = value)))}
    \deqn{ + geom_{freqpoly}(binwidth, aes(color = factor(stress))).}
}
\description{
Plots the histogram of a stressed model component (random variable)
    under the scenario weights.
}
\examples{
## example with a stress on VaR
set.seed(0)
x <- data.frame("gamma" = rgamma(10^5, shape = 2))
res1 <- stress(type = "VaR", x = x,
  alpha = c(0.75, 0.95), q_ratio = 1.1)
plot_hist(res1, xCol = "gamma", wCol = 1:2, base = TRUE, binwidth = 0.4)
plot_hist(res1, xCol = "gamma", wCol = 1:2, base = TRUE, binwidth = 0.4, displLines = TRUE)

}
\seealso{
See \code{\link{cdf}} and \code{\link{plot_cdf}} for
    values and plotting of the empirical distribution
    function of a stressed model, respectively, and
    \code{\link{quantile_stressed}} for sample quantiles of
    a stressed model.
}
