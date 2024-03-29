% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importance_rank.R
\name{importance_rank}
\alias{importance_rank}
\title{Importance Ranking for a Stressed Model}
\usage{
importance_rank(
  object,
  xCol = "all",
  wCol = "all",
  type = c("Gamma", "Wasserstein", "reverse", "all"),
  f = NULL,
  k = NULL,
  s = NULL
)
}
\arguments{
\item{object}{A \code{SWIM} or \code{SWIMw} object.}

\item{xCol}{Numeric or character vector, (names of) the columns
of the underlying data of the \code{object}
(\code{default = "all"}). If \code{xCol = NULL}, only
the transformed data \code{f(x)} is considered.}

\item{wCol}{Vector, the columns of the scenario weights 
of the \code{object} corresponding to different 
stresses (\code{default = "all"}).}

\item{type}{Character, one of \code{"Gamma", "Wasserstein", "Kolmogorov",
"reverse", "all"}.}

\item{f}{A function, or list of functions, that, applied to
\code{x}, constitute the transformation of the data
for which the sensitivity is calculated.}

\item{k}{A vector or list of vectors, same length as \code{f},
indicating which columns of \code{x} each function
in \code{f} operates on.\cr
When \code{f} is a list, \code{k[[i]]} corresponds
to the input variables of \code{f[[i]]}.}

\item{s}{A function that, applied to \code{x}, defines the reverse
sensitivity measure. If \code{type = "reverse"} and 
\code{s = NULL}, defaults to \code{type = "Gamma"}.}
}
\value{
A data.frame containing the importance ranks of the
    stressed model for different sensitivity measures. Small values
    correspond to large sensitivities. Different rows correspond
    to different random variables. The first two rows specify the
    \code{stress} and \code{type} of the sensitivity measure on
    which the ranking is calculated.
}
\description{
Provides the importance ranks of the components (random variables)
    of a stressed model for different sensitivity measures.
}
\details{
For the definition of the sensitivity
    measures (\code{type}), see \code{\link{sensitivity}}.
}
\examples{
## example with a stress on VaR
set.seed(0)
x <- as.data.frame(cbind(
  "log-normal" = rlnorm(1000),
  "gamma" = rgamma(1000, shape = 2)))
res1 <- stress(type = "VaR", x = x,
  alpha = c(0.9, 0.95), q_ratio = 1.05)

importance_rank(res1, wCol = 1:2, type = "Gamma")
## sensitivity of log-transformed data
importance_rank(res1, wCol = 1, type = "all",
  f = list(function(x)log(x), function(x)log(x)), k = list(1, 2))

}
\seealso{
See \code{\link{sensitivity}} for the values of the
    sensitivity measures, \code{\link{plot_sensitivity}} for plotting
    sensitivity measures and \code{\link{summary}} for a
    summary statistic of a stressed model.
}
\author{
Silvana M. Pesenti
}
