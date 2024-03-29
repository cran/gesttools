% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gestboot.R
\name{gestboot}
\alias{gestboot}
\title{Percentile Based Bootstrap Confidence Intervals}
\usage{
gestboot(
  gestfunc,
  data,
  idvar,
  timevar,
  Yn,
  An,
  Cn,
  outcomemodels,
  propensitymodel,
  censoringmodel = NULL,
  type,
  EfmVar = NA,
  cutoff,
  bn,
  alpha = 0.05,
  onesided = "twosided",
  seed = NULL,
  ...
)
}
\arguments{
\item{gestfunc}{Name (without quotations) of the g-estimation function to run.
One of \code{gestSingle} or \code{gestMultiple}.}

\item{data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel, censoringmodel, type, EfmVar, cutoff}{Same arguments as in gest functions, to be input into gestfunc.}

\item{bn}{Number of bootstrapped datasets.}

\item{alpha}{Confidence level of confidence intervals.}

\item{onesided}{Controls the type of confidence interval generated. Takes one of three inputs, \code{"upper"} for upper one-sided confidence intervals,
\code{"lower"} for lower one-sided confidence intervals, and \code{"twosided"} for two-sided confidence intervals. Defaults to \code{"twosided"}.}

\item{seed}{Integer specifying the random seed for generation of bootstrap samples.}

\item{...}{additional arguments.}
}
\value{
Returns a list of the following four elements.
\item{original }{The value of the causal parameters estimated on the original data \code{data}.}
\item{mean.boot }{The average values of the causal parameters estimated on the bootstrapped datasets.}
\item{conf }{The upper and/or lower bounds of \eqn{1-\alpha} confidence intervals for each element of \eqn{\psi}.
For example, if \code{type=2}, and \eqn{\psi=(\psi_0,\psi_1)}, a separate confidence interval is fitted for \eqn{\psi_0} and \eqn{\psi_1}.}
\item{conf.Bonferroni }{The upper and/or lower bounds of Bonferroni corrected confidence
intervals for \eqn{\psi}, used for multiple comparisons.}
\item{boot.results}{A tibble containing the result for each bootstrapped dataset}
}
\description{
Generates percentile based confidence intervals for the causal parameters
of a fitted SNMM. Bonferroni corrected confidence
intervals are also reported for multiple comparisons.
}
\examples{
datas <- dataexamples(n = 1000, seed = 123, Censoring = FALSE)
data <- datas$datagest
data <- FormatData(
  data = data, idvar = "id", timevar = "time", An = "A",
  varying = c("A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
)
idvar <- "id"
timevar <- "time"
Yn <- "Y"
An <- "A"
Cn <- NA
outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A")
propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
censoringmodel <- NULL
type <- 1
EfmVar <- NA
bn <- 5
alpha <- 0.05
gestfunc <- gestSingle
gestboot(gestfunc, data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
  censoringmodel = NULL, type = 1, EfmVar,
  bn = bn, alpha = alpha, onesided = "twosided", seed = 123
)
}
