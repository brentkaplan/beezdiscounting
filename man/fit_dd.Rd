% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_dd}
\alias{fit_dd}
\title{Fit Delay-Discounting Model}
\usage{
fit_dd(dat, equation, method)
}
\arguments{
\item{dat}{A data frame containing delay (\code{x}) and indifference point (\code{y}) data. For "two stage" methods, the data must include an \code{id} column to identify participants.}

\item{equation}{A character string specifying the delay-discounting equation to use. Options include:
\itemize{
\item \code{"mazur"} or \code{"hyperbolic"}: Hyperbolic delay-discounting model (\eqn{y = 1 / (1 + k \cdot x)}).
\item \code{"exponential"}: Exponential delay-discounting model (\eqn{y = \exp(-k \cdot x)}).
}}

\item{method}{A character string specifying the method for fitting the model. Options include:
\itemize{
\item \code{"pooled"} or \code{"agg"}: Fits the model using pooled data.
\item \code{"mean"}: Fits the model using the mean of indifference points at each delay.
\item \code{"ts"} or \code{"two stage"}: Fits the model separately for each participant (requires an \code{id} column in \code{dat}).
}}
}
\value{
A list object of class \code{"fit_dd"}, containing:
\itemize{
\item The fitted model(s).
\item The original dataset (\code{dat}).
\item The specified method (\code{method}).
}
}
\description{
This function fits a delay-discounting model to the given dataset using the specified equation and method.
}
\examples{
data <- data.frame(
  id = rep(1:2, each = 6),
  x = rep(c(1, 7, 30, 90, 180, 365), 2),
  y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
)
fit_dd(data, equation = "mazur", method = "two stage")
}
