% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{coef.sclr}
\alias{coef.sclr}
\alias{vcov.sclr}
\alias{confint.sclr}
\alias{model.matrix.sclr}
\alias{model.frame.sclr}
\alias{logLik.sclr}
\title{ML estimate components}
\usage{
\method{coef}{sclr}(object, ...)

\method{vcov}{sclr}(object, ...)

\method{confint}{sclr}(object, parm, level = 0.95, ...)

\method{model.matrix}{sclr}(object, ...)

\method{model.frame}{sclr}(formula, ...)

\method{logLik}{sclr}(object, ...)
}
\arguments{
\item{object, formula}{An object returned by \code{\link{sclr}}.}

\item{...}{Not used. Needed to match generic signature.}

\item{parm}{Parameter name, if missing, all parameters are considered.}

\item{level}{Confidence level.}
}
\description{
\code{coef} returns MLE's.
\code{vcov} returns the estimated variance-covariance matrix at MLE's. 
\code{confint} returns the confidence interval.
\code{model.matrix} returns the model matrix (x).
\code{model.frame} returns the model frame (x and y in one matrix).
}
