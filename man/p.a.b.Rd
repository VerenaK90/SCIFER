% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clone_size_distribution.R
\name{p.a.b}
\alias{p.a.b}
\title{Clone size distribution in a noncritical birth-death process (approximate).}
\usage{
p.a.b(lambda, delta, t, a, b, mode = "cumulative", approx = "highnumbers")
}
\arguments{
\item{lambda}{proliferation rate}

\item{delta}{loss rate}

\item{t}{time}

\item{a}{clone size at t=0}

\item{b}{clone size at t=t}

\item{mode}{"density" if density distribution is to be returned , "cumulative" if cumulative distribution is to be returned. Defaults to "cumulative"}

\item{approx}{Approximation to be used. Defaults to "highnumbers"; i.e. the distribution is approximated with a gamma distribution if \code{a} and \code{b} are large.}
}
\value{
The probability of growing from size a to size b within t. The Function switches between the exact solution and an approximate solution according to a parametrized gamma distribution.
}
\description{
Probability of a clone of size "a" to grow to size "b" within "t" according to a noncritical linear birth-death process.
}
\details{
If \code{approx="highnumbers"}, the function is approximated with a \eqn{\Gamma}-distribution if \code{a+b>100} and \code{mode="density"} or if \code{a+b>10} and \code{mode="cumulative"}. The \eqn{\Gamma}-distribution is parametrized with
\eqn{shape = \mu^2/\sigma, scale = \sigma/\mu}, where
\eqn{\mu = a e^{(\lambda - \delta)t}, \sigma = a \frac{\lambda + \delta}{\lambda - \delta}e^{(\lambda - \delta)t}(e^{(\lambda - \delta)t}-1)}
}
\references{
Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
}
