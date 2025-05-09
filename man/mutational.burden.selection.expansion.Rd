% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clone_size_distribution.R
\name{mutational.burden.selection.expansion}
\alias{mutational.burden.selection.expansion}
\title{Mutation accumulation during exponential expansion with clonal selection.}
\usage{
mutational.burden.selection.expansion(mu, lambda, delta, s, t.s, t.end, b)
}
\arguments{
\item{mu}{mutation rate per cell division}

\item{lambda}{proliferation rate}

\item{delta}{loss rate}

\item{s}{selective advantage}

\item{t.s}{time point at which selective advantage is acquired.}

\item{t.end}{end point}

\item{b}{minimal clone size of interest. Number or vector.}

\item{N}{population size}
}
\value{
This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least \code{b} cells.
}
\description{
Mutation accumulation during exponential expansion with clonal selection.
}
