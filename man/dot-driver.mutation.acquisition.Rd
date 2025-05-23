% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulate_trees.R
\name{.driver.mutation.acquisition}
\alias{.driver.mutation.acquisition}
\title{Driver mutation acquisition}
\usage{
.driver.mutation.acquisition(mut.rate.D, s.shape, s.rate)
}
\arguments{
\item{mut.rate.D}{integer; driver mutation rate per cell division}

\item{s.shape}{shape parameter of the gamma distribution from which the selective advantage of a new driver is drawn}

\item{s.rate}{rate parameter of the gamma distribution from which the selective advantage of a new driver is drawn}
}
\value{
the additional selective advantage acquired during the cell division
}
\description{
Simulate the acquisition of new driver mutations in both daughter cells
}
