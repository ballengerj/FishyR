% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VBLRT_Function.R
\name{VBLRT}
\alias{VBLRT}
\title{Von Bertalanffy Growth Function Likelihood Ratio Test}
\usage{
VBLRT(
  len = NULL,
  age = NULL,
  group = NULL,
  error = 1,
  select = 1,
  Linf = NULL,
  k = NULL,
  t0 = NULL
)
}
\arguments{
\item{len}{Length data}

\item{age}{Age data}

\item{group}{Grouping Variable (only 2 groups allowed)}

\item{error}{Error Type}

\item{select}{TBD}

\item{Linf}{Starting estimate for Linf}

\item{k}{Starting estimate for k}

\item{t0}{Starting estimate for t0}
}
\value{
Individual von Bertalanffy growth function model fits, based on data
  and specifications defined, along with results of likelihood ratio tests
  and information criterion comparing the fits of the different models.
}
\description{
\code{VBLRT} constructs a series of von Bertalanffy growth function models
  looking at likelihood ratio tests to determine if there are differences
  between defined groups.  Here we only allow for a single grouping variable
  with two levels (i.e., sex (Male and Female)).
}
