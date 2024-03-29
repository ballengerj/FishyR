% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_functions.R
\name{ZI.Index.Boot}
\alias{ZI.Index.Boot}
\title{Zero-Inflated Model Bootstrap Function}
\usage{
ZI.Index.Boot(data, i, bestmod, dist = "negbin", link = "logit", grid)
}
\arguments{
\item{data}{see \code{\link[boot]{boot}}}

\item{i}{see \code{\link[boot]{boot}} for \code{stype}}

\item{bestmod}{character vector representing saved object representing the best
(i.e., final) model}

\item{dist}{character vector used to specify the distribution to use in the
\code{\link[pscl]{zeroinfl}} function being implemented}

\item{link}{character vector used to specify the link to use in the
\code{\link[pscl]{zeroinfl}} function being implemented}

\item{grid}{prediction data frame for which you will get your predicted
response variable at.  Generally will be created via a call to
\code{\link[FishyR]{Index.Grid}}}
}
\value{
numeric vector with four points:
  \describe{
    \item{Points}{}
    \item{Group}{}
    \item{Factor}{}
    \item{Value}{}
  }
}
\description{
\code{ZI.Index.Boot} calculates a raw, a normalized, a centered, and a scaled
  relative abundance index for each bootstrap iteration
}
\seealso{
\code{\link[boot]{boot}}

Other Index Bootstrap Functions: 
\code{\link{Boot.Diag}()},
\code{\link{Stat.Calc}()}
}
\concept{Index Bootstrap Functions}
