% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precision_functions.R
\name{bias.plot}
\alias{bias.plot}
\title{Bias Plot}
\usage{
bias.plot(
  Data,
  R1.label = "Reader 1",
  R2.label = "Reader 2",
  a = 0,
  b = 1,
  pch = 16,
  col.points = "black",
  lty.ab = 2,
  lty.col = "black",
  err.col = "black"
)
}
\arguments{
\item{Data}{Numeric data frame having columns titled the following
\describe{
    \item{Age}{numeric vector representing the unique ages assigned to
        indvidual fish by Reader 1}
    \item{Avg}{numeric vector representing the average age given by
        Reader 2 for fish aged "x' by Reader 1}
    \item{CI}{numeric vector representing the confidence interval of ages
        given by Reader 2 for fish aged "x" by Reader 1}
    }
Other columns may be present in the data frame.}

\item{R1.label}{Text label for the x-axis (defaults to Reader 1) for the
resultant plot}

\item{R2.label}{Text label for the y-axis (defaults to Reader 2) for the
resultant plot}

\item{a}{numeric; y-axis intercept for comparison line}

\item{b}{numeric; slope estimate for comparison line}

\item{pch}{see \code{\link{pch}}; symbol for points used in plot}

\item{col.points}{A specification for the plotting color of the points.  See
description of \emph{Color Specificaiton} in the \code{\link{par}} for
more information.}

\item{lty.ab}{specify line type for the "ab" line.  See description of
\emph{lty} in the \code{\link{par}} for more information}

\item{lty.col}{specify the line color for the "ab" line.  See
description of \emph{Color Specificaiton} in the \code{\link{par}} for
more information.}

\item{err.col}{specify the line color for the error bars.  See
description of \emph{Color Specificaiton} in the \code{\link{par}} for
more information.}
}
\value{
A \code{\link{xyplot}} using the package lattice depicting the
   bias plot between Reader 2 and Reader 1 ages.
}
\description{
\code{bias.plot} creates a bias plot (see Campana 2001) that helps the
   researcher assess whether bias is present between two sets of age reads
   on the same individual fish.  The plot is constructed using lattice
   graphics.
}
\examples{
# Example with RedPorgy data
data(RedPorgy)
df <- with(RedPorgy, aggregate(R4.Age ~ R1.Age, FUN = describe, digits = 4 ))
df <- cbind(df[-ncol(df)], df[[ncol(df)]])
colnames(df) <- c("Age" ,"n", "Avg", "StDev", "StErr", "CI", "CV","RSE",
   "RMSD", "MAE")
bias.plot(Data = df)
}
\seealso{
Other Bias Plot Functions: 
\code{\link{dist.plot}()}
}
\concept{Bias Plot Functions}
