% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Francis_bomb_functions.R
\name{sim.hdist}
\alias{sim.hdist}
\title{Simulate h-Statistic Distribution}
\usage{
sim.hdist(
  ref,
  test,
  nsim = 5000,
  diag = 0,
  age = "age",
  age.err = "age.err",
  cy = "catch.year",
  y = "C14",
  y.err = "C14.err",
  s.age = "samp.age",
  b.yr = "C14year"
)
}
\arguments{
\item{ref}{reference data used to develop a reference line; a data frame in
which each row corresponds to one reference datum.}

\item{test}{data frame containing the data from the species whose ages are
being validated.}

\item{nsim}{number of simulated data sets.}

\item{diag}{integer, with possible values of 0, 1 and 2 that controls what,
if any, diagnostic plots are plotted:
\describe{
  \item{0}{no diagnostic plots are constructed}
  \item{1}{plot actual data (refrence and test) and also the data
     considered 'true' for the simulations; does no simulations; returns
     the 'true test data.}
  \item{2}{plots one panel per simulation, 16 panels per page, with each
     panel showing the simulated (solid line) and true (broken line)
     reference lines and the simulated test data, with the estimated
     h-statistic given above the panel.  (Note - use of this option is
     not sensible when nsim is 'large')}
  }}

\item{age}{TBD}

\item{age.err}{TBD}

\item{cy}{TBD}

\item{y}{TBD}

\item{y.err}{TBD}

\item{s.age}{TBD}

\item{b.yr}{TBD}
}
\value{
A dataframe with components (each with length nrow(samp)) with
   components X, Y - x & y coordinates of closest points - and optionally
   (if incl.dist = TRUE) dist - signed distance from reference curve.
   Optionally plots (if plot = TRUE) the data, including lines joining each
   sample point to its closest point.
}
\description{
\code{sim.hdist} Uses simulation to generate the distribution of the
   statistic h, which measures how well we can expect test and reference
   bomb radiocarbon data sets to agree, given the extent of sampling error
   associated with each.  A positive (or negative) value of h implies that
   the test data is, on average, to the right (or left) of the reference
   line.
}
\examples{
data(snapper)
snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
data(bluenose)
bluenose.red <- bluenose[bluenose$C14year >= 1955 &
  bluenose$C14year <= 1972, ]
bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
  bluenose.red$C14 <= 94.0, ]
ref <- get.ref.line(snapper.red$C14year, snapper.red$C14)
sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 16, diag = 0)
sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 16, diag = 1)
sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 16, diag = 2)
}
\seealso{
Other Bomb Radiocarbon Analyses: 
\code{\link{bias.plot.BR}()},
\code{\link{calc.bias}()},
\code{\link{calc.h}()},
\code{\link{closest.pt}()},
\code{\link{get.h.from.bias}()},
\code{\link{get.nearest}()},
\code{\link{get.ref.line}()}
}
\concept{Bomb Radiocarbon Analyses}
