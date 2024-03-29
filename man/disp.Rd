% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Index_functions.R
\name{disp}
\alias{disp}
\title{Calculate the dispersion estimate from a GLM model
\code{disp}}
\usage{
disp(mod = NA)
}
\arguments{
\item{mod}{Name of a fitted model of various types from which Pearson
residuals and a residual degrees of freedom can be extracted from the
model object}
}
\description{
Calculate the dispersion estimate from a GLM model
\code{disp}
}
\examples{
X <- seq(1, 100)
Y <- 100 + 0.2 * X + rnorm(1, mean = 0, sd = 5)
mod <- glm(Y ~ X)
disp(mod) # Dispersion parameter estimated using functio
sigma(mod)^2 # Dispersion parameter estimated using the calculated sigma
# paramter of the model object
}
\seealso{
Other Model Evaluation: 
\code{\link{Index.Grid}()},
\code{\link{Index.Summary}()}
}
\concept{Model Evaluation}
