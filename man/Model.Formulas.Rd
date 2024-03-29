% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Index_functions.R
\name{Model.Formulas}
\alias{Model.Formulas}
\title{Develop formulas for a variety of different types of modles}
\usage{
Model.Formulas(
  Predictors,
  Response,
  model,
  offset = NULL,
  VarType,
  names,
  count.n,
  basis
)
}
\arguments{
\item{Predictors}{see help pages for \\code{\link{GLM.Formulas}} and
\code{\link{ZI.Formulas}} for description.}

\item{Response}{see help pages for \code{\link{GAM.Formulas}},
\code{\link{GLM.Formulas}}, and \code{\link{ZI.Formulas}} for description.}

\item{model}{character string specifiying the type of model formula to be
constructed.  Currently only "GAM", "GLM", and "ZI" are acceptable values,
with all others giving an error.}

\item{offset}{see help pages for \code{\link{GLM.Formulas}} and
\code{\link{ZI.Formulas}} for description.}

\item{VarType}{see help pages for \code{\link{GAM.Formulas}},
\code{\link{GLM.Formulas}}, and \code{\link{ZI.Formulas}} for description.}

\item{names}{see help pages for \code{\link{GAM.Formulas}},
\code{\link{GLM.Formulas}}, and \code{\link{ZI.Formulas}} for description.}

\item{count.n}{see help page for \code{\link{ZI.Formulas}} for description.}

\item{basis}{see help pages for \code{\link{GAM.Formulas}} for description.}
}
\value{
Character string representing a potential model formula.
}
\description{
\code{Model.Formulas} generates a character vector that can be used to
  represent formulas for use in individual GAM, GLM, or zero-inflation
  model runs.  In essence it is a wrapper for the individual functions
  \code{\link{GAM.Formulas}}, \code{\link{GLM.Formulas}}, and
  \code{\link{ZI.Formulas}}.  The intent is to quickly generate model
  formulas that incorporates information from several potential covariates
  and the type of covariate utilized.
}
\examples{
data(iris)
# GAM Model
form.GAM <- GAM.Formulas(Response = "Sepal.Length", names =
  c("Petal.Length", "Petal.Width", "Sepal.Width", "Species"), VarType =
  c(rep("C", 3), "F"), basis = "tp")
form.GAM
# GLM Model
form.GLM <- GLM.Formulas(Response = "Sepal.Length", Predictors =
c(1, 2, 0, 3), names = c("Petal.Length", "Petal.Width", "Sepal.Width",
"Species"), VarType = c(rep("C", 3), "F"))
form.GLM
# Zero-Inflation Model
form.ZI <- ZI.Formulas(Response = "Sepal.Length", Predictors =
c(1, 2, 0, 3, 1, 0, 0, 0), names = rep(c("Petal.Length", "Petal.Width",
"Sepal.Width", "Species"), 2), VarType = rep(c(rep("C", 3), "F"), 2),
count.n = 4)
form.ZI
}
\seealso{
Other Model Formulas: 
\code{\link{GAM.Formulas}()},
\code{\link{GLM.Formulas}()},
\code{\link{ZI.Formulas}()}
}
\concept{Model Formulas}
