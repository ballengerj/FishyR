% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/completeFun.R
\name{complete.data}
\alias{complete.data}
\title{Extract Complete Cases from a Data Frame}
\usage{
complete.data(data, Cols)
}
\arguments{
\item{data}{Data frame}

\item{Cols}{vector, which may represent a character vector or numeric vector,
of columns for which we want to determine the complete cases
(\code{\link{complete.cases}}) for}
}
\value{
A restricted data frame that contains all complete cases for the columns specified.
}
\description{
\code{complete.data} allows one to extract complete cases (i.e., no NA values),
based on a potential subset of the total columns in the data frame, from a
data frame.
}
\examples{
(df <- data.frame(w = c(1, 1, NA, 2, 3, NA), x = c(2, 6, NA, 7, 8, 10),
y = c(2, NA, 10, 5, 6, 1), z = c(NA, 1, 3, NA, 4, 5))) # Original Data
complete.data(df, Cols = c("w", "x")) # Using named columns
complete.data(df, Cols = c(1, 2))
}
