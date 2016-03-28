#' Root Mean Square Deviation
#'
#' \code{rmsd} is a function that allows the calculation of the
#'  \strong{root mean square deviation (RMSD)} for a vector of numbers.  This
#'  is also referred to as the '\emph{root mean squared error}' in the
#'  literature.
#' @param x numeric vector.
#' @param expected expectation, given the model, for each value of \strong{x}.
#' @return Based on the assumed underlying model, the estimated RMSD given the
#'  observed data in \code{x}.
#' @examples
#' # Examples
#' # Generated Data
#' set.seed(1234)
#' x <- rnorm(100)
#' mean(x)
#' rmsd(x, mean(x))
#' # Observed Data (Red Porgy data set)
#' data(RedPorgy)
#' x <- na.exclude(RedPorgy[, 1])
#' rmsd(x = x, expected = mean(x))
#' @family Model Error Functions
#' @export
rmsd <- function(x, expected) {
  sqrt(mean( (x - expected) ^ 2)) # Calculate root mean square deviation
}

# ------------------------------------------------------------------------------
#' Mean Absolute Error
#'
#' \code{mae} is a function that allows the calculation of the
#'  \strong{mean absolute error (MAE)} for a vector of numbers.
#' @param x numeric vector.
#' @param expected expectation, given the model, for each value of \strong{x}.
#' @return Based on the assumed underlying model, the estimated MAE given the
#'  observed data in \code{x}.
#' @examples
#' # Examples
#' set.seed(1234)
#' x <- rnorm(100)
#' mean(x)
#' mae(x, mean(x))
#' # Observed Data (Red Porgy data set)
#' data(RedPorgy)
#' x <- na.exclude(RedPorgy[, 1])
#' mae(x = x, expected = mean(x))
#' @family Model Error Functions
#' @export
mae <- function(x, expected) {
  mean( abs(x - expected)) # Calculate mean absolute error
}

# ------------------------------------------------------------------------------
#' Standard Error
#'
#' \code{SE} is a function that allows the calculation of the
#'  \strong{standard error} for a vector of numbers.
#' @param x numeric vector.
#' @return Standard error estimate for x.
#' @family Error Estimates
#' @examples
#' # Examples
#' set.seed(1234)
#' SE(rnorm(100))
#' # Observed Data (Red Porgy data set)
#' data(RedPorgy)
#' SE(x = RedPorgy$R1.Age)
#' @export
SE <- function(x) {
  sd(x, na.rm=TRUE) / sqrt(length(x[is.na(x) == FALSE])) # Calculate standard
} # error

# ------------------------------------------------------------------------------
#' Coefficient of Variation (CV)
#'
#' \code{CV} is a function that allows the calculation of the
#'  \strong{coefficient of variation} for a vector of numbers.
#' @param x numeric vector.
#' @return CV estimate for x.
#' @family Error Estimates
#' @examples
#' # Examples
#' set.seed(1234)
#' CV(rnorm(100))
#' # Observed Data (Red Porgy data set)
#' data(RedPorgy)
#' CV(x = RedPorgy$R1.Age)
#' @export
CV <- function(x) {
  abs(sd(na.exclude(x)) / mean(na.exclude(x))) # Calculate the CV
  # around the mean
}

# ------------------------------------------------------------------------------
#' Relative (Percent) Standard Error
#'
#' \code{RSE} is a function that allows the calculation of the
#'  \strong{relative standard error (RSE)} for a vector of numbers.  RSE is also known as \strong{Percent Standard Error (PSE)}.
#' @param x numeric vector.
#' @return RSE estimate for x.
#' @family Error Estimates
#' @examples
#' # Examples
#' set.seed(1234)
#' RSE(rnorm(100))
#' # Observed Data (Red Porgy data set)
#' data(RedPorgy)
#' RSE(x = RedPorgy$R1.Age)
#' @export
RSE <- function(x) {
  SE(na.exclude(x)) / mean(na.exclude(x)) # Calculate the relative (percent)
    # standard error for a sample x
}

# ------------------------------------------------------------------------------
#' Add Error Bars to Base or Lattice Plot
#'
#' \code{errbar} is a function that allows the rapid addition of error bars to
#' base or Lattice graphics.
#' @param x vector of x values corresponding to x-values of data plotted in a
#'  plot.  There should be paired y-value data for each of these x values. Must
#'  be the same length as y and error.
#' @param y vector of y values corresponding to y-values of data plotted in a
#'  plot.  There should be paired x-value data for each of these y values. Must
#'  be the same length as x and error.
#' @param error vector of "\strong{error}" (e.g., could represent, standard
#'  errors, confidence intervals, etc.) estimates for each of the x and y values
#'  provided in x and y.  Must be the same length as x and y.
#' @param color A specification for the plotting color of the error bars.  See
#'  description of \emph{Color Specificaiton} in the \code{\link{par}} for more
#'  information.
#' @param type Used to specify the type of plot the error bars are to be used
#'  in.  If type=1, error bars and constructed for use in base graphics plots.
#'  If type=2, error bars are constructed for use in lattice graphics plots.
#'  Any other specificaiton of type will cause an error.
#' @return This function is used to add error bars to an already created
#' graph.  The actual way to use this function will depend on the type of plot
#' being created.
#' @examples
#' # Generate Data
#' set.seed(1234)
#' x <- seq(1, 50, 1)
#' y <- (2 + 2 * x) * rnorm(50, 1, 0.2)
#' error <- rnorm(50, 2, 0.5)
#' # Example - Base Graphics
#' plot(y ~ x, xlab = "X Variable", ylab = "Y Variable")
#' errbar(x = x, y = y, error = error, type = 1)
#' # Example - Lattice Graphics
#' lattice::xyplot(y ~ x, xlab = "X Variable", ylab = "Y Variable",
#'  panel = function(x, y) {
#'    panel.xyplot(x, y)
#'    errbar(x = x, y = y, error = error, type = 2)
#'  })
#' @export
errbar <- function(x, y, error, color = "black", type = 2) {
  if (type == 1) {
    # Base Plots ------------------------------------------------------
    arrows(x, y - error, x, y + error, # Define type of "arrow" to draw in a
           angle = 90, length = 0.05, code = 3, col = color) # standard graphics
           # plot
  } else {
    if (type == 2) {
      # Lattice Plots -------------------------------------------------
      lattice::panel.arrows(x0 = x, y0 = y - error, x1 = x, y1 = y + error,
                            length = 0.05, angle = 90, # Define type of "arrow"
                            code = 3, col = color) # to draw in a lattice
        # graphics plot
    } else {
      if (type != 1 | type != 2) stop("Invalid type=? provided") # Indicates an
      # invalid "type" value was provided to the function
    }
  }
}

# ------------------------------------------------------------------------------
#' Calculate Smoothed Running Mean
#'
#' \code{errbar} is a function that allows the calculation of smoothed running
#'  means for a numeric vector.  It uses the \code{\link[zoo]{rollmean}}
#'  function
#'  from the package \code{\link{zoo}} to compute rolling means for the numeric
#'  vector and then the \code{\link{smoothEnds}} function from the package
#'  \code{\link{stats}} to smooth the end points of the resulting vector.
#' @param x numeric vector representing a series of observations.
#' @param window integer width of the rolling window.
#' @param align character specifying whether the index of the result should be
#'  left- or right-aligned or centered (default) compared to the rolling window
#'  of observations.
#' @return An object of the same class as x with the rolling mean with smoothed
#'  ends.
#' @seealso \code{\link[zoo]{rollmean}}, \code{\link{smoothEnds}}.
#' @examples
#' set.seed(1234)
#' x <- seq(1, 50, 1)
#' y <- (2 + 2 * x) * rnorm(50, 1, 0.2)
#' runmean(x = y, window=3)
#' @export
runmean <- function(x, window, align = "center") {
  ori <- x
  new <- zoo::rollmean(x, window, fill = NA, align = align)
  new[is.na(new)] <- ori[is.na(new)]
  smoothEnds(new, window)
}

# ------------------------------------------------------------------------------
#' Descriptive statistics
#'
#' \code{describe} constructs a named vector that contains various descriptive statistics of numeric data.
#' @param x Numeric vector.
#' @param confidence confidence level (default of 0.05) for resulting
#'    confidence interval calculation.
#' @param digits integer (defaults to 2) specifying the number of significant
#'    digits for the resulting descriptive statistics
#' @return A numeric vector containig the following information:
#'    \describe{
#'      \item{n}{sample size}
#'      \item{mean}{mean}
#'      \item{StDev}{standard deviation}
#'      \item{StErr}{standard error of the mean}
#'      \item{CI}{confidence interval width}
#'      \item{RSE}{relative standard error}
#'      \item{RMSE}{root mean square deviation}
#'      \item{MAE}{mean absolute error}
#'      }
#' @examples
#' # Example
#' data(BSB)
#' describe(BSB$TL, digits = 4)
#' # Example 2
#' with(BSB, aggregate(TL ~ Age, FUN = describe, digits = 4))
#' @export
describe <- function(x, confidence = 0.05, digits = 2) {
  N <- sum(is.na(x))
  n <- length(x) - N
  avg <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  se <- SE(na.exclude(x))
  cv <- CV(x)
  RSE <- RSE(x)
  RMSD <- rmsd(x, mean(x))
  MAE <- mae(x, mean(x))
  CI <- se * qt(1 - confidence / 2, n - 1)
  out <- c(n, avg, sd, se, CI, cv, RSE, RMSD, MAE)
  names(out) <- c("n", "mean", "StDev", "StErr", "CI", "CV","RSE", "RMSD",
                  "MAE")
  round(out, digits)
}
