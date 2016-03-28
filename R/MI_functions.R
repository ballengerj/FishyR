#' General Random Normal Noise About a Vector
#'
#' \code{rnorm.noise} generates n random normal deviates about a mean value with a
#'    specified standard deviation.
#'
#' @param Data frame with two columns.  1st column identifies the group and the
#'    2nd identifies the mean for that group
#' @param n integer specifying the number of random normal deviates to generate
#'    for each value of x.
#' @param sd standard deviation of the random normal deviates to be created.
#' @return An object of class 'data frame' that contains two columns, one
#'    representing the "group" defined by the \code{x} parameter and one
#'    representing individual random normal deviates about the given "group".
#' @examples
#' # Example
#' rnorm.noise(Data=c(1, 10), n = 10, sd = 0.5)
#' # Example 2
#' # Signal data
#' x <- seq(0, 4 * pi, length = 241)
#' y <- scales::rescale(cos(seq(0, 4 * pi, length = 241)), to=c(0.2, 0.8))
#' plot(y ~ seq(0, 24, 0.1), type = "l", xaxt = "n", xlab="Month")
#' axis(1, at = seq(0, 24, 1), labels = c (rep(month.abb, 2),"Jan"))
#' # Noise Data
#' noise <- apply(data.frame(x=seq(0, 24, 0.1), y = y), 1, rnorm.noise, n = 10,
#'    sd = 0.1)
#' noise.df <- Reduce(function(...) merge(..., all = TRUE), noise)
#' # Plot Noise
#' plot(noise.df[, 2] ~ noise.df[, 1], pch = 20, ylab = "Proportion",
#'    xlab = "Month", xaxt="n")
#' axis(1, at = seq(0, 24, 1), labels = c (rep(month.abb, 2),"Jan"))
#' @export
rnorm.noise <- function(Data, n, sd) {
  data.frame(x = rep(Data[1], n), y = rnorm(n, mean = Data[2], sd = sd))
}

# ------------------------------------------------------------------------------
#' Marginal Increment Analysis Plot
#'
#' \code{MIA} creates a marginal increment analysis plot.
#' @param Data Data frame containing, at a minimum, two columns.  The first
#'    column should represnt the x-variable.  The 2nd column should represent
#'    the y-variable.  Other columns may be present.
#' @param Statistics Data frame containing, at a minimum, four colummns.  The
#'    first column should represent the x-variable.  The 2nd column should
#'    represent the y-variable, typically the mean of y at a given x.  The 3rd
#'    column should represent a measure of variability, typically the confidence
#'    interval of y given x.  The 4th column should represent the number of
#'    samples of y at a given level of x.  Other columns may be present in the
#'    data frame.
#' @param x.label Text label for the x-axis of the plot
#' @param pch.raw symbol to use for the raw x and y data found in \code{Data}
#' @param col.raw color to use for the raw x and y data found in \code{Data}
#' @param alpha.raw amount of transparency of the raw x and y data found in
#'    \code{Data}
#' @param cex.raw symbol size to use for the raw x and y data found in
#'    \code{Data}
#' @param pch.fit symbol to use for the summarized marginal increment data
#'    found in \code{Statistics}
#' @param col.fit color to use for the summarized marginal increment data
#'    found in \code{Statistics}
#' @param cex.fit symbol size to use for the summarized marginal increment
#'    data found in \code{Statistics}
#' @param col.text color of text to use for depiciting the sample size at each x
#'    variable level
#' @param cex.text text size to use for depiciting the sample size at each x
#'    variable level
#'    bias plot between Reader 2 and Reader 1 ages.
#' @return A \code{\link{xyplot}} using the package lattice depicting the
#'    results of a marginal increment analysis.
#' @seealso \code{\link{MIA.Edge}}
#' @examples
#' # Example with Sheepshead Data
#' data(Sheepshead)
#' Data <- data.frame(Month = Sheepshead$Month,
#'                    Perc.Comp = Sheepshead$MI/Sheepshead$Prev.Inc)
#' Data$Perc.Comp <- with(Data, ifelse(Perc.Comp >1, 1, Perc.Comp))
#' Statistics <- with(Data, aggregate(Perc.Comp ~ Month, FUN = describe,
#'    digits = 4 ))
#' Statistics <- cbind(Statistics[-ncol(Statistics)],
#'    Statistics[[ncol(Statistics)]])
#' MIA(Data = Data, Statistics = Statistics[, c(1, 3, 6, 2)])
#' @export
MIA <- function(Data, Statistics, x.label = "Month", pch.raw = 20,
                col.raw = "#0080ff", alpha.raw = 0.5, cex.raw = 0.8,
                pch.fit = 16, col.fit = "black", cex.fit = 1,
                col.text = "Black", cex.text = 1) {
    grid <- lattice::xyplot(Data[, 2] ~ Data[, 1],
                            xlab = x.label,
                            ylab = "Percent Completion",
                            type = c("g","p"),
                            pch = pch.raw,
                            col= col.raw,
                            alpha = alpha.raw,
                            cex = cex.raw)
    points <- lattice::xyplot(Statistics[, 2] ~ Statistics[, 1],
                              panel = function(x, y) {
                                  lattice::panel.xyplot(x, y, type = "o",
                                                        pch = pch.fit,
                                                        col = col.fit,
                                                        cex = cex.fit)
                                  errbar(x = Statistics[, 1],
                                         y = Statistics[, 2],
                                         error = Statistics[, 3],
                                         color = col.fit)
                                  lattice::panel.text(x,
                                                      y = max(Data[, 2]),
                                                      labels = Statistics[, 4],
                                                      col = col.text,
                                                      cex = cex.text)
                              })
    grid + latticeExtra::as.layer(points)
}

# ------------------------------------------------------------------------------
#' Marginal Increment Analysis Plot - Uses only Categorical Edge Codes
#'
#' \code{MIA.Edge} creates a marginal increment analysis plot.
#' @param x vector representing the grouping variable for the marginal
#'    increment analysis.
#' @param y vector representing the edge type for each value of "x". Must be
#'    the same length as "x".
#' @param x.label label for the x-axis
#' @param y.label label for the y-axis
#' @return A \code{\link{barchart}} using the package lattice depicting the
#'    proportional edge type distribution by unique "x".
#' @seealso \code{\link{MIA}}
#' @examples
#' # Example with Sheepshead Data
#' data(Sheepshead)
#' Month <- Sheepshead$Month
#' Edge <- with(Sheepshead, ifelse(MI/Prev.Inc>=0.6667, "Wide", "Narrow"))
#' MIA.Edge(x = Month, y = Edge)
#' @export
MIA.Edge <- function(x, y, x.label = "Month", y.label = "Proportion") {
  distr <- prop.table(table(x, y), 1)
  lattice::barchart(distr,
                    horizontal = FALSE,
                    auto.key = list(columns = dim(distr)[2],
                                    title = "Edge Type"),
                    ylab = y.label,
                    xlab = x.label)
}
