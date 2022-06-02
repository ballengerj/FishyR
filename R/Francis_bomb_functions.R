#' Find Nearest Point
#'
#' \code{get.nearest} for each point (x,y), returns the x-coordinate of the
#'    point on the line that is closest to it in the x-direction.  The lie is
#'    a linear spline defined by joining a series of poitns by straight lines.
#' @param line a list, with components x and y (each an n-vector),
#'    containing the series of n points definign the line
#' @param x,y - m-vectors defining m points
#' @return Returns a vecto of lenght m.  This has value NA for any point whose
#'    y value lies outside the range of the y values in the line.
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' ref <- get.ref.line(snapper.red$C14year, snapper.red$C14)
#' get.nearest(ref, bluenose.red$C14year, bluenose.red$C14)
#' @family Bomb Radiocarbon Analyses
#' @export
get.nearest <- function(line, x, y) {
    npt <- length(x)
    xclose <- rep(NA, npt)
    line.len <- length(line$x)
    for(i in 1:npt) {
        above.pt <- ifelse(line$y < y[i], 0, 1)
        i1 <- (1:(line.len - 1))[diff(above.pt) != 0]
        i2 <- i1 + 1
        xposs <- line$x[i1] +(y[i] - line$y[i1]) * (line$x[i2] - line$x[i1]) /
            (line$y[i2] - line$y[i1])
        if(length(xposs) > 1) {
            closest <- order(abs(xposs - x[i]))[1]
            xclose[i] <- xposs[closest]
        }
        else if(length(xposs) == 1) xclose[i] <- xposs
    }
    xclose
}
#-------------------------------------------------------------------------------
#' Isotonic/Montone Regression Best Fit
#'
#' \code{get.ref.line} generates a list containing \code{\link{isoreg}} function
#'    model fits to \code{x} and \code{y} data.  The function 1st checks to see
#'    if multiple \code{y} values are present for a given \code{x}. If there
#'    are, it determines the \code{\link{mean}} \code{y} for each variable of
#'    \code{x}" before using the function \code{\link{isoreg}} to compute the
#'    isotonic (montonely increasing nonparametric) least squares regression
#'    which is piecewise constant that  bests fits the data.
#' @param x numeric vector representing the \code{x}-coordinates of the
#'    regression points.
#' @param y numeric vector representing the \code{y}-coordinates of the
#'    regression points.
#' @return Object of class \code{\link{list}} containing the results of the
#'    \code{link{isoreg}} function.  The \code{\link{list}} will be of length
#'    two:
#'    \describe{
#'      \item{x}{original (constructed) abscissa values \code{x}, sorted in
#'          ascending order}
#'      \item{y}{fitted vlaues correspondinb to \code{ordered x} values}
#'      }
#' @examples
#' set.seed(1234)
#' x <- 1:100
#' y <- x + rnorm(100, 0, 5)
#' get.ref.line(x, y)
#' @seealso \code{\link{isoreg}}
#' @family Bomb Radiocarbon Analyses
#' @export
get.ref.line <- function(x, y) {
    if(length(unique(x)) < length(x)) {
        y <- tapply(y, x, mean)
        x <- tapply(x, x, mean)
    }
    tmp <- isoreg(x, y)
    list(x = sort(tmp$x), y = tmp$yf)
}
#-------------------------------------------------------------------------------
#' Calculate the Median Horizontal Distance Between Linear Spline and Test
#'   Points
#'
#' \code{calc.h} calculates the median horizontal distance between
#'    the linear spline line and the set of points defined by test, where the
#'    horizontal distance is measured in standard errors.  Points whose
#'    y-values lie outside the range of the y-values for the line are excluded
#'    from teh calculation (because their horizontal distance from the line is
#'    undefined).
#' @param line a list, with components x and y (each an n-vector),
#'    containing the series of n points defining the line
#' @param test dataframe representing test data set as descrivbed in fucntion
#'    \code{\link[FishyR]{calc.bombcarbon.bias}}
#' @param age.err,test.x,test.y character strings representing the names of the
#'    columns in \code{test} containing the aging error estimates, the
#'    x-variable, and the y-variable, respectively.
#' @return Returns a numeric value representing the median horizontal distance
#'    from the reference line give the test data set.
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' ref <- get.ref.line(snapper.red$C14year, snapper.red$C14)
#' calc.h(line = ref, test = bluenose.red, age.err = "age.err", test.x =
#'   "C14year", test.y = "C14")
#' @family Bomb Radiocarbon Analyses
#' @export
calc.h <-
    function(line, test, age.err, test.x, test.y){
        yref <- get.nearest(line, test[, test.x],test[, test.y])
        median((test[, test.x] - yref) / test[, age.err], na.rm = TRUE)
    }
#-------------------------------------------------------------------------------
#' Determine Closest Point
#'
#' \code{closest.pt} finds, for each sample point in \code{samp}, the point on
#'    the reference curve (continuous linear spline defined by ref) that is
#'    'closest' to it.  Distance (and thus the definition of 'closest') is
#'    defined, differently for each sample point, by dividing the x and y
#'    values by the sampling errors for that point.
#' @param ref list, data frame, or matrix containing x- and y-coordinates for a
#'    reference line one wishes to compare sample points to.  Typically this
#'    will be the output from a call to \code{\link[FishyR]{get.ref.line}}.
#' @param samp data frame, with a minimum of 4 columns, representing your
#'    'test' data set.  Individual columns should represent x-cordinates,
#'    y-cordinates, aging error, and y error.  Typically x will represnt
#'    sample year or birth year and y will rerpresent either D14C or F14C.
#' @param ref.x character string giving the name of the x-value in the provided
#'    \code{ref} data set.
#' @param ref.y character string giving the name of the y-value in the provided
#'    \code{ref} data set.
#' @param samp.x character string giving the name of the x-value in the
#'    provided \code{samp} data set.
#' @param samp.y character string giving the name of the y-value in the
#'    provided \code{samp} data set.
#' @param age.err character string giving the name of the column representing
#'    aging error provided in the \code{samp} data set.
#' @param y.err character string giving the name of the column representing
#'    the y-variable error provided in the \code{samp} data set.
#' @param plot logical indicating whether to generate plot
#' @param incl.dist logical indicating whether to include distance calculations
#'    in resultant output
#' @return A dataframe with components (each with length nrow(samp)) with
#'    components X, Y - x & y coordinates of closest points - and optionally
#'    (if incl.dist = TRUE) dist - signed distance from reference curve.
#'    Optionally plots (if plot = TRUE) the data, including lines joining each
#'    sample point to its closest point.
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' ref <- get.ref.line(snapper.red$C14year, snapper.red$C14)
#' closest.pt(ref = ref, samp = bluenose.red, samp.x = "C14year",
#'   samp.y = "C14", age.err = "age.err", y.err = "C14.err", plot = TRUE,
#'   incl.dist = TRUE)
#' @family Bomb Radiocarbon Analyses
#' @export
closest.pt <- function(ref, samp, ref.x = "x", ref.y = "y",
                       samp.x = "Sample.Year", samp.y = "F14C",
                       age.err = "SE.Age", y.err = "F14C.SD", plotit = FALSE,
                       incl.dist = FALSE) {
    ## The points in ref may be defined by components ref$x, ref$y
    ## or ref$C14year, ref$C14
    nref <- length(ref[[ref.y]])
    nsamp <- length(samp[, samp.y])

    closest <- list(rep(NA, nsamp), rep(NA, nsamp), rep(Inf, nsamp))
    names(closest) <- c(samp.x, samp.y, "dist")
    for(j in 1:nsamp) {
        xref <- ref[[ref.x]] / samp[, age.err][j]
        yref <- ref[[ref.y]] / samp[, y.err][j]
        xsamp <- samp[, samp.x][j] / samp[, age.err][j]
        ysamp <- samp[, samp.y][j] / samp[, y.err][j]
        for(i in 1:(nref - 1)) {
            m <- (yref[i + 1] - yref[i]) / (xref[i + 1] - xref[i])
            if(m == Inf)# ith segment of ref line vertical
                intcpt <- c(ysamp, yref[i], yref[i + 1])
            else if(m == 0)# ith segment of ref line horizontal
                intcpt <- c(xsamp, xref[i], xref[i + 1])
            else{ # ith segment of ref line neither horizontal or vertical
                mprime <- -1 / m
                intcpt <- c(ysamp - mprime * xsamp, yref[i] - mprime * xref[i],
                            yref[i + 1] - mprime * xref[i + 1])
            }
            if((intcpt[1] - intcpt[2]) * (intcpt[1] - intcpt[3]) < 0){
                if(m == Inf){
                    xclose <- xref[i]
                    yclose <- ysamp
                }
                else if(m == 0){
                    xclose <- xsamp
                    yclose <- yref[i]
                }
                else{
                    xclose <- ((ysamp - mprime * xsamp)-
                                   (yref[i] -m * xref[i])) / (m - mprime)
                    yclose <- ysamp + mprime * (xclose - xsamp)
                }
            }
            else if((intcpt[2] == intcpt[1]) |
                    (intcpt[2] - intcpt[1]) * (intcpt[2] - intcpt[3]) < 0){
                xclose <- xref[i]
                yclose <- yref[i]
            }
            else{
                xclose <- xref[i + 1]
                yclose <- yref[i + 1]
            }
            ds <- sqrt((xsamp - xclose) ^ 2 + (ysamp - yclose) ^ 2) *
                ifelse(xsamp < xclose, -1, 1)
            if(abs(ds) < abs(closest$dist[j])){
                closest[[samp.x]][j] <- xclose * samp[, age.err][j]
                closest[[samp.y]][j] <- yclose * samp[, y.err][j]
                closest[["dist"]][j] <- ds
            }
        }
    }
    if(plotit){
        plot(ref[[ref.x]], ref[[ref.y]], type = 'b', pch = '+', xlab = '',
             ylab = '', xlim = range(c(samp[, samp.x],ref[[ref.x]])),
             ylim = range(c(samp[, samp.y], ref[[ref.y]])))
        points(samp[, samp.x], samp[, samp.y], col = 2, pch = 'x')
        sel <- closest$dist < 0
        if(any(sel))segments(samp[, samp.x][sel], samp[, samp.y][sel],
                             closest[[samp.x]][sel], closest[[samp.y]][sel],
                             col = 3, lty = 2)
        if(any(!sel))segments(samp[, samp.x][!sel], samp[, samp.y][!sel],
                              closest[[samp.x]][!sel], closest[[samp.y]][!sel],
                              col = 4, lty = 2)
    }
    if(incl.dist) as.data.frame(closest) else as.data.frame(closest[1:2])
}
#-------------------------------------------------------------------------------
#' Simulate h-Statistic Distribution
#'
#' \code{sim.hdist} Uses simulation to generate the distribution of the
#'    statistic h, which measures how well we can expect test and reference
#'    bomb radiocarbon data sets to agree, given the extent of sampling error
#'    associated with each.  A positive (or negative) value of h implies that
#'    the test data is, on average, to the right (or left) of the reference
#'    line.
#' @param ref reference data used to develop a reference line; a data frame in
#'    which each row corresponds to one reference datum.
#' @param test data frame containing the data from the species whose ages are
#'    being validated.
#' @param nsim number of simulated data sets.
#' @param diag integer, with possible values of 0, 1 and 2 that controls what,
#'    if any, diagnostic plots are plotted:
#'    \describe{
#'      \item{0}{no diagnostic plots are constructed}
#'      \item{1}{plot actual data (refrence and test) and also the data
#'         considered 'true' for the simulations; does no simulations; returns
#'         the 'true test data.}
#'      \item{2}{plots one panel per simulation, 16 panels per page, with each
#'         panel showing the simulated (solid line) and true (broken line)
#'         reference lines and the simulated test data, with the estimated
#'         h-statistic given above the panel.  (Note - use of this option is
#'         not sensible when nsim is 'large')}
#'      }
#' @return A dataframe with components (each with length nrow(samp)) with
#'    components X, Y - x & y coordinates of closest points - and optionally
#'    (if incl.dist = TRUE) dist - signed distance from reference curve.
#'    Optionally plots (if plot = TRUE) the data, including lines joining each
#'    sample point to its closest point.
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' ref <- get.ref.line(snapper.red$C14year, snapper.red$C14)
#' sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 16, diag = 0)
#' sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 16, diag = 1)
#' sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 16, diag = 2)
#' @family Bomb Radiocarbon Analyses
#' @export
sim.hdist <- function(ref, test, nsim = 5000, diag = 0, age = "age",
                      age.err = "age.err", cy = "catch.year", y = "C14",
                      y.err = "C14.err", s.age = "samp.age", b.yr = "C14year")
{
    nref <- nrow(ref)
    ntest <- nrow(test)
    ## Set up test data set
    if(is.na(match(s.age, names(test)))) test[, s.age] <- rep(0, ntest)
    bad.rows <- (1:ntest)[test[, b.yr] !=
                              (test[, cy] - (test[, age] - test[, s.age]))]
    if(length(bad.rows) > 0)
        stop(paste('Y-varaible wrongly calculated at rows',
                   paste(bad.rows, collapse = ' '), 'of test data\n'))
    ## Set up reference data set and calculate true.ref.line
    if(is.na(match(s.age, names(ref)))) ref[, s.age] <- rep(0, nref)
    bad.rows <- (1:nref)[ref[, b.yr] !=
                             (ref[, cy] - (ref[, age] -ref[, s.age]))]
    if(length(bad.rows) > 0)
        stop(paste('Y-variable wrongly calculated at rows',
                   paste(bad.rows, collapse = ' '), 'of ref data\n'))
    true.ref.line <- get.ref.line(ref[, b.yr], ref[, y])
    ## Calculate true.ref and true.test by finding the point on
    ## true.ref.line that is 'closest' to each reference and test data
    ## point, where 'closest' is in terms of s.e.s in both x- and
    ## y-directions
    true.ref <- ref
    true.ref[c(b.yr, y)] <- closest.pt(ref = true.ref.line, samp = ref,
                                       samp.x = b.yr, samp.y = y,
                                       age.err = age.err, y.err = y.err)
    true.ref[, age] <- true.ref[, cy] + true.ref[, s.age] - true.ref[, b.yr]
    true.test <- test
    true.test[c(b.yr, y)] <- closest.pt(ref = true.ref.line, samp = test,
                                        samp.x = b.yr, samp.y = y,
                                        age.err = age.err, y.err = y.err)
    true.test[, age] <- true.test[, cy] + true.test[, s.age] -
        true.test[, b.yr]
    ## Diagnostic plot for diag==1: actual and true data
    if(diag == 1){
        palette(c('black','red'))
        Range <- function(x) range(x[!is.na(x)])
        xlm <- Range(c(test[, b.yr], true.test[, b.yr], ref[, b.yr],
                       true.ref[, b.yr], true.ref.line$x))
        ylm <- Range(c(test[, y], true.test[, y], ref[, y], true.ref[, y],
                       true.ref.line$y))
        plot(ref[, b.yr], ref[, y], pch = 'x', xlab = '', ylab = '',
             xlim = xlm, ylim = ylm)
        lines(true.ref.line,col=2)
        points(true.ref[, b.yr], true.ref[, y], pch = 'x', col = 2)
        points(test[, b.yr], test[, y], pch = 'o')
        points(true.test[, b.yr], true.test[, y], pch = 'o', col = 2)
        pr <- par()
        for(i in 1:4)
            text(xlm[1], ylm[2] - c(0.5, 1.5, 3.5, 4.5)[i] * pr$cxy[2],
                 c("'True'", "Actual", "x - Reference", "o - Test")[i],
                 adj = 0, col = ifelse(i==1, 2, 1))
        return(true.test)
    }
    if(diag == 2){
        par(mfrow = c(4, 4), mar = c(1, 2, 2, 1) + .1, mgp = c(0, 0.5, 0),
            las = 1)
    }
    h <- rep(0, nsim)
    for(isim in 1:nsim){  ## simulation loop
        sim.ref.age <- round(rnorm(nref, true.ref[, age], ref[, age.err]))
        sim.ref.age[sim.ref.age < 0] <- 0
        sim.ref <- list(Y = rnorm(nref, true.ref[, y], ref[, y.err]),
                        Year = true.ref[, cy] -
                            (sim.ref.age - ref[, s.age]))
        sim.ref.line <- get.ref.line(sim.ref$Year, sim.ref$Y)
        sim.test.age <- round(rnorm(ntest, true.test[, age],
                                    test[, age.err]))
        sim.test.age[sim.test.age < 0] <- 0
        sim.test <- data.frame(Y = rnorm(ntest, true.test[, y],
                                         test[, y.err]),
                               Year = true.test[, cy] -
                                   (sim.test.age-test[, s.age]),
                               age.err = test[, age.err])
        h[isim] <- calc.h(line = sim.ref.line, test = sim.test,
                          age.err = "age.err", test.x = "Year",
                          test.y = "Y")
        if(diag == 2){
            plot(sim.test$Year, sim.test$Y, pch = '+', xlab = '',
                 ylab = '', xlim = range(c(sim.ref$Year,sim.test$Year,
                                           sim.ref.line$x)),
                 ylim=range(c(sim.ref$Y, sim.test$Y, sim.ref.line$y)))
            lines(sim.ref.line)
            lines(true.ref.line,lty=2)
            mtext(paste('h =', round(h[isim], 2)))
            if((isim < nsim) & (isim %% 16) == 0) locator(1)
        }
    }
    h
}
#-------------------------------------------------------------------------------
#' Calculate Relationship Between "h" and Bias
#'
#' \code{get.h.from.bias} calculates the relationship between \code{h} and
#'    bias, returning a list with components \code{hval} and
#'    \code{bias.conf.int}
#' @param ref,test dataframes representing reference and test data sets,
#'    respectively, as descrivbed in fucntion \code{\link[FishyR]{calc.bias}}
#' @param hdist simulated set of h values as generaed by function
#'    \code{\link[FishyR]{sim.hdist}}
#' @param bias.step,bias.range the function calculates \code{h} for each value
#'    of bias in \code{seq(bias.range[1], bias.range[2], bias.step)}.  If bias
#'    .range == NULL, then bias.range is initially set to c(-100, 100) but then
#'    is truncated to be just wide enough to include the 95% confidence
#'    interval for h.
#' @param ref.x,ref.y character string giving the name of the x and y-values in
#'    the provided \code{ref} data set, respectively.
#' @param age,age.err,cy,s.age,test.x,test.y character strings giving the
#'    names of the columns representing ages, age error estimates, capture
#'    year, sample ages, x-variable, and y-variable in the \code{test}
#'    dataframe.
#' @return A list.
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' ref <- get.ref.line(snapper.red$C14year, snapper.red$C14)
#' hdist <- sim.hdist(ref = snapper.red, test = bluenose.red, nsim = 5000,
#'    diag = 0)
#' get.h.from.bias(ref = ref, test = bluenose.red, hdist = hdist, bias.step = 5,
#'    bias.range = NULL, age = "age", cy = "catch.year", age.err = "age.err",
#'    s.age = "samp.age", test.x = "C14Year", test.y = "C14")
#' @family Bomb Radiocarbon Analyses
#' @export
get.h.from.bias <- function(ref, test, hdist, bias.step = 5, bias.range = NULL,
                            ref.x = "x", ref.y = "y", age = "Fish.Age",
                            cy = "Catch.Year", age.err = "SE.Age",
                            s.age = "Sample.Age", test.x = "Sample.Year",
                            test.y = "F14C") {
    h.conf.int <- quantile(hdist,c(0.025,0.975))
    if(is.null(bias.range)){
        trunc.bias <- TRUE
        bias.range <- c(-100,100)
    }
    else trunc.bias <- FALSE
    biasval <- seq(bias.range[1], bias.range[2], bias.step)
    hval <- rep(NA, length(biasval))
    done <- FALSE
    bs <- 0
    bs.inc <- bias.step
    ref.line <- get.ref.line(ref[[ref.x]], ref[[ref.y]])
    while(bs <= max(biasval) & bs >= min(biasval) & !done){
        test$age.corrected <- round(test[, age] / (1 + bs / 100))
        test$new.C14year <- test[, cy] - (test$age.corrected - test[, s.age])
        hval[biasval == bs] <-
            calc.h(ref.line, test = test, age.err = age.err,
                   test.x = "new.C14year", test.y = test.y)
        if(bs > 0 & ((trunc.bias & hval[biasval == bs] > h.conf.int['97.5%']) |
                     bs == bias.range[2])){
            bs <- 0
            bs.inc <- -bias.step
        }
        else if(bs < 0 &
                (trunc.bias & hval[biasval == bs] < h.conf.int['2.5%'])) {
            done <- TRUE
        }
        bs <- bs + bs.inc
    }
    sel <- !is.na(hval)
    biasval <- biasval[sel]
    hval <- hval[sel]
    names(hval) <- paste(biasval)

    ## Deal with repeated values of hval, if any
    if(length(unique(hval)) < length(hval)) {
        biasval <- tapply(biasval, hval, mean)
        hval <- tapply(hval, hval, mean)
        names(hval) <- paste(biasval)
    }
    bias.conf.int <- round(approx(hval, biasval, h.conf.int)$y)
    names(bias.conf.int) <- names(h.conf.int)
    list(hval = hval, bias.conf.int = bias.conf.int)
}
#-------------------------------------------------------------------------------
#' Calculate Relationship Between "h" and Bias
#'
#' \code{calc.bias} calculates a 95% confidence interval for the percentage
#'    bias in zone-count ages given appropriate bomb-carbon data.
#' @param ref,test data frames representing reference and test data sets,
#'    respectively.  Each row corresponds to one ref (test) datum and the
#'    columns represent data needed to perform the bomb carbon analyses.
#'    One column should represent each of the following input data in each data
#'    frame:
#'    \describe{
#'      \item{age}{estimated age}
#'      \item{age.err}{standard deviation of estimated age}
#'      \item{s.age}{estimated age of the fish when the material taken from
#'         the otolith/aging structure for bomb carbon analysis was
#'         incorporated into the otolith (if this column is omitted it will
#'         be assumed that \code{s.age = 0} for all samples)}
#'      \item{cy}{year of sampling, i.e. catch year}
#'      \item{x}{inferred year in which material taken from the otolith/aging
#'         structure for bomb carbon analysis was incorporated into the
#'         otolith (A warning will be output if \code{x} is inconsistent with
#'         other variables; should be \code{x = cy - (age - s.age)})}
#'      \item{y}{bomb carbon value (this will generally either be in C14 or
#'         fraction modern values)}
#'      \item{y.err}{standard deviation of bomb carbon value}
#'      }
#'    Note, the columns in both the \code{ref} and \code{test} data frames
#'    representing each variable should be named the same.
#' @param nsim number of simulated data sets used in generating the
#'    distribution of h values
#' @param plot.type type of plot to produce (if 0, no plot is made).  See
#'    function \code{\link[FishyR]{plot.bias}} for available plot types
#' @param age,age.err,s.age,cy,x,y,y.err Character strings representing the
#'    column names for each of the required elements (see above) in the ref and
#'    test data frames.
#' @param bias.step,bias.range See \code{\link[FishyR]{get.h.from.bias}} for
#'    description
#' @return A named list of length 4 with the following elements:
#'    \describe{
#'      \item{hdist}{obersved h-statistic value for each of the simulations}
#'      \item{h.conf.int}{estimated 95% confidence interval for the h-statistic
#'         based on the simulations}
#'      \item{hval}{observed h-statistic at bias levels defined by the range of
#'         biases considered}
#'      \item{bias.conf.int}{95% confidence interval about aging bias}
#'      }
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' calc.bias(snapper.red, bluenose.red, nsim = 5000, plot.type = 0,
#'   x = "C14year", y = "C14", age = "age", cy = "catch.year",
#'   age.err = "age.err", s.age = "samp.age", y.err = "C14.err")
#' calc.bias(snapper.red, bluenose.red, nsim = 5000, plot.type = 1,
#'   x = "C14year", y = "C14", age = "age", cy = "catch.year",
#'   age.err = "age.err", s.age = "samp.age", y.err = "C14.err")
#' calc.bias(snapper.red, bluenose.red, nsim = 5000, plot.type = 2,
#'   x = "C14year", y = "C14", age = "age", cy = "catch.year",
#'   age.err = "age.err", s.age = "samp.age", y.err = "C14.err")
#' @family Bomb Radiocarbon Analyses
#' @export
calc.bias <- function(ref, test, nsim = 5000, plot.type = 0, x, y, age, cy,
                      age.err, s.age, y.err, bias.step = 5,
                      bias.range = NULL) {
    nref <- nrow(ref)
    ntest <- nrow(test)
    ## Set up test data set
    if(is.na(match(s.age, names(test))))
        test[, s.age] <- rep(0, ntest)
    bad.rows <- (1:ntest)[test[, x] !=
                              (test[, cy] - (test[, age] - test[, s.age]))]
    if(length(bad.rows) > 0)
        stop(paste('Sample year wrongly calculated at rows',
                   paste(bad.rows, collapse = ' '), 'of test data\n'))
    ## Set up reference data set and calculate ref.line
    if(is.na(match(s.age, names(ref))))
        ref[, s.age] <- rep(0, nref)
    bad.rows <- (1:nref)[ref[, x] !=
                             (ref[, cy] - (ref[, age] - ref[, s.age]))]
    if(length(bad.rows) > 0)
        stop(paste('Sample year wrongly calculated at rows',
                   paste(bad.rows, collapse = ' '), 'of ref data\n'))
    ## Generate simulated distribution of h and associated confidence int.
    hdist <- sim.hdist(ref, test, nsim = nsim, age = age,
                       age.err = age.err, cy = cy, y = y, y.err = y.err,
                       s.age = s.age, b.yr = x)
    h.conf.int <- quantile(hdist, c(0.025, 0.975))
    ## Calculate relationship between h and b to cover range of h.conf.int
    tmp <- get.h.from.bias(ref, test, hdist, bias.step = bias.step,
                           bias.range = bias.range, ref.x = x, ref.y = y,
                           age = age, cy = cy, age.err = age.err,
                           s.age = s.age, test.x = x, test.y = y)
    bias.out <- list(hdist = hdist, h.conf.int = h.conf.int,
                     hval = tmp$hval, bias.conf.int = tmp$bias.conf.int)
    #plot.bias(bias.out, plot.type)
    bias.out
}











#-------------------------------------------------------------------------------
#' Bomb Radiocarbon Plots
#'
#' \code{plot.bias} plots output from function \code{\link[FishyR]{calc.bias}}
#' @param bias.out list containing output from function
#'    \code{\link[FishyR]{calc.bias}}
#' @param type type of plot
#'    \describe{
#'      \item{0}{no plot}
#'      \item{1}{the simulated h distribution is plotted, along with arrows
#'         showing the h values for different values of percentage bias and
#'         broken lines showing the 95% confidence interval for h}
#'      \item{2}{a cumulative distribution of percentage bias is plotted, with
#'         broken lines showing the 95% confidence interval for bias}
#'      }
#' @return A plot.
#' @examples
#' data(snapper)
#' snapper.red <- snapper[snapper$C14year >= 1955 & snapper$C14year <= 1972, ]
#' data(bluenose)
#' bluenose.red <- bluenose[bluenose$C14year >= 1955 &
#'   bluenose$C14year <= 1972, ]
#' bluenose.red <- bluenose.red[bluenose.red$C14 >= -506 &
#'   bluenose.red$C14 <= 94.0, ]
#' out <- calc.bias(snapper.red, bluenose.red, nsim = 5000, plot.type = 0,
#'   x = "C14year", y = "C14", age = "age", cy = "catch.year",
#'   age.err = "age.err", s.age = "samp.age", y.err = "C14.err",
#'   bias.range = c(-50, 50), bias.step = 1)
#' #plot.bias(out, type = 1)
#' #plot.bias(out, type = 2)
#' @family Bomb Radiocarbon Analyses
#' @export
plot.bias <- function(bias.out, type) {
    hval <- bias.out$hval
    nbias <- length(hval)
    hdist <- bias.out$hdist
    biaslab <- names(hval)
    biasval <- as.numeric(biaslab)
    if(type == 1){
        hdist.line <- density(hdist)
        ymax <- 1.05 * max(hdist.line$y)
        par(mar = c(3, 3, 1, 1) + .1, mgp = c(2, 0.5, 0), las = 1)
        plot(hdist.line, type = 'l', ylim = c(0, ymax), yaxs = 'i', xlab = 'h',
             ylab = '', main = '', yaxt = 'n')
        arrows(hval, rep(0.5 * ymax, nbias), hval, rep(0, nbias), length = 0.1)
        text(hval, rep(0.5 * ymax, nbias), biaslab, srt = 90, adj = 0)
        abline(v = bias.out$h.conf.int, lty = 2)
    }
    else if(type==2){
        pval <- rep(0, nbias)
        for(i in 1:nbias)pval[i] <- mean(hdist <= hval[i])
        par(mar = c(3, 3, 1, 1) + .1, mgp = c(2, 0.5, 0), las = 1)
        plot(biasval, pval, type = 'b', pch = 'x', xlab = '% bias',
             ylab = 'Cumulative probability', ylim = c(0, 1), yaxs = 'i')
        abline(v = bias.out$bias.conf.int, lty = 2)
    }
}
