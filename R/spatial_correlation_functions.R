#' Moran's I Test for Spatial Autocorrelation Among Residuals
#'
#' \code{Moran.I.calc} computes Moran's I autocorrelation coefficient of *resid*
#'   given a matrix of weights defined as the inverse of a distance matrix of
#'   spatial coordinates.
#' @param data numeric data frame with three columns.  First column should
#'   represent the "x"-coordinates, 2nd column represents the "y"-coordinates,
#'   and the 3rd column represents a numeric vector of model residuals
#' @return A list containing the elements \code{observed}, \code{expected},
#'   \code{sd}, & \code{p.value}
#' @seealso \code{\link[ape]{Moran.I}}
#' @export
Moran.I.calc <- function(data) {
    Dist <- as.matrix(dist(data[, c(1, 2)]))
    Dist.inv <- 1 / Dist
    diag(Dist.inv) <- 0
    Dist.inv[is.infinite(Dist.inv)] <- 0
    test <- ape::Moran.I(x = data[, 3], weight = Dist.inv, na.rm = TRUE)
    return(test)
}
# ------------------------------------------------------------------------------
#' Calculate maximum distance among spatial points
#'
#' \code{dist.max} computes the maximum distance among a set of spatial points
#' @param data numeric data frame with two columns.  First column should
#'   represent the "x"-coordinates and the 2nd column represents the
#'   "y"-coordinates
#' @param prop numeric variable
#' @return Character string representing a model name which represents the
#'   proposed model structure.
#' @export
dist.max <- function(data, prop) {
    Dist <- as.matrix(dist(data[, c(1, 2)]))
    diag(Dist) <- 0
    max(Dist) / 1000 * prop
}
# ------------------------------------------------------------------------------
#' Calculate variogram summary statistics and plots
#'
#' \code{Variogram_Plot} is a wrapper funtion that allows the efficient
#'   computation of numerous variagram statistics and variogram plots over a
#'   range of considered "plausible" distances.
#' @param data numeric data frame with two columns.  First column should
#'   represent the "x"-coordinates and the 2nd column represents the
#'   "y"-coordinates.  These columns must be titled "X" and "Y", respectively
#' @param mindist numeric value representing the minimum distance over which
#'   variogram statistics should be calculated
#' @param maxdist numeric value representing the maximum distance over which
#'   variogram statistics should be calculated
#' @param length numeric value representing how many points between
#'   \code{mindist} and \code{maxdist} variogram statistics should be calculated
#' @param type logical; if \code{TRUE}, use Cressie's robust variogram
#'   estimate; if \code{FALSE} use the classical method of moments variogram
#'   estimate.  See \code{\link[gstat]{variogram}}
#' @param map logical; See \code{\link[gstat]{variogram}}
#' @param cloud logical; See \code{\link[gstat]{variogram}}
#' @param bins numeric value used to calculate the \code{width} (See
#'   \code{\link[gstat]{variogram}}) of subsequent distance intervals
#' @return Data frame with results statistics and a series of variogram plots
#'   at different threshold distances
#' @export
Variogram_Plot <- function(data, mindist, maxdist, length, type, map = FALSE,
                           cloud = FALSE, bins = 100) {
    Vario.data <- data
    sp::coordinates(Vario.data) <- c("X","Y")
    pltList <- list()
    value <- NULL
    variable <- NULL
    Results <- data.frame(Max.Dist = NA, Min = NA, Max = NA, Mean = NA,
                           Median = NA, Model = NA, psill = NA, range = NA,
                           singular = NA)
    for (i in seq(mindist, maxdist, length.out = length)) {
        Vario1 <- gstat::variogram(Resid ~ 1,
                                   data = Vario.data,
                                   cutoff = i,
                                   width = i / bins,
                                   cressie = type,
                                   map = map,
                                   cloud=cloud)
        tmp <- data.frame(Max.Dist = round(i / 1000, 2),
                          Min = round(min(Vario1$np), 0),
                          Max = round(max(Vario1$np), 0),
                          Mean = round(mean(Vario1$np), 0),
                          Median = round(median(Vario1$np), 0))
        tmp <- tmp[rep(seq_len(nrow(tmp)), each = 3),]
        Vario1.fit.Exp <- gstat::fit.variogram(Vario1,
                                               gstat::vgm(1.5, "Exp", i))
        Vario1.fit.Sph <- gstat::fit.variogram(Vario1,
                                               gstat::vgm(1.5, "Sph", i))
        Vario1.fit.Gau <- gstat::fit.variogram(Vario1,
                                               gstat::vgm(1.5, "Gau", i))
        tmp2 <- rbind(Vario1.fit.Exp, Vario1.fit.Sph, Vario1.fit.Gau)
        test.Exp <- ifelse(attributes(Vario1.fit.Exp)[4] == TRUE,
                           "TRUE",
                           "FALSE")
        test.Sph <- ifelse(attributes(Vario1.fit.Sph)[4] == TRUE,
                           "TRUE",
                           "FALSE")
        test.Gau <- ifelse(attributes(Vario1.fit.Gau)[4] == TRUE,
                           "TRUE",
                           "FALSE")
        tmp$Model <- tmp2[, 1]
        tmp$psill <- round(tmp2[, 2], 2)
        tmp$range <- round(tmp2[, 3], 2)
        tmp$singular <- c(test.Exp, test.Sph, test.Gau)
        Results <- rbind(Results, tmp)
        Fitted <- data.frame(dist = seq(0.01, max(Vario1$dist), length=300))
        Fitted$gamma.Exp <- gstat::variogramLine(Vario1.fit.Exp,
                                                 dist_vector =
                                                     Fitted$dist)$gamma
        Fitted$gamma.Sph <- gstat::variogramLine(Vario1.fit.Sph,
                                                 dist_vector =
                                                     Fitted$dist)$gamma
        Fitted$gamma.Gau <- gstat::variogramLine(Vario1.fit.Gau,
                                                 dist_vector =
                                                     Fitted$dist)$gamma
        Empirical <- reshape::melt(Vario1,
                                   id.vars = "dist",
                                   measure.vars = c("gamma"))
        Modeled <- reshape::melt(Fitted,
                                 id.vars = "dist",
                                 measure.vars = c("gamma.Exp",
                                                  "gamma.Sph",
                                                  "gamma.Gau"))
        pltName <- paste0("VarioPlot.", round(i / 1000, 2))
        pltList[[pltName]] <<- ggplot2::ggplot(Empirical,
                                               ggplot2::aes(x = dist,
                                                            y = value,
                                                            colour = variable))
        +
            ggplot2::geom_point() +
            ggplot2::geom_line(data = Modeled) +
            ggplot2::ylab("Semivariance") +
            ggplot2::xlab("Distance (m)") +
            ggplot2::ggtitle(paste(ifelse(type == "TRUE",
                                          "Cressie Semivariance \n",
                                          "Semivariance \n"),
                                   "Maximum Distance=",
                                   round(i / 1000, 2),
                                   "km")) +
            ggplot2::theme(legend.title = ggplot2::element_blank()) +
            ggplot2::scale_colour_discrete(breaks = c("gamma",
                                                      "gamma.Exp",
                                                      "gamma.Sph",
                                                      "gamma.Gau"),
                                           labels = c("Observed",
                                                      "Exponential",
                                                      "Spherical",
                                                      "Gaussian")) +
            ggplot2::stat_smooth()
    }
    Results <- Results[-1, c(1, 6, 2:5, 7:9)]
    Moran.test <- Moran.I.calc(data)
    Results$observed <- round(rep(Moran.test$observed, length(Results[, 1])), 4)
    Results$expected <- round(rep(Moran.test$expected, length(Results[, 1])), 4)
    Results$sd <- round(rep(Moran.test$sd, length(Results[, 1])), 4)
    Results$p.value <- round(rep(Moran.test$p.value, length(Results[, 1])), 4)
    Results <<- Results
    return(Results)
}
# ------------------------------------------------------------------------------
