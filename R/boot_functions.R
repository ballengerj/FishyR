#' Calculate bootstrap diagnostics from relative abundance index results
#'
#' \code{Boot.Diag} generates a data frame containing summary statistics for
#'   various number of bootstrap runs.  Intent is to inform when uncertainty
#'   estimates of a relative abundance index stabilize.
#' @param Rows TBD
#' @param Columns TBD
#' @param By TBD
#' @param Statistics TBD
#' @return Data frame with five columns:
#'   \describe{
#'     \item{Group}{}
#'     \item{Points}{}
#'     \item{Factor}{}
#'     \item{Frame}{}
#'   }
#' @family Index Bootstrap Functions
#' @export
Boot.Diag <- function(Rows, Columns, By, Statistics=c("CV","Var")) {
    Group <- rep(sort(rep(Columns, length(seq(1, Rows, by=By)))),
                 length(Statistics))
    Points <- rep(seq(By, Rows, by = By), length(Columns) * length(Statistics))
    Factor <- sort(rep(Statistics, length(Columns) * length(seq(1, Rows,
                                                                by = By))))
    Frame <- data.frame(Points = Points, Group = Group, Factor = Factor,
                        stringsAsFactors = TRUE)
    return(Frame)
}
#-------------------------------------------------------------------------------
#' Calculate bootstrap statistics
#'
#' \code{Stat.Calc} calculates statistics (i.e., variance and cv estimates)
#'   about a relative abundance index from a bootstrap analysis.
#' @param X TBD
#' @param data TBD
#' @param index TBD
#' @return numeric vector with four points:
#'   \describe{
#'     \item{Points}{}
#'     \item{Group}{}
#'     \item{Factor}{}
#'     \item{Value}{}
#'   }
#' @family Index Bootstrap Functions
#' @export
Stat.Calc <- function(X, data, index) {
    if(X[, 3] == "Var") {
        Value <- var(data[1 : X[, 1], X[, 2]])
        vector <- c(Points = X[, 1], Group = X[, 2], Factor = X[, 3],
                    Value = Value)
    }
    else {
        Value <- sd(data[1 : X[, 1], X[, 2]]) / index[X[, 2]]
        vector <- c(Points = X[, 1], Group = X[, 2], Factor = X[, 3],
                    Value = Value)
    }
    return(vector)
}
#-------------------------------------------------------------------------------
#' Zero-Inflated Model Bootstrap Function
#'
#' \code{ZI.Index.Boot} calculates a raw, a normalized, a centered, and a scaled
#'   relative abundance index for each bootstrap iteration
#' @param data see \code{\link[boot]{boot}}
#' @param i see \code{\link[boot]{boot}} for \code{stype}
#' @param bestmod character vector representing saved object representing the best
#'   (i.e., final) model
#' @param dist character vector used to specify the distribution to use in the
#'   \code{\link[pscl]{zeroinfl}} function being implemented
#' @param link character vector used to specify the link to use in the
#'   \code{\link[pscl]{zeroinfl}} function being implemented
#' @param grid prediction data frame for which you will get your predicted
#'   response variable at.  Generally will be created via a call to
#'   \code{\link[FishyR]{Index.Grid}}
#' @return numeric vector with four points:
#'   \describe{
#'     \item{Points}{}
#'     \item{Group}{}
#'     \item{Factor}{}
#'     \item{Value}{}
#'   }
#' @seealso \code{\link[boot]{boot}}
#' @family Index Bootstrap Functions
#' @export
ZI.Index.Boot <- function(data, i, bestmod, dist = "negbin", link = "logit",
                          grid) {
    fit <- try(pscl::zeroinfl(formula(get(bestmod)), data = data[i, ],
                              dist = dist, link = link))
    Raw <- try(aggregate(predict(object = fit, newdata = grid) ~ grid[, 1],
                         FUN = mean)[, 2])
    Normalized <- try(Raw / mean(Raw))
    Centered <- try(scale(Raw, scale = FALSE))
    Scaled <- try(scale(Raw))
    c(Raw, Normalized, Centered, Scaled)
}
