#' Red Drum catches from the South Carolina Department of Natural Resources
#' Trammel Net Survey
#'
#' @format A tibble with 20,366 rows and 15 variables:
#' \describe{
#'   \item{Collection}{Factor Collection identifier}
#'   \item{Catch}{dbl Catch of red drum}
#'   \item{Random}{Factor Denoting whether site was randomly (1) or
#'   non-randomly (0) selected for sampling}
#'   \item{Year}{Factor Year of collection}
#'   \item{Month}{Factor Month of collection}
#'   \item{Day}{dbl Day of collection}
#'   \item{Date}{POSICct Date of collection}
#'   \item{Site}{Factor Site of collection; \code{site} is nested within
#'   \code{stratum} and \code{area}}
#'   \item{Stratum}{Factor Stratum of collection with seven possible values:
#'   Port Royal Sound, ACE Basin, Ashley River, Charleston Harbor, Wando River,
#'   Cape Romain, and Winyah Bay. Stratum is nested within \code{Area}}
#'   \item{Area}{Factor Area (Estuary) of collection with five possible
#'   values: Port Royal Sound, ACE Basin, Charleston Harbor, Cape Romain and
#'   Winyah Bay}
#'   \item{Temp_Water}{dbl Water temperature (degrees Celsius) at site of
#'   collection}
#'   \item{Salinity}{dbl Salinity (PSU) at site of collection}
#'   \item{DO}{dbl Dissolved oxygen (in mg per liter) in the water at site of
#'   collection}
#'   \item{Tide}{Factor Tidal stage at time of collection}
#'   \item{Day_of_Year}{dbl Day of year of collection}
#' }
"RedDrum"
