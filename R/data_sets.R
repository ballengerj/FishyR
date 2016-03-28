#' Age Determinations of Red Porgy
#'
#' A dataset containing multiple age determinations for Red Porgy made by four
#'    different readers.  Each reader did not age each individual fish.
#'
#' @format A data frame with 1511 rows and 4 variables:
#' \describe{
#'  \item{R1.Age}{biological age, in years, determiend by reader 1}
#'  \item{R2.Age}{biological age, in years, determiend by reader 2}
#'  \item{R3.Age}{biological age, in years, determiend by reader 3}
#'  \item{R4.Age}{biological age, in years, determiend by reader 4}
#' }
#' @source South Carolina Deparment of Natural Resources, Coastal Finfish Section
#'  (ballengerj@dnr.sc.gov)
"RedPorgy"

#' Black Sea Bass Life History Data
#'
#' A dataset containing black sea bass life history data, as collected by the
#'    South Carolina Department of Natural Resources, Marine Resources Research
#'    Institute, Coastal Finfish Section, Reef Fish Survey.
#'
#' @format A data frame with 33,297 rows and 13 variables
#' \describe{
#'  \item{PCGSS}{unique identifier identifying individual fish}
#'  \item{Year}{year of capture}
#'  \item{Month}{month of capture}
#'  \item{Latitude}{latitude, in decimal degrees, where captured}
#'  \item{Depth}{depth, in meters, where captured}
#'  \item{TL}{total length, in mm}
#'  \item{SL}{standard length, in mm}
#'  \item{WholeWt}{whole weight, in grams}
#'  \item{GonadWt}{gonad weight, in grams}
#'  \item{Sex}{sex}
#'  \item{Mat}{maturity stage}
#'  \item{Age}{biological age}
#'  \item{Year.Class}{year born}
#' }
#' @source South Carolina Deparment of Natural Resources, Coastal Finfish Section
#'  (ballengerj@dnr.sc.gov)
"BSB"

#' Sheepshead Marginal Increment Data
#'
#' A dataset containing Sheepshead marginal increment measurements.
#'
#' @format A data frame with 416 rows and 4 variables
#' \describe{
#'  \item{Month}{numeric, month of capture of individual fish}
#'  \item{MI}{numeric, marginal increment measurement, in micrometers}
#'  \item{Prev.Inc}{numeric, width of previous increment, in micrometers}
#'  \item{Oto.Age}{numeric, age, in years, of individual fish}
#' }
#' @source Joseph C. Ballenger (\email{ballengerj@dnr.sc.gov})
"Sheepshead"
