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

#' Bluenose (Hyperoglyphe antarctica) Bomb Radiocarbon Data
#'
#' A dataset containing Bluenose Bomb Radiocarbon data, as provided in Francis
#'   et al. 2010.
#'
#' @format A data frame with 12 rows and 9 variables
#' \describe{
#'  \item{fish}{factor, variable identifying individual fish}
#'  \item{age}{integer, age estimate}
#'  \item{catch.date}{factor, catch date}
#'  \item{C14year}{integer, estimated eyar of birth based on age and capture
#'     year}
#'  \item{C14}{numeric, C-C14 estimate}
#'  \item{C14.err}{numeric, C-14 estimate error estimate}
#'  \item{age.err}{numeric, age estimate error estimate}
#'  \item{samp.age}{integer, age estimate of analyzed sample}
#'  \item{catch.year}{capture year of fish}}
#'  @source Francis, R.I.C., Campana, S.E., and H.L. Neil.  2010.  Validation
#'     of fish ageing methods should involve bias estimation rather than
#'     hypothesis testing: a proposed approach for bomb radiocarbon
#'     validations.  Can. J. Fish. Aquat. Sci. 67: 1398-1408.
"bluenose"

#' Snapper (Pagrus auratus) Bomb Radiocarbon Data
#'
#' A dataset containing Snapper Bomb Radiocarbon data, as provided in Kalish
#'   1993.
#'
#' @format A data frame with 16 rows and 7 variables
#' \describe{
#'  \item{age}{integer, age estimate}
#'  \item{catch.year}{capture year of fish}
#'  \item{C14year}{integer, estimated birth year based on age and capture year}
#'  \item{C14}{numeric, C-14 estimate}
#'  \item{C14.err}{numeric, C-14 estimate error estimate}
#'  \item{age.err}{numeric, age estimate error estimate}
#'  \item{samp.age}{integer, age estimate of analyzed sample}}
#' @source Kalish, J.M.  1993.  Pre- and post-bomb radiocarbon in fish
#'    otoliths.  Earth Planet. Sci. Lett. 114(4): 549-554.
"snapper"

#' Adult Haddock (Melanogrammus aeglefinus) Bomb Radiocarbon Data
#'
#' A dataset containing adult Haddock Bomb Radiocarbon data, as provided in
#'   Campana 1997 (test data in table 1).
#'
#' @format A data frame with 11 rows and 7 variables
#' \describe{
#'  \item{age}{integer, age estimate}
#'  \item{catch.year}{capture year of fish}
#'  \item{C14year}{integer, estimated birth year based on age and capture year}
#'  \item{C14}{numeric, C-14 estimate}
#'  \item{C14.err}{numeric, C-14 estimate error estimate}
#'  \item{age.err}{numeric, age estimate error estimate}
#'  \item{samp.age}{integer, age estimate of analyzed sample}}
#' @source Campana, S.E.  1997.  Use of radiocarbon from nuclear fallout as a
#'    dated marker in the otoliths of haddock Melanogrammus aeglefinus.  Mar.
#'    Ecol. Prog. Ser. 150: 49-56.
"haddock"

#' Northwest Atlantic Bomb Radiocarbon Reference Data
#'
#' A dataset containing bomb radiocarbon data derived from juvenile Haddock,
#'    redfish (Sebastes spp.) and yellowtail flounder collected in the
#'    Northwest Atlntic, adult Haddock Bomb Radiocarbon data, as provided in
#'    Campana et al. 2008 (Table 2).
#'
#' @format A data frame with 11 rows and 8 variables
#' \describe{
#'  \item{age}{integer, age estimate}
#'  \item{catch.year}{capture year of fish}
#'  \item{C14year}{integer, estimated birth year based on age and capture year}
#'  \item{C14}{numeric, C-14 estimate}
#'  \item{C14.err}{numeric, C-14 estimate error estimate}
#'  \item{age.err}{numeric, age estimate error estimate}
#'  \item{samp.age}{integer, age estimate of analyzed sample}
#'  \item{YC}{integer, capture year}}
#' @source Campana, S.E., Casselman, J.M., and Jones, C.M.  2008.  Bomb
#'    radiocarbon chronologies in the Arctic, with implications for the age
#'    validation of lake trout (Salvelinus namaycush) and other Arctic
#'    species.  Can. J. Fish. Aquat. Sci. 65(4): 733-743.
"nwa"
