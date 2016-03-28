#' Rounding
#'
#' \code{roundup} is a function that allows for traditional rounding of numbers.
#'  The built in R rounding function (\code{\link{round}}) from the R base
#'  package follows the standard of the IEC 60559 when determining how to round
#'  off the number 5. This rule states that you '\emph{go to the even digit}'.
#'  This is different from how most people learn to round off 5 in school (and
#'  different from how Excel rounds numbers), with one rounding 5 off to the
#'  next highest absolute value integer. To illustrate the difference, where
#'  the base function rounds 0.5 and -2.5  to 0 and -2, respectively, roundup
#'  rounds 0.5 and -2 to  1 and -3, respectively.
#' @param x a numeric vector.
#' @param numdigits positive integer indicating the number of decimal places to
#'  be used.
#' @return numeric vector of \code{x} rounded to the number of digits specified
#'  by \code{numdigits}.
#' @examples
#' # Examples
#' roundup(-2.25, numdigits=1)
#' round(.5 + -2:4) # IEEE rounding
#' roundup(.5+ -2:4) # Traditional rounding
#' # Comparing IEEE Rounding vs. Traditional Rounding
#' data(BSB)
#' IEEE <- with(BSB, round(TL / 10))
#' Trad <- with(BSB, roundup(TL / 10))
#' table(IEEE,Trad)
#' plot(Trad ~ jitter(IEEE), xlab = "IEEE Total Length (cm)",
#'  ylab = "Traditional Total Length (cm)")
#' @family Rounding functions
#' @seealso \code{\link{round}} for base rounding function.
#' @export
roundup <- function(x, numdigits = 0) {
    x <- x * 10 ^ numdigits # Multiplying original values by some power of 10,
      # as defined using "numdigits
    x <- ifelse(x < 0, -trunc(abs(x) + 0.5), trunc(x + 0.5))
    x / 10 ^ numdigits
}

# ------------------------------------------------------------------------------
#' Rounding to a Power of 10
#'
#' \code{roundpow} is a function that allows for traditional rounding of
#'  numbers to a power of 10.  The built in R rounding function
#'  (\code{\link{round}}) from the R base package follows the standard of the
#'  IEC 60559 when determining how to round off the number 5. This rule states
#'  that you '\emph{go to the even digit}'. This is different from how most
#'  people learn to round off 5 in school (and different from how Excel rounds
#'  numbers), with one rounding 5 off to the next highest absolute value
#'  integer. The function \code{\link{roundup}} in this package allows one to
#'  apply the traditinal rounding measure.  This function extends
#'  \code{\link{roundup}} by allowing one to round numbers to a power of 10.
#' @param x a numeric vector.
#' @param power positive integer indicating the power of ten to round too.
#' @return numeric vector of \code{x} rounded to the power of 10 specified by
#' \code{power}.
#' @examples
#' # Examples
#' roundpow(125, power=1)
#' roundpow(125, power=2)
#' roundpow(10555, power=1)
#' roundpow(10555, power=2)
#' roundpow(10555, power=3)
#' roundpow(10555, power=4)
#' # Comparing IEEE Rounding vs. Traditional Rounding
#' data(BSB)
#' IEEE <- with(BSB, round(TL / 10) * 10)
#' Trad <- with(BSB, roundpow(TL, 1))
#' plot(Trad ~ jitter(IEEE), xlab = "IEEE Total Length (mm)",
#'  ylab = "Traditional Total Length (mm)")
#' @family Rounding functions
#' @export
roundpow <- function(x, power = 1) {
    (roundup(roundup(x) / 10 ^ power)) * 10 ^ power
}
