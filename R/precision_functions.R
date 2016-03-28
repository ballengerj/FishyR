#' Average Percent Error (APE) Estimation
#'
#' \code{APE} calculates the average percent error (APE) for multiple age
#'    determinations on the same fish.  This function use the mean of the
#'    observed ages as the basis for APE calculation.
#' @param x Numeric vector of independent age determinations from a single
#'    fish.
#' @return The APE of the vector of ages \code{x} from an individual fish.
#' @examples
#' # Example giving an estimate with an APE of 0.0
#' APE(c(1, 1, 1))
#'
#' # Example using missing data
#' APE(c(1, NA, 2))
#'
#' # Another example with missing data
#' APE(c(3, 5, 3, 4, NA, 6))
#'
#' #Example showing how APE.Med can be used to calculate mean APE for many fish
#' w <- seq(1, 10, 1)
#' x <- seq(2, 11, 1)
#' y <- seq(3, 12, 1)
#' z <- seq(4, 13, 1)
#' df <- data.table::rbindlist(apply(data.frame(w, x, y, z), 1, APE)); df
#' mean(data.frame(df)$APE)
#' # Example with RedPorgy data
#' data(RedPorgy)
#' df <- data.table::rbindlist(apply(na.exclude(RedPorgy[, c(1, 2)]), 1,
#'  APE))
#' df
#' mean(data.frame(df)$APE)
#' @family APE functions
#' @export
APE <- function(x) {
  r <- length(na.exclude(x)) # length of x
  avg <- mean(na.exclude(x)) # mean of x
  Diff <- abs(na.exclude(x) - avg) # absolute value of the difference between x
  # and avg
  PE <- Diff / (r * avg) # perecent error for each x value
  APE <- sum(PE) # get the total percent error for all x values
  data.frame(r   = r, # combine the length of x, the mean of x, and the average
             avg = avg, # percent error of x into a data frame for output
             APE = APE)
}

# ------------------------------------------------------------------------------
#' Average Percent Error (APE) Estimation
#'
#' \code{APE.Med} calculates the average percent error (APE) for multiple age
#'    determinations on the same fish.  This function use the median of the
#'    observed ages as the basis for APE calculation.
#' @param x Numeric vector of independent age determinations from a single
#'    fish.
#' @return The APE of the vector of ages \code{x} from an individual fish.
#' @examples
#' #Example giving an estimate with an APE.Med of 0.0
#' APE.Med(c(1, 1, 1))
#'
#' #Example using missing data
#' APE.Med(c(1, NA, 2))
#'
#' #Another example with missing data
#' APE.Med(c(3, 5, 3, 4, NA, 6))
#'
#' #Example showing how APE.Med can be used to calculate mean APE for many fish
#' w <- seq(1, 10, 1)
#' x <- seq(2, 11, 1)
#' y <- seq(3, 12, 1)
#' z <- seq(4, 13, 1)
#' df <- data.table::rbindlist(apply(na.exclude(RedPorgy[, c(1, 2)]), 1,
#'  APE.Med))
#' df
#' mean(data.frame(df)$APE)
#' # Example with RedPorgy data
#' data(RedPorgy)
#' df <- data.table::rbindlist(apply(na.exclude(RedPorgy[, c(1,2)]), 1,
#'  APE.Med))
#' df
#' mean(data.frame(df)$APE)
#' @family APE functions
#' @seealso \code{\link{APE}} for APE calculation using the mean age as
#'    the basis APE calculation.
#' @export
APE.Med <- function(x) {
  r <- length(na.exclude(x)) # length of x
  med <- median(na.exclude(x)) # median of x
  B.j <- sum(abs(na.exclude(x) - med) / r) # calculate the difference between x
  # and median for each x, divide this value by the length of x, and sum those
  # for every x
  APE <- B.j / med # calculate the averge percent error for x
  data.frame(r      = r, # combine the length of x, the median value of x, and
             median = med, # the averge percent error of x into a data frame
             APE    = APE) # for output
}

# ------------------------------------------------------------------------------
#' Percent Agreement Estimation
#'
#' \code{perc.agr} calculates the percent agreement of two sets of age reads on
#'    the same fish.  This function will calculate the percent agreement from
#'    +/- 0 years to +/- the number of years defined by \code{levels},
#'    by 1 year increments.
#' @param Data Numeric data frame, with no \code{\link{NA}} values,
#'    representing paired age reads on the same individual fish.
#' @param levels Numeric integer representing the maximum number of years you
#'    want to consider for the percent agreement analysis
#' @return Object of class \code{\link{data.table}} and
#'    \code{\link{data.frame}} containing the percent agreement metrics for
#'    each year from +/- years to \code{levels}.  The \code{\link{data.table}}
#'    will contain four columns:
#'    \describe{
#'      \item{Scale}{+/- number of years considered for percent agreement
#'          calculation}
#'      \item{n}{total sample size for the analysis}
#'      \item{agree}{# of samples agreed upon by the readers, given that scale}
#'      \item{perc.agree}{proportinal agreement among the individual age reads}
#'      }
#' @examples
#' #Example - RedPorgy Dataset
#' data(RedPorgy)
#' perc.agr(na.exclude(RedPorgy[, c(1, 2)]), levels = 3)
#' @export
perc.agr <- function(Data, levels= 2) {
  # Setup analysis list --------------------------------------
  list.names <- seq(0, levels, 1) # create names for list
  list <- vector("list", length(list.names)) # create an empty list
  names(list) <- list.names # names the elements of the empty list
  # Finding total sample size --------------------------------
  total.n <- length(Data[, 1])
  # Calculating percent agreement ----------------------------
  # Calculates percent agreement from +/- 0 years to +/- (levels) years by 1
  # year increments
  for(i in seq(0, levels, 1)) {
    Diff <- abs(Data[, 1] - Data[, 2]) # calculate absolute value of difference
    Diff <- subset(Diff, Diff <= i) # between reads and only retains those
      # difference values if they are > i
    agree.n <- length(Diff) # determine numbers of times "agreed"
    PA <- length(Diff) / total.n # cacluate percent agreement
    list[[as.character(i)]] <- data.frame(Scale = i, n = total.n, # saving
                                          agree = agree.n, perc.agree = PA)
    # results to the list created earlier
  }
  data.table::rbindlist(list) # converting the list of data frames to a single
  # data frame for output
}

# ------------------------------------------------------------------------------
#' Bowker's Symmetry Test
#'
#' \code{symmetry} Performs Bowker's symmetry test comparing two sets of age
#'    reads on the same fish.
#' @param data Numeric data frame, with no \code{\link{NA}} values,
#'    representing paired age reads on the same individual fish.
#' @return Object of class \code{\link{data.frame}} containing the results of
#'    the Bowker's symmetry test.  The \code{\link{data.frame}} will contain
#'    three columns:
#'    \describe{
#'      \item{chisq}{Chi-square statistic}
#'      \item{df}{estimated degrees of freedom for the chi-square test}
#'      \item{p.value}{p-value of the symmetry test}
#'      }
#' @examples
#' #Example - RedPorgy Dataset
#' data(RedPorgy)
#' symmetry(na.exclude(RedPorgy[, c(1, 2)]))
#' @export
symmetry <- function(data) {
  table <- table(data[, 1], data[, 2])
  # Creating a vector of ages from minimum to maximum age recorded for either
  # reader
  Ages <- seq( min( data, na.rm = TRUE), max(data, na.rm = TRUE) )
  x <- as.numeric( unlist( dimnames(table)[1] ) )
  y <- as.numeric( unlist( dimnames(table)[2] ) )
  # Creating an empty matrix with a row and column for each and every age found
  # in the above vector
  temp.matrix <- test.temp <- matrix( rep(NA, length(Ages) ^ 2),
                                      nrow = length(Ages),
                                      ncol = length(Ages),
                                      dimnames = list(Ages, Ages) )
  # Placing observed frequency counts into temporary matrix
  for (i in x) {
    k <- as.character(i)
    for (j in y) {
      l <- as.character(j)
      temp.matrix[k, l] <- table[k, l]
    }
  }

  # Replacing NA values in temp.matrix with 0
  temp.matrix <- ifelse( is.na(temp.matrix == TRUE), 0, temp.matrix)

  # Fill the dummy contingency table with symmetry calcualtions and set all
  # diagonal elements equal to 0
  for (i in 1:dim(temp.matrix)[1]) {
    for (j in 1:dim(temp.matrix)[2]) {
      test.temp[i, j] <- ( (temp.matrix[i, j] - temp.matrix[j, i]) ^ 2) /
        (temp.matrix[i, j] + temp.matrix[j, i])
    }
  }
  diag(test.temp) <- NA

  # Performing the Chi-squared symmetry test
  symmetry.test <- sum(test.temp, na.rm = TRUE) # Calculating observed chi-
      # square value
  test2.temp <- as.vector(test.temp)
  symmetry.df <- length( na.exclude(test2.temp) ) # Calculating the degrees of
      # freedom of the chi-square test
  p.value <- 1 - pchisq(symmetry.test, symmetry.df) # Calculating the p-value
      # of the chi-square test

  # Outputting Results
  data.frame(chisq = symmetry.test, df = symmetry.df, p.value = p.value)
}

# ------------------------------------------------------------------------------
#' Distribution Plot
#'
#' \code{dist.plot} creates a distribution plot that helps the researcher
#'    assess whether bias is present between two sets of age reads on the same
#'    individual fish.  The plot is constructed using lattice graphics.  It
#'    provides three unique pieces of information:
#'    \enumerate{
#'        \item 1:1 line showing the expected distribution of reads given no
#'            bias in age determinations between the paired reads
#'        \item linear regression line fit to the observed paired data for the
#'            individual reads
#'        \item text providing the number of times reader 2 assigned a fish
#'            an age, given reader 1 assigned that fish age "x"
#'    }
#' @param R1 Numeric vector of independent age determinations from a single
#'      fish from Reader 1.  This must be the same length as R2.
#' @param R2 Numeric vector of independent age determinations from a single
#'    fish from Reader 2.  This must be the same length as R1.
#' @param R1.label Text label for the x-axis (defaults to Reader 1) for the
#'    resultant plot
#' @param R2.label Text label for the y-axis (defaults to Reader 2) for the
#'    resultant plot
#' @return A \code{\link{xyplot}} using the package lattice depicting the
#'    distribution of reader 2 ages, given reader 1 assigned a fish a certain
#'    age.
#' @examples
#' # Example with Red Porgy data
#' data(RedPorgy)
#' dist.plot(R1 = RedPorgy$R3.Age, R2 = RedPorgy$R4.Age, R1.label = "Reader 3",
#'    R2.label = "Reader 4")
#' @family Bias Plot Functions
#' @export
dist.plot <- function(R1, R2, R1.label = "Reader 1", R2.label = "Reader 2") {
  data.table <- data.frame( table(R1, R2) )
  df <- data.frame(R1 = as.integer( as.character( data.table[, 1]) )
                   [data.table[, 3] != 0],
                   R2 = as.integer( as.character( data.table[, 2]) )
                   [data.table[, 3] != 0],
                   Freq = paste0( data.table[, 3][data.table[, 3] != 0]) )
  m1 <- lm(R2~R1)
  lattice::xyplot(R2 ~ R1,
                  data = df,
                  type = c("n", "g"),
                  xlab = R1.label,
                  ylab = R2.label,
                  xlim = c( min( df[, c(1 : 2)]) - 0.5,
                            max( df[, c(1 : 2)]) + 0.5 ),
                  ylim = c( min( df[, c(1 : 2)]) - 0.5,
                            max( df[, c(1 : 2)]) + 0.5 ),
                  panel = function(x, y) {
                    lattice::panel.xyplot(x, y, type="n")
                    lattice::panel.abline(a = 0, b = 1, col = "black")
                    lattice::panel.abline(a = m1$coef[1], b = m1$coef[2],
                                          col = "black", lty = 2)
                    lattice::panel.text(x = x, y = y, labels = df$Freq)
                  } )
}

# ------------------------------------------------------------------------------
#' Bias Plot
#'
#' \code{bias.plot} creates a bias plot (see Campana 2001) that helps the
#'    researcher assess whether bias is present between two sets of age reads
#'    on the same individual fish.  The plot is constructed using lattice
#'    graphics.
#' @param Data Numeric data frame having columns titled the following
#'    \describe{
#'        \item{Age}{numeric vector representing the unique ages assigned to
#'            indvidual fish by Reader 1}
#'        \item{Avg}{numeric vector representing the average age given by
#'            Reader 2 for fish aged "x' by Reader 1}
#'        \item{CI}{numeric vector representing the confidence interval of ages
#'            given by Reader 2 for fish aged "x" by Reader 1}
#'        }
#'    Other columns may be present in the data frame.
#' @param R1.label Text label for the x-axis (defaults to Reader 1) for the
#'    resultant plot
#' @param R2.label Text label for the y-axis (defaults to Reader 2) for the
#'    resultant plot
#' @param a numeric; y-axis intercept for comparison line
#' @param b numeric; slope estimate for comparison line
#' @param pch see \code{\link{pch}}; symbol for points used in plot
#' @param col.points A specification for the plotting color of the points.  See
#'    description of \emph{Color Specificaiton} in the \code{\link{par}} for
#'    more information.
#' @param lty.ab specify line type for the "ab" line.  See description of
#'    \emph{lty} in the \code{\link{par}} for more information
#' @param lty.col specify the line color for the "ab" line.  See
#'    description of \emph{Color Specificaiton} in the \code{\link{par}} for
#'    more information.
#' @param err.col specify the line color for the error bars.  See
#'    description of \emph{Color Specificaiton} in the \code{\link{par}} for
#'    more information.
#' @return A \code{\link{xyplot}} using the package lattice depicting the
#'    bias plot between Reader 2 and Reader 1 ages.
#' @examples
#' # Example with RedPorgy data
#' data(RedPorgy)
#' df <- with(RedPorgy, aggregate(R4.Age ~ R1.Age, FUN = describe, digits = 4 ))
#' df <- cbind(df[-ncol(df)], df[[ncol(df)]])
#' colnames(df) <- c("Age" ,"n", "Avg", "StDev", "StErr", "CI", "CV","RSE",
#'    "RMSD", "MAE")
#' bias.plot(Data = df)
#' @family Bias Plot Functions
#' @export
bias.plot <- function(Data, R1.label = "Reader 1", R2.label = "Reader 2",
                      a = 0, b = 1, pch = 16, col.points = "black", lty.ab = 2,
                      lty.col = "black", err.col = "black") {
  grid <- lattice::xyplot( (Avg + CI) + (Avg - CI) ~ (Age + CI) + (Age - CI),
                          data = Data,
                          xlab = R1.label,
                          ylab = R2.label,
                          type = "g")
  points <- lattice::xyplot(Avg ~ Age,
                            data = Data,
                            panel = function(x, y) {
                              lattice::panel.xyplot(x, y, type = "p",
                                                    pch = pch, col = col.points)
                              lattice::panel.abline(a = a, b = b, lty = lty.ab,
                                                    col = lty.col)
                              errbar(x = Data$Age, y = Data$Avg,
                                     error = Data$CI, color = err.col)
                   })
  grid + latticeExtra::as.layer(points)
}
