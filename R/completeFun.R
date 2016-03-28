#' Extract Complete Cases from a Data Frame
#'
#' \code{complete.data} allows one to extract complete cases (i.e., no NA values),
#' based on a potential subset of the total columns in the data frame, from a
#' data frame.
#' @param data Data frame
#' @param Cols vector, which may represent a character vector or numeric vector,
#' of columns for which we want to determine the complete cases
#' (\code{\link{complete.cases}}) for
#' @return A restricted data frame that contains all complete cases for the columns specified.
#' @examples
#' (df <- data.frame(w = c(1, 1, NA, 2, 3, NA), x = c(2, 6, NA, 7, 8, 10),
#' y = c(2, NA, 10, 5, 6, 1), z = c(NA, 1, 3, NA, 4, 5))) # Original Data
#' complete.data(df, Cols = c("w", "x")) # Using named columns
#' complete.data(df, Cols = c(1, 2))
#' @export
complete.data <- function(data, Cols) {
  complete.vector <- complete.cases(data[, Cols])
  return(data[complete.vector, ])
}
