#' Example sim123 Dataset
#'
#' A simple dataset for demonstration purposes.
#'
#' @format A list containing 3 data frames:
#' \describe{
#'   \item{smat}{A 1089 x 3 numeric matrix of the true source signals. Each column is an 33 x 33 image.}
#'   \item{mmat}{A 3 x 50 numeric mixing matrix of the true time series. Each row is a time series of corresponding column in \code{smat}.}
#'   \item{xmat}{A 1089 x 50 numeric matrix of the simulated data. Each column is the simulated mixed signal at a time point.}
#' }
#' @examples
#' data(example_sim123)
#' str(example_sim123)
"example_sim123"
