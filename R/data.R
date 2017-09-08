#' Australian Institute of Sport data
#'
#' Data on 102 male and 100 female athletes collected at the Australian Institute of Sport, courtesy of Richard Telford and Ross Cunningham.
#' @format A data frame with 202 observations on the following 13 variables:\cr
#' \itemize{
#' \item{[,1] \code{sex}}{ - categorical, levels = \code{female, male}}
#' \item{[,2] \code{sport}}{ - categorical, levels = \code{B_Ball, Field, Gym, Netball, Row, Swim, T_400m, Tennis, T_Sprnt, W_Polo}}
#' \item{[,3] \code{RCC}}{ - red cell count (numeric)}
#' \item{[,4] \code{WCC}}{ - white cell count (numeric)}
#' \item{[,5] \code{Hc}}{ - Hematocrit (numeric)}
#' \item{[,6] \code{Hg}}{ - Hemoglobin (numeric)}
#' \item{[,7] \code{Fe}}{ - plasma ferritin concentration (numeric)}
#' \item{[,8] \code{BMI}}{ - body mass index: \code{Wt/(Ht)^2} (numeric)}
#' \item{[,9] \code{SSF}}{ - sum of skin folds (numeric)}
#' \item{[,10] \code{Bfat}}{ - body fat percentage (numeric)}
#' \item{[,11] \code{LBM}}{ - lean body mass (numeric)}
#' \item{[,12] \code{Ht}}{ - height, cm (numeric)}
#' \item{[,13] \code{Wt}}{ - weight, kg (numeric)}
#' }
#' @details The data have been made publicly available in connection with the book by Cook and Weisberg (1994).
#' @references Cook, R. D. and Weisberg, S. (1994), \emph{An Introduction to Regression Graphics}. John Wiley & Sons, New York.
#' @examples
#' data(ais, package="MoEClust")
#' pairs(ais[,c(3:4, 10:13)], col=as.numeric(ais[,1]), main = "AIS data")
#' @docType data
#' @usage data(ais)
"ais"

#' GNP and CO2 Data Set
#'
#' This data set gives the gross national product (GNP) per capita in 1996 for various countries as well as their estimated carbon dioxide (CO2) emission per capita for the same year.
#'
#' @format This data frame consists of 28 countries and the following variables:\cr
#' \itemize{
#' \item{\code{GNP}}{ - The gross product per capita in 1996.}
#' \item{\code{CO2}}{ - The estimated carbon dioxide emission per capita in 1996.}
#' \item{\code{country}}{ - An abbreviation pertaining to the country measures (e.g. \code{"GRC"} = Greece and \code{"CH"} = Switzerland).}
#' }
#' @references Hurn, M., Justel, A. and Robert, C. P. (2003) Estimating Mixtures of Regressions, \emph{Journal of Computational and Graphical Statistics}, 12(1): 55-79.
#' @docType data
#' @usage data(CO2data)
"CO2data"
