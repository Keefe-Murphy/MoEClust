#' Australian Institute of Sport data
#'
#' Data on 102 male and 100 female athletes collected at the Australian Institute of Sport, courtesy of Richard Telford and Ross Cunningham.
#' @format A data frame with 202 observations on the following 13 variables:\cr
#' \describe{
#' \item{\code{sex}}{categorical, levels = \code{female, male}}
#' \item{\code{sport}}{categorical, levels = \code{B_Ball, Field, Gym, Netball, Row, Swim, T_400m, Tennis, T_Sprnt, W_Polo}}
#' \item{\code{RCC}}{red cell count (numeric)}
#' \item{\code{WCC}}{white cell count (numeric)}
#' \item{\code{Hc}}{Hematocrit (numeric)}
#' \item{\code{Hg}}{Hemoglobin (numeric)}
#' \item{\code{Fe}}{plasma ferritin concentration (numeric)}
#' \item{\code{BMI}}{body mass index: \code{Wt/(Ht)^2} (numeric)}
#' \item{\code{SSF}}{sum of skin folds (numeric)}
#' \item{\code{Bfat}}{body fat percentage (numeric)}
#' \item{\code{LBM}}{lean body mass (numeric)}
#' \item{\code{Ht}}{height, cm (numeric)}
#' \item{\code{Wt}}{weight, kg (numeric)}
#' }
#' @details The data have been made publicly available in connection with the book by Cook and Weisberg (1994).
#' @references Cook, R. D. and Weisberg, S. (1994), \emph{An Introduction to Regression Graphics}. Volume 405 of \emph{Wiley Series in Probability and Statistics}, New York, NY, USA: John Wiley & Sons.
#' @examples
#' data(ais, package="MoEClust")
#' pairs(ais[,c(3:7)], col=as.numeric(ais$sex), main = "AIS data")
#' apply(ais[,c(3:7)], 2, summary)
#' @docType data
#' @keywords datasets
#' @usage data(ais)
"ais"

#' GNP and CO2 Data Set
#'
#' This data set gives the gross national product (GNP) per capita in 1996 for various countries as well as their estimated carbon dioxide (CO2) emission per capita for the same year.
#' @format This data frame consists of 28 countries and the following variables:\cr
#' \describe{
#' \item{\code{GNP}}{The gross product per capita in 1996.}
#' \item{\code{CO2}}{The estimated carbon dioxide emission per capita in 1996.}
#' \item{\code{country}}{An abbreviation pertaining to the country measures (e.g. \code{"GRC"} = Greece and \code{"CH"} = Switzerland).}
#' }
#' @references Hurn, M., Justel, A. and Robert, C. P. (2003) Estimating mixtures of regressions, \emph{Journal of Computational and Graphical Statistics}, 12(1): 55-79.
#' @examples 
#' data(CO2data, package="MoEClust")
#' plot(CO2data$GNP, CO2data$CO2, type="n", ylab=expression('CO'[2]))
#' text(CO2data$GNP, CO2data$CO2, CO2data$country)
#' @docType data
#' @keywords datasets
#' @usage data(CO2data)
"CO2data"
