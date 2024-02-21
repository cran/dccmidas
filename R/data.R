#' S&P 500 data
#'
#' Daily data on S&P 500 collected from the realized library of
#' the Oxford-Man Institute \insertCite{heber_2009}{dccmidas}. 
#'
#' @docType data
#'
#' @usage data(sp500)
#'
#' @format An object of class \code{"xts"}.
#' @details sp500 includes the open price (open_price), the realized variance 
#' (rv5), and the close price (close_price). The realized variance has been calculated
#' using intradaily intervals of five minutes \insertCite{andersen_boll_1998}{dccmidas}. 
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @references
#' \insertAllCited{} 
#'
#' @source Realized library of the Oxford-Man Institute
#'
#' @examples
#' head(sp500)
#' summary(sp500)
"sp500"

#' FTSE 100 data
#'
#' Daily data on FTSE 100 collected from the realized library of
#' the Oxford-Man Institute \insertCite{heber_2009}{dccmidas}. 
#'
#' @docType data
#'
#' @usage data(ftse100)
#'
#' @format An object of class \code{"xts"}.
#' @details ftse100 includes the open price (open_price), the realized variance 
#' (rv5), and the close price (close_price). The realized variance has been calculated
#' using intradaily intervals of five minutes \insertCite{andersen_boll_1998}{dccmidas}. 
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @references
#' \insertAllCited{} 
#'
#' @source Realized library of the Oxford-Man Institute
#'
#' @examples
#' head(ftse100)
#' summary(ftse100)
"ftse100"

#' NASDAQ data
#'
#' Daily data on NASDAQ collected from the realized library of
#' the Oxford-Man Institute \insertCite{heber_2009}{dccmidas}. 
#'
#' @docType data
#'
#' @usage data(nasdaq)
#'
#' @format An object of class \code{"xts"}.
#' @details nasdaq includes the open price (open_price), the realized variance 
#' (rv5), and the close price (close_price). The realized variance has been calculated
#' using intradaily intervals of five minutes \insertCite{andersen_boll_1998}{dccmidas}. 
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @references
#' \insertAllCited{} 
#'
#' @source Realized library of the Oxford-Man Institute
#'
#' @examples
#' head(nasdaq)
#' summary(nasdaq)
"nasdaq"

#' Monthly U.S. Industrial Production
#'
#' Monthly data on the U.S. Industrial Production index (IP, index 2012=100, seasonally adjusted) collected from the 
#' Federal Reserve Economic Data (FRED) archive. The IP has been used as MIDAS term in different contributions  
#' (see, for instance, \insertCite{engle_ghysels_sohn_2013;textual}{rumidas}, \insertCite{conrad_lock_2015;textual}{rumidas}, and
#' \insertCite{amendola_candila_scognamillo_2017;textual}{rumidas}).
#'
#' @docType data
#'
#' @usage data(indpro)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @references
#' \insertAllCited{} 
#'
#' @source Archive of the Federal Reserve Economic Data \href{https://fred.stlouisfed.org/series/INDPRO}{(FRED)}
#'
#' @examples
#' head(indpro)
#' summary(indpro)
#' plot(indpro)
"indpro"


