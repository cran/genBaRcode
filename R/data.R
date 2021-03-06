#' Barcode distribution of an example experiment.
#'
#' A dataset containing an example BCdat object which consists of 98 barcode sequences and with no error correction yet.
#'
#'BC_dat:
#'
#' @format A S4 data object with the following slots:
#' \describe{
#'   \item{class}{sequence overview}
#'   \item{barcode read counts}{a data frame consisting of read counts and barcode sequences}
#'   \item{results dir}{path to a directory for any kind of results}
#'   \item{barcode backbone}{a string clarifying the barcode backbone structure}
#'   \item{label}{character string, used as label for file names etc.}
#' }

"BC_dat"

#' Barcode distribution of an example experiment.
#'
#' A dataset containing an example BCdat object after error-correction which consists of 10 barcode sequences.
#'
#'BC_dat_EC:
#'
#' @format A S4 data object with the following slots:
#' \describe{
#'   \item{class}{sequence overview}
#'   \item{barcode read counts}{a data frame consisting of read counts and barcode sequences}
#'   \item{results dir}{path to a directory for any kind of results}
#'   \item{barcode backbone}{a string clarifying the barcode backbone structure}
#'   \item{label}{character string, used as label for file names etc.}
#' }

"BC_dat_EC"
