
#' @title Generating Time Series Data Object
#' @description Generates a matrix containing barcodes sequences as rows and consecutive measurements at columns. It serves as
#' the necessary data object for the plotting function `plotTimeSeries`.
#'
#' @param BC_dat_list a list of BCdat objects.
#'
#' @return a data.frame containing every identified barcode and its read count per time point/measurement.
#'
#' @export

generateTimeSeriesData <- function(BC_dat_list) {

  if (length(BC_dat_list) < 2) {
    stop("# data object needs to be at least a list of length = 2")
  }

  for(i in BC_dat_list) {
    if (methods::is(i)[1] != "BCdat") {
      stop("# All list elements need to be of type BCdat.")
    }
  }

  BCs <- NULL
  for (i in 1:length(BC_dat_list)) {
      BCs <- unique(c(BCs, as.character(methods::slot(BC_dat_list[[i]], "reads")$"barcode")))
  }

  dat <- matrix(0, ncol = length(BC_dat_list), nrow = length(BCs), dimnames = list(BCs, unlist(lapply(BC_dat_list, function(x) { methods::slot(x, "label") }))))

  for (i in 1:length(BC_dat_list)) {
      index <- match(BCs, as.character(methods::slot(BC_dat_list[[i]], "reads")$"barcode"))
      dat[!is.na(index), i] <- methods::slot(BC_dat_list[[i]], "reads")$"read_count"[index[!is.na(index)]]
  }

  return(dat)
}

