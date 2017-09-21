
#' @title Error Correction
#'
#' @description Corrects a list of equally long (barcode) sequences. Based on calculated hamming distances
#' as a measure of similarity, highly similar sequences are clustered together and the cluster label will
#' be the respective sequence with the highest read count.
#'
#' @param BC_dat a BCdat object, containing the necessary sequences.
#' @param maxHD an integer value representing the maximal hamming distance for which it is allowed to
#' cluster two sequences together.
#' @param save_it a logical value. If TRUE the data will be saved as csv-file.
#'
#' @export
#' @examples errorCorrection(BC_dat, maxHD = 8, save_it = FALSE)

errorCorrection <- function(BC_dat, maxHD, save_it = FALSE) {

    dat <- methods::slot(BC_dat, "reads")[, 2:3]
    dat <- dat[order(dat[, 1], decreasing = FALSE), ]

    merging_barcode <- 1
    stop <- dim(dat)[1]

    while(merging_barcode < stop) {

      HD <- 1
      hammDist <- stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method="h")

      while(HD <= maxHD) {
        HDs <- c(rep(FALSE, merging_barcode), hammDist == HD)

        if(sum(HDs) > 0) {
          index <- which(HDs)[1]
          dat[index,1] <- dat[index,1] + dat[merging_barcode, 1]

          dat <- dat[-merging_barcode, ]
          dat <- dat[order(dat[, 1]), ]

          stop <- dim(dat)[1]
          HD <- maxHD + 1
        } else {
          if(HD < maxHD) {
            HD <- HD + 1
          } else {
            merging_barcode <- merging_barcode + 1
            HD <- maxHD + 1
          }
        }
      }
    }

    dat <- dat[order(dat[, 1], decreasing=TRUE), ]
    if(save_it) {
      utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_ecBCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    methods::slot(BC_dat, "reads") <- data.frame(pos = 1:dim(dat)[1], read_count = dat[, 1], barcode = dat[, 2])

    return(BC_dat)
}
