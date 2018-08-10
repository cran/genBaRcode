
#' @title Error Correction
#'
#' @description Corrects a list of equally long (barcode) sequences. Based on calculated hamming distances
#' as a measure of similarity, highly similar sequences are clustered together and the cluster label will
#' be the respective sequence with the highest read count.
#'
#' @param BC_dat one or a list of BCdat objects, containing the necessary sequences.
#' @param maxDist an integer value representing the maximal hamming distance for which it is allowed to
#' cluster two sequences together.
#' @param save_it a logical value. If TRUE the data will be saved as csv-file.
#' @param cpus an integer value, in case multiple BCdat objects are provided a CPU number greater than one
#' would allow for a parallelized calculation (one CPU per BCdat object).
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information)
#'
#' @export
#' @examples errorCorrection(BC_dat, maxDist = 8, save_it = FALSE, m = "hamming")

errorCorrection <- function(BC_dat, maxDist, save_it = FALSE, cpus = 1, m = "hamming") {

  if(cpus < 1) {
    warning("# cpus needs to be >= 1")
    cpus <- 1
  }

  if(cpus > parallel::detectCores()) {
    warning(paste("# available number of CPUs is", parallel::detectCores()))
    cpus <- parallel::detectCores() - 1
  }

  if(is.list(BC_dat)) {
      return(errorCorrection_multiple(BC_dat, maxDist, save_it, cpus, m))
  } else {
      if(cpus > 1) {
        warning("parallelization only possible for mutliple BCdat objects.")
      }
      return(errorCorrection_single(BC_dat, maxDist, save_it, m))
  }
}



#' @title Error Correction for multiple objects
#'
#' @description Corrects a list of equally long (barcode) sequences. Based on calculated hamming distances
#' as a measure of similarity, highly similar sequences are clustered together and the cluster label will
#' be the respective sequence with the highest read count.
#'
#' @param BC_dat a BCdat object, containing the necessary sequences.
#' @param maxDist an integer value representing the maximal (hamming) distance for which it is allowed to
#' cluster two sequences together.
#' @param save_it a logical value. If TRUE the data will be saved as csv-file.
#' @param cpus an integer value, in case multiple BCdat objects are provided a CPU number greater than one
#' would allow for a parallelized calculation (one CPU per BCdat object).
#' @param m a character string. Method for distance calculation, default is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information)
#'
#' @return a list of BCdat objects.

errorCorrection_multiple <- function(BC_dat, maxDist, save_it = FALSE, cpus = 1, m = "hamming") {

  if(length(BC_dat) < cpus) {
    cpus <- length(BC_dat)
  }

  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)

  # there needs to be a separately defined variable 'i' because the package check function seems to not recognise the foreach i
  i <- NULL
  tmp <- foreach::foreach(i = 1:length(BC_dat)) %dopar% {
                    if(is.list(BC_dat[[i]])) {
                      tmp <- list()
                      for(j in 1:length(BC_dat[[i]])) {
                        tmp[[j]] <- errorCorrection_single(BC_dat[[i]][[j]], maxDist = maxDist, save_it = save_it, m = m)
                      }
                      tmp
                    } else {
                      errorCorrection_single(BC_dat[[i]], maxDist = maxDist, save_it = save_it, m = m)
                    }
  }

  parallel::stopCluster(cl)
  return(tmp)

}

errorCorrection_single <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming") {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")[, 2:3]
  dat <- dat[order(dat[, 1], decreasing = FALSE), ]

  merging_barcode <- 1
  stop <- dim(dat)[1]

  while(merging_barcode < stop) {

      HD <- 1
      hammDist <- stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method = m)

      while(HD <= maxDist) {
        HDs <- c(rep(FALSE, merging_barcode), hammDist == HD)

        if(sum(HDs) > 0) {
          index <- which(HDs)[1]
          dat[index,1] <- dat[index,1] + dat[merging_barcode, 1]

          dat <- dat[-merging_barcode, ]
          dat <- dat[order(dat[, 1]), ]

          stop <- dim(dat)[1]
          HD <- maxDist + 1
        } else {
          if(HD < maxDist) {
            HD <- HD + 1
          } else {
            merging_barcode <- merging_barcode + 1
            HD <- maxDist + 1
          }
        }
      }
  }

  dat <- dat[order(dat[, 1], decreasing=TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  methods::slot(BC_dat, "reads") <- data.frame(pos = 1:dim(dat)[1], read_count = dat[, 1], barcode = dat[, 2])

  return(BC_dat)
}
