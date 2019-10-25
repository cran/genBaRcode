
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
#' @param strategy since the future package is used for parallelisation a strategy has to be stated, the default is "sequential"  (cpus = 1) and "multiprocess" (cpus > 1). It is not necessary to chose a certain strategy, since it will be adjusted accordingly to the number of cpus which were choosen. For further information please read future::plan() R-Documentation.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information)
#' @param type there are different error correction strategies avalable ("standard", "connectivity based", "graph based", "clustering").
#' @param EC_analysis a logical value. If TRUE additional error correction details will be returned, which can also be visualised with the respective "error correction" plots.
#' @param only_EC_BCs a logical value. If TRUE only informations about barcodes which are still present after error correction will be saved. Only meaningful if EC_analysis is set to TRUE.
#' @param start_small a logical value. If TRUE, the error correcton type "standard" will cluster always the smallest highly similar BC with the BC of interest. IF FALSE, the error correcton type "standard" will adapt its cluster strategy and cluster always BC of interest with the most frequent highly similar BC.
#'
#' @export
#' @examples
#' data(BC_dat)
#' BC_dat_EC <- errorCorrection(BC_dat, maxDist = 8, save_it = FALSE, m = "hamming")

errorCorrection <- function(BC_dat, maxDist, save_it = FALSE, cpus = 1, strategy = "sequential", m = "hamming", type = "standard", only_EC_BCs = TRUE, EC_analysis = FALSE, start_small = TRUE) {

  if(length(type) > 1) {
    warning("# errorCorrection() allows only for one type statement - therefore type will be set to your first choice.")
    type <- type[1]
  }

  if(cpus < 1) {
    warning("# the cpus-parameter needs to be >= 1")
    cpus <- 1
    strategy <- "sequential"
  }

  if(cpus > future::availableCores()) {
    warning(paste("# available number of CPUs is", future::availableCores()))
    cpus <- future::availableCores() - 1
  }

  if(strategy == "sequential" & cpus > 1) {
    strategy <- "multiprocess"
  }
  if(strategy != "sequential" & cpus == 1) {
    strategy <- "sequential"
  }

  params <- list(maxDist = maxDist,
                 save_it = save_it,
                 m = m,
                 EC_analysis = EC_analysis,
                 nt = ifelse(strategy == "sequential", ifelse(is.null(getOption("sd_num_thread")), 1, getOption("sd_num_thread")), 1))

  if(type == "standard") {
    ecf <- errorCorrection_single_variation
    params <- c(params, start_small = start_small)
  } else {
    if(type == "connectivity based") {
      ecf <- errorCorrection_single_connections
    } else {
      if(type == "graph based") {
        ecf <- errorCorrection_single_graphComp
      } else {
        if(type == "clustering") {
          ecf <- errorCorrection_single_clustering_absolute
        } else {
          warning("# Unknown type chosen, proceed with 'standard' type.")
          ecf <- errorCorrection_single_variation
          params <- c(params, start_small = start_small)
        }
      }
    }
  }

  if(is.list(BC_dat) & cpus > 1) {
      return(errorCorrection_multiple(BC_dat, cpus = cpus, strategy = strategy, func = ecf, params = params))
  } else {
      # if(cpus > 1) {
      #   warning("# parallelization only possible for mutliple BCdat objects.")
      # }
  # if processingRawData did run with multiple input NGS files and at the same time with multiple backbones the resulting object will be
  # a list of lists hence the two is.list if-loops
      if(is.list(BC_dat)) {
        if(is.list(BC_dat[[1]])) {
            tmp <- lapply(BC_dat, function(x) {
                          tmp <- lapply(x, function(y) {
                                    do.call(ecf, c(BC_dat = y, params))
                                  })
            })
            return(tmp)
        } else {
            tmp <- lapply(BC_dat, function(x) {
                            do.call(ecf, c(BC_dat = x, params))
            })
            return(tmp)
        }
      } else {
            tmp <- do.call(ecf, c(BC_dat = BC_dat, params))
        return(tmp)
      }
  }
}

errorCorrection_multiple <- function(BC_dat, cpus, strategy, func, params) {

  if(length(BC_dat) < cpus) {
    cpus <- length(BC_dat)
  }

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(strategy = strategy, workers = cpus)

  tmp <- future.apply::future_lapply(1:length(BC_dat), function(i, BC_dat, params, func) {
    if(is.list(BC_dat[[i]])) {
      tmp <- list()
      for(j in 1:length(BC_dat[[i]])) {
        tmp[[j]] <- do.call(func, c(BC_dat = BC_dat[[i]][[j]], params))
      }
      tmp
    } else {
      do.call(func, c(BC_dat = BC_dat[[i]], params))
    }
  }, BC_dat, params, func)

  return(tmp)
}


errorCorrection_multiple_old <- function(BC_dat, maxDist, save_it = FALSE, cpus = 1, strategy, m = "hamming", EC_analysis = FALSE, func, params) {

  if(length(BC_dat) < cpus) {
    cpus <- length(BC_dat)
  }

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(strategy = strategy, workers = cpus)

  tmp <- future.apply::future_lapply(1:length(BC_dat), function(i, BC_dat, maxDist, save_it, m) {
                if(is.list(BC_dat[[i]])) {
                      tmp <- list()
                      for(j in 1:length(BC_dat[[i]])) {
                            tmp[[j]] <- func(BC_dat = BC_dat[[i]][[j]], maxDist = maxDist, save_it = save_it, m = m, EC_analysis = EC_analysis)
                      }
                      tmp
                } else {
                      func(BC_dat = BC_dat[[i]], maxDist = maxDist,
                                                 save_it = save_it, m = m, EC_analysis = EC_analysis)
                }
  }, BC_dat, maxDist, save_it, m, EC_analysis = EC_analysis)

  return(tmp)
}

errorCorrection_single_variation <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, start_small = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  if (start_small) {
    func <- min
  } else {
    func <- max
  }

  # dat could be dismissed
  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = FALSE), ]

  if (EC_analysis) {
    datEC <- as.list(as.character(dat$barcode))
    datEC_reads <- as.list(as.numeric(dat$read_count))
    datEC_index_list <- as.character(dat$barcode)
    from_part <- to_part <- NULL
  }

  bcs <- dat[, 1]
  seqs <- as.character(dat[, 2])

  merging_barcode <- 1

  basesOccurences <- t(matrix(unlist(lapply(seqs, function(x) {
    c(sum(unlist(strsplit(x, split = "")) == "A") -
        sum(unlist(strsplit(x, split = "")) == "T"),
      sum(unlist(strsplit(x, split = "")) == "G") -
        sum(unlist(strsplit(x, split = "")) == "C"),
      sum(unlist(strsplit(x, split = "")) == "N"))
  })), ncol = 3, byrow = TRUE))
  basesOccurences <- list(basesOccurences[1,], basesOccurences[2,], basesOccurences[3,])

  keep <- rep(TRUE, length(bcs))

  while (merging_barcode < length(bcs)) {
    candidates <- (abs(basesOccurences[[1]][(merging_barcode + 1):length(bcs)] - basesOccurences[[1]][merging_barcode]) +
                     abs(basesOccurences[[2]][(merging_barcode + 1):length(bcs)] - basesOccurences[[2]][merging_barcode]) +
                     abs(basesOccurences[[3]][(merging_barcode + 1):length(bcs)] - basesOccurences[[3]][merging_barcode])) / 2

    distances <- rep(NA, maxDist)
    index <- NA
    for (d in 0:maxDist) {
      candIndices <- which(candidates == d)
      if (length(candIndices) > 0) {
        candDists <- stringdist::stringdist(seqs[merging_barcode], seqs[merging_barcode + candIndices], method = m, nthread = nt)
        distances[min(candDists)] <- func(distances[min(candDists)], func(candIndices[candDists == min(candDists)], na.rm = TRUE), na.rm = TRUE)
      }
      if (d != 0 && !is.null(distances[d])) {
        index <- distances[d]
        break
      }
    }

    if (is.na(index)) {
      merging_barcode <- merging_barcode + 1
      next()
    }

    nIndex <- index + merging_barcode
    if (EC_analysis) {

      ECindex_merg <- which(as.character(seqs[merging_barcode]) == datEC_index_list)
      ECindex_target <- which(as.character(seqs[nIndex]) == datEC_index_list)

      datEC[[ECindex_target]] <- c(datEC[[ECindex_target]], datEC[[ECindex_merg]])
      datEC_reads[[ECindex_target]] <- c(datEC_reads[[ECindex_target]], datEC_reads[[ECindex_merg]])

      to_part <- c(to_part, seqs[merging_barcode])
      from_part <- c(from_part, seqs[nIndex])
    }
    bcs[nIndex] <- bcs[nIndex] + bcs[merging_barcode]

    ## new order, in case the read count of the newly summed up BC is greater then the right-hand side BC
    if (nIndex < length(bcs) && bcs[nIndex] > bcs[nIndex + 1]) {
      tmpBcs <- bcs[nIndex]
      tmpSeqs <- seqs[nIndex]
      tmpBasesOccurences <- c(basesOccurences[[1]][nIndex], basesOccurences[[2]][nIndex], basesOccurences[[3]][nIndex])

      sIndex <- nIndex + which.min(c(bcs[(nIndex + 1):length(bcs)], Inf) < bcs[nIndex]) - 1

      fis <- nIndex:(sIndex - 1)
      lis <- (nIndex + 1):sIndex

      bcs[fis] <- bcs[lis]
      bcs[sIndex] <- tmpBcs

      seqs[nIndex:(sIndex - 1)] <- seqs[(nIndex + 1):sIndex]
      seqs[sIndex] <- tmpSeqs

      for (i in 1:length(basesOccurences)) {
        basesOccurences[[i]][fis] <- basesOccurences[[i]][lis]
        basesOccurences[[i]][sIndex] <- tmpBasesOccurences[i]
      }

    }

    keep[merging_barcode] <- FALSE
    merging_barcode <- merging_barcode + 1
  }

  dat <- data.frame(read_count = rev(bcs[keep]), barcode = rev(seqs[keep]))
  if (save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep = ""), sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {
    if(only_EC_BCs) {
      index <- from_part %in% dat[, 2]
      from_part <- from_part[index]
      to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = rev(dat[, 2])),
                    data.frame(from = rev(from_part), to = rev(to_part)))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = datEC,
                  EC_reads = datEC_reads)

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
}


errorCorrection_single_biggest <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = FALSE), ]

  merging_barcode <- 1
  stop <- dim(dat)[1]

  if(EC_analysis) {
    datEC <- as.list(as.character(dat$barcode))
    datEC_reads <- as.list(as.numeric(dat$read_count))
    datEC_index_list <- as.character(dat$barcode)

    from_part <- to_part <- NULL
  }

  while(merging_barcode < stop) {

    hammDist <- stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method = m, nthread = nt)
    index <- max(which(hammDist == min(hammDist)))

    if(hammDist[index] <= maxDist) {

        if(EC_analysis) {
          ECindex_merg <- which(as.character(dat[merging_barcode, 2]) == datEC_index_list)
          ECindex_target <- which(as.character(dat[index + merging_barcode, 2]) == datEC_index_list)

          datEC[[ECindex_target]] <- c(datEC[[ECindex_target]], datEC[[ECindex_merg]])
          datEC_reads[[ECindex_target]] <- c(datEC_reads[[ECindex_target]], datEC_reads[[ECindex_merg]])

          to_part <- c(to_part, as.character(dat[merging_barcode, 2]))
          from_part <- c(from_part, as.character(dat[merging_barcode + index, 2]))
        }

        dat[index + merging_barcode, 1] <- dat[index + merging_barcode, 1] + dat[merging_barcode, 1]
        dat <- dat[-merging_barcode, ]
        dat <- dat[order(dat[, 1]), ]
        stop <- dim(dat)[1]
    } else {
        merging_barcode <- merging_barcode + 1
    }
  }

  dat <- dat[order(dat[, 1], decreasing = TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {

    if(only_EC_BCs) {
        index <- from_part %in% dat[, 2]
        from_part <- from_part[index]
        to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = rev(dat[, 2])),
                    data.frame(from = rev(from_part), to = rev(to_part)))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = datEC,
                  EC_reads = datEC_reads)

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
  return(BC_dat)
}

errorCorrection_single <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = FALSE), ]

  merging_barcode <- 1
  stop <- dim(dat)[1]

  if(EC_analysis) {
    datEC <- as.list(as.character(dat$barcode))
    datEC_reads <- as.list(as.numeric(dat$read_count))
    datEC_index_list <- as.character(dat$barcode)

    from_part <- to_part <- NULL
  }

  while(merging_barcode < stop) {

    hammDist <- stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method = m, nthread = nt)
    index <- which.min(hammDist)

    if(hammDist[index] <= maxDist) {

      if(EC_analysis) {
        ECindex_merg <- which(as.character(dat[merging_barcode, 2]) == datEC_index_list)
        ECindex_target <- which(as.character(dat[index + merging_barcode, 2]) == datEC_index_list)

        datEC[[ECindex_target]] <- c(datEC[[ECindex_target]], datEC[[ECindex_merg]])
        datEC_reads[[ECindex_target]] <- c(datEC_reads[[ECindex_target]], datEC_reads[[ECindex_merg]])

        to_part <- c(to_part, as.character(dat[merging_barcode, 2]))
        from_part <- c(from_part, as.character(dat[merging_barcode + index, 2]))
      }

      dat[index + merging_barcode, 1] <- dat[index + merging_barcode, 1] + dat[merging_barcode, 1]
      dat <- dat[-merging_barcode, ]
      dat <- dat[order(dat[, 1]), ]
      stop <- dim(dat)[1]
    } else {
      merging_barcode <- merging_barcode + 1
    }
  }

  dat <- dat[order(dat[, 1], decreasing = TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {

    if(only_EC_BCs) {
      index <- from_part %in% dat[, 2]
      from_part <- from_part[index]
      to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = rev(dat[, 2])),
                    data.frame(from = rev(from_part), to = rev(to_part)))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = datEC,
                  EC_reads = datEC_reads)

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
  return(BC_dat)
}

errorCorrection_single_connections <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = FALSE), ]

  if(EC_analysis) {
    datEC <- as.list(as.character(dat$barcode))
    datEC_reads <- as.list(as.numeric(dat$read_count))
    datEC_index_list <- as.character(dat$barcode)

    from_part <- to_part <- NULL
  }

  connections <- colSums(stringdist::stringdistmatrix(dat[, 2], dat[, 2], method = m, nthread = nt) == 1)
  dat <- cbind(dat, connections)
  d <- matrix(NA, ncol = 3, nrow = 1)
  colnames(d) <- colnames(dat)
  for(i in sort(unique(connections))) {
    tmp <- dat[connections == i, ]
    tmp <- tmp[order(tmp[, 1]), ]
    d <- rbind(d, tmp)
  }
  d <- d[-1, ]

  merging_barcode <- which.max(d$connections != 0)
  stop <- dim(dat)[1]

  dat <- d[, 1:2]
  while(merging_barcode < stop) {

    hammDist <- stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method = m, nthread = nt)
    index <- which.min(hammDist)

    if(hammDist[index] <= maxDist) {

      if(EC_analysis) {
        ECindex_merg <- which(as.character(dat[merging_barcode, 2]) == datEC_index_list)
        ECindex_target <- which(as.character(dat[index + merging_barcode, 2]) == datEC_index_list)

        datEC[[ECindex_target]] <- c(datEC[[ECindex_target]], datEC[[ECindex_merg]])
        datEC_reads[[ECindex_target]] <- c(datEC_reads[[ECindex_target]], datEC_reads[[ECindex_merg]])

        to_part <- c(to_part, as.character(dat[merging_barcode, 2]))
        from_part <- c(from_part, as.character(dat[merging_barcode + index, 2]))
      }

      dat[index + merging_barcode, 1] <- dat[index + merging_barcode, 1] + dat[merging_barcode, 1]
      dat <- dat[-merging_barcode, ]
      #dat <- dat[order(dat[, 1]), ]
      stop <- dim(dat)[1]
    } else {
      merging_barcode <- merging_barcode + 1
    }
  }

  dat <- dat[order(dat[, 1], decreasing = TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {

    if(only_EC_BCs) {
      index <- from_part %in% dat[, 2]
      from_part <- from_part[index]
      to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = rev(dat[, 2])),
                    data.frame(from = rev(from_part), to = rev(to_part)))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = datEC,
                  EC_reads = datEC_reads)

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
  return(BC_dat)
}

errorCorrection_single_clustering_absolute <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = TRUE), ]

  if(EC_analysis) {
    datEC <- as.list(as.character(dat$barcode))
    datEC_reads <- as.list(as.numeric(dat$read_count))
    datEC_index_list <- as.character(dat$barcode)

    from_part <- to_part <- NULL
  }

  merging_barcode <- 1
  while(merging_barcode < dim(dat)[1]) {

      index <- which(stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method = m, nthread = nt) <= maxDist)

      if(length(index) != 0) {

        if(EC_analysis) {
          ECindex_target <- merging_barcode
          ECindex_merg <- which(datEC_index_list %in% as.character(dat[index + merging_barcode, 2]))

          datEC[[ECindex_target]] <- c(datEC[[ECindex_target]], unlist(datEC[ECindex_merg]))
          datEC_reads[[ECindex_target]] <- c(datEC_reads[[ECindex_target]], unlist(datEC_reads[ECindex_merg]))

          from_part <- c(from_part, rep(as.character(dat[merging_barcode, 2]), length(index)))
          to_part <- c(to_part, as.character(dat[merging_barcode + index, 2]))
        }

        dat[merging_barcode, 1] <- sum(dat[index + merging_barcode, 1], dat[merging_barcode, 1])
        dat <- dat[-(index + merging_barcode), ]
        #dat <- dat[order(dat[, 1]), ]
        stop <- dim(dat)[1]
      }
      merging_barcode <- merging_barcode + 1
  }

  dat <- dat[order(dat[, 1], decreasing = TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {

    if(only_EC_BCs) {
      index <- from_part %in% dat[, 2]
      from_part <- from_part[index]
      to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = rev(dat[, 2])),
                    data.frame(from = rev(from_part), to = rev(to_part)))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = rev(datEC),
                  EC_reads = rev(datEC_reads))

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
  return(BC_dat)
}


errorCorrection_single_clustering_stepwise <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = TRUE), ]

  if(EC_analysis) {
    datEC <- as.list(as.character(dat$barcode))
    datEC_reads <- as.list(as.numeric(dat$read_count))
    datEC_index_list <- as.character(dat$barcode)

    from_part <- to_part <- NULL
  }

  for(i in 1:maxDist) {

    merging_barcode <- 1
    stop <- dim(dat)[1]

    dat <- dat[order(dat[, 1], decreasing = TRUE), ]

    while(merging_barcode < stop) {

      index <- which(stringdist::stringdist(dat[merging_barcode, 2], dat[(merging_barcode+1):dim(dat)[1], 2], method = m, nthread = nt) == i)

      if(length(index) != 0) {

        if(EC_analysis) {
          ECindex_target <- which(as.character(dat[merging_barcode, 2]) == datEC_index_list)
          ECindex_merg <- which(datEC_index_list %in% as.character(dat[index + merging_barcode, 2]))

          datEC[[ECindex_target]] <- c(datEC[[ECindex_target]], unlist(datEC[ECindex_merg]))
          datEC_reads[[ECindex_target]] <- c(datEC_reads[[ECindex_target]], unlist(datEC_reads[ECindex_merg]))

          from_part <- c(from_part, rep(as.character(dat[merging_barcode, 2]), length(index)))
          to_part <- c(to_part, as.character(dat[merging_barcode + index, 2]))
        }

        dat[merging_barcode, 1] <- sum(dat[index + merging_barcode, 1], dat[merging_barcode, 1])
        dat <- dat[-(index + merging_barcode), ]
        #dat <- dat[order(dat[, 1]), ]
        stop <- dim(dat)[1]
      }
      merging_barcode <- merging_barcode + 1
    }

  }

  dat <- dat[order(dat[, 1], decreasing = TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {

    if(only_EC_BCs) {
      index <- from_part %in% dat[, 2]
      from_part <- from_part[index]
      to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = rev(dat[, 2])),
                    data.frame(from = rev(from_part), to = rev(to_part)))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = rev(datEC),
                  EC_reads = rev(datEC_reads))

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
  return(BC_dat)
}

errorCorrection_single_graphComp <- function(BC_dat, maxDist, save_it = FALSE, m = "hamming", EC_analysis = FALSE, only_EC_BCs = TRUE, nt) {

  methods::slot(BC_dat, "label") <- paste0(methods::slot(BC_dat, "label"), "_EC")

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(BC_dat)
  }

  dat <- methods::slot(BC_dat, "reads")
  dat <- dat[order(dat[, 1], decreasing = TRUE), ]

  if(EC_analysis) {
    datEC <- list()
    datEC_reads <- list()

    from_part <- to_part <- NULL
  }

  adj_dat <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m, nthread = nt) <= maxDist
  colnames(adj_dat) <- methods::slot(BC_dat, "reads")$"barcode"
  g <- igraph::graph.adjacency(adj_dat, mode = c("undirected"), diag = FALSE, add.colnames = "names")
  comps <- igraph::components(g)

  seqs <- reads <- NULL

  for(i in 1:comps$no) {

    index <- which(comps$membership == i)

    seqs <- c(seqs, as.character(dat[index[1], 2]))
    reads <- c(reads, sum(as.numeric(dat[index, 1])))

    if(EC_analysis) {
          datEC[[index[1]]] <- as.character(dat[index, 2])
          datEC_reads[[index[1]]] <- as.numeric(dat[index, 2])

          to_part <- c(to_part, rep(as.character(dat[index[1], 2]), length(index[-1])))
          from_part <- c(from_part, as.character(dat[index[-1], 2]))
    }
  }

  dat <- data.frame(read_count = reads, barcode = seqs)
  dat <- dat[order(dat[, 1], decreasing = TRUE), ]
  if(save_it) {
    utils::write.table(dat, paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  if(EC_analysis) {

    if(only_EC_BCs) {
      index <- from_part %in% dat[, 2]
      from_part <- from_part[index]
      to_part <- to_part[index]
    }

    fromTo <- rbind(data.frame(from = "origin", to = dat[, 2]),
                    data.frame(from = rev(from_part), to = to_part))

    vertices <- data.frame(barcodes = c("origin", as.character(methods::slot(BC_dat, "reads")[, 2])),
                           read_counts = c(1, as.numeric(methods::slot(BC_dat, "reads")[, 1])))

    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    final <- list(BC_dat = BC_dat,
                  edges = fromTo,
                  vertices = vertices,
                  EC_seqs = rev(datEC),
                  EC_reads = rev(datEC_reads))

    if(save_it) {
      save(final, file = paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), "_stats.rds", sep = ""))
    }
    return(final)
  } else {
    methods::slot(BC_dat, "reads") <- data.frame(read_count = dat[, 1], barcode = dat[, 2])
    return(BC_dat)
  }
  return(BC_dat)
}
