
#' @title Data processing
#'
#' @description Reads the corresponding fast(q) file(s), extracts the defined barcode constructs and counts them. Optionally,
#' a Phred-Score based quality filtering will be conducted and the results will be saved within a csv file.
#' @param file_name a character string or a character vector, containing the file name(s).
#' @param source_dir a character string which contains the path to the source files.
#' @param results_dir a character string which contains the path to the results directory.
#' @param mismatch an positive integer value, default is 0, if greater values are provided they indicate the number of allowed mismtaches when identifying the barcode constructes.
#' @param indels a logical value. If TRUE the chosen number of mismatches will be interpreted as edit distance and allow for insertions and deletions as well.
#' @param label a character string which serves as a label for every kind of created output file.
#' @param bc_backbone a character string describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param bc_backbone_label a character vector, an optional list of barcode backbone names serving as additional identifier within file names and BCdat labels. If not provided ordinary numbers will serve as alternative.
#' @param min_score a positive integer value, all fastq sequence with an average score smaller
#'  then min_score will be excluded, if min_score = 0 there will be no quality score filtering
#' @param min_reads positive integer value, all extracted barcode sequences with a read count smaller than min_reads will be excluded from the results
#' @param save_it a logical value. If TRUE, the raw data will be saved as a csv-file.
#' @param seqLogo a logical value. If TRUE, the sequence logo of the entire NGS file will be generated and saved.
#' @param cpus an integer value, indicating the number of available cpus.
#' @param full_output a logical value. If TRUE, additional output files will be generated.
#'
#' @return a BCdat object which includes reads, seqs, directories, masks.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' bc_backbone <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
#'
#' source_dir <- system.file("extdata", package = "genBaRcode")
#'
#' processingRawData(file_name = "test_data.fastq", source_dir, results_dir = getwd(), mismatch = 0,
#' label = "test", bc_backbone, min_score = 30, indels = FALSE,
#' min_reads = 2, save_it = TRUE, seqLogo = FALSE)
#'
#' }

processingRawData <- function(file_name, source_dir, results_dir, mismatch = 0, indels = FALSE, label = "", bc_backbone, bc_backbone_label = "", min_score = 30, min_reads = 2, save_it = TRUE, seqLogo = FALSE, cpus = 1, full_output = FALSE) {

  # if(.Platform$OS.type != "unix" & unix) {
  #   warning("# no unix based system recognized, set parameter 'unix' to FALSE")
  #   unix <- FALSE
  # }

  source_dir <- .testDirIdentifier(source_dir)
  results_dir <- .testDirIdentifier(results_dir)

  if(length(bc_backbone_label) != length(bc_backbone)) {
    warning(paste(length(bc_backbone), "BC patterns but", length(bc_backbone_label),"BC pattern labels"))
    bc_backbone_label = 1:length(bc_backbone)
  }

 if(length(file_name) > 1) {
   dat <- processingRawData_multiple(file_name, source_dir, results_dir, mismatch, indels, label, bc_backbone, bc_backbone_label, min_score, min_reads, save_it, seqLogo, cpus, full_output)
 } else {
   if(length(file_name) == 1) {
      dat <- processingRawData_single(file_name, source_dir, results_dir, mismatch, indels, label, bc_backbone, bc_backbone_label, min_score, min_reads, save_it, seqLogo, full_output, cpus)
    } else {
      stop("# file_name needs to be a character string (processingRawData)")
    }
 }

  return(dat)
}

processingRawData_single <- function(file_name, source_dir, results_dir, mismatch = 0, indels = FALSE, label = "", bc_backbone, bc_backbone_label = 1:length(bc_backbone), min_score = 30, min_reads = 2, save_it = TRUE, seqLogo = FALSE, full_output, cpus) {

    if(label == "") {
        label <- strsplit(file_name, split = "[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
      }

      ending <- strsplit(file_name, split="[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
      if(min_score > 0) {
        if(ending == "fastq" | ending == "FASTQ" | ending == "fq") {
          dat <- qualityFiltering(file_name, source_dir, min_score)
          if(length(dat) == 0) {
            warning("There are no sequences after quality filtering!")
          }
        } else {
          warning("NGS quality filtering only available for FASTQ file formats! Proceeding without filtering!", call. = FALSE, immediate. = TRUE)
        }
      } else {
        if(ending == "fasta") {
          dat <- ShortRead::readFasta(dirPath = source_dir, pattern = file_name)
        } else {
          dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)
        }
      }

      if(seqLogo) {
        l<- nchar(as.character(ShortRead::sread(dat)[1]))
        ggplot2::ggsave(filename = paste("seqLogo_", label, "_NGS.png", sep = ""),
                        plot = suppressMessages(plotSeqLogo(as.character(ShortRead::sread(dat))) + ggplot2::scale_x_continuous(breaks = c(1, round(l/2), l))),
                        device = "png",
                        path = results_dir,
                        width = l / 2,
                        height = l / 16)
      }

      dat <- extractBarcodes(dat, label, results_dir, mismatch = mismatch, indels = indels, bc_backbone, full_output, cpus)

      if(length(bc_backbone) == 1) {
        dat <- prepareDatObject(dat, results_dir, label, bc_backbone, min_reads, save_it)
      } else {
        for(b in 1:length(bc_backbone)) {
          dat[[b]] <- prepareDatObject(dat[[b]], results_dir, label = paste(label, bc_backbone_label[b], sep = "_"), bc_backbone[b], min_reads, save_it)
        }
      }

      return(dat)
}

#' @importFrom foreach %dopar%

processingRawData_multiple <- function(file_name, source_dir, results_dir, mismatch = 0, indels = FALSE, label = "", bc_backbone, bc_backbone_label = 1:length(bc_backbone), min_score = 30, min_reads = 2, save_it = TRUE, seqLogo = FALSE, cpus = 1, full_output) {

      if(cpus < 1) {
        warning("# cpus needs to be >= 1")
        cpus <- 1
      }

      if(cpus > parallel::detectCores()) {
        warning(paste("# available number of CPUs is", parallel::detectCores()))
        cpus <- parallel::detectCores() - 1
      }

      if(length(label) != length(file_name)) {
        if(label != "") {
          message(paste0("# there are ", length(file_name), " files but ", length(label), "labels"))
        }

        label <- unlist(lapply(strsplit(file_name, split = "[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE), function(x) { x[1] }))
      }

      if(cpus > length(file_name)) {
        cpus <- length(file_name)
      }

      cl <- parallel::makeCluster(cpus)
      doParallel::registerDoParallel(cl)

      # due to dopar problems
      i <- NULL
      tmp <- foreach::foreach(i = 1:length(file_name)) %dopar% {
        processingRawData_single(file_name[i], source_dir, results_dir, mismatch, indels, label[i], bc_backbone, bc_backbone_label, min_score, min_reads, save_it, seqLogo, full_output, cpus)
      }

      parallel::stopCluster(cl)
      return(tmp)
}



#' Data Object Preparation
#'
#' generates BCdat object after barcode backbone identification.
#'
#' @param dat a tbl_df object (e.g. created by dplyr::count)
#' @param results_dir a character string which contains the path to the results directory.
#' @param label a character string which serves as a label for every kind of created output file.
#' @param bc_backbone a character string describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param min_reads positive integer value, all extracted barcode sequences with a read count smaller than min_reads will be excluded from the results
#' @param save_it a logical value. If TRUE, the raw data will be saved as a csv-file.
#'
#' @return a BCdat object.
#' @importFrom methods new
#'
prepareDatObject <- function(dat, results_dir, label, bc_backbone, min_reads, save_it) {

  if(dim(dat)[1] == 0) {
      return(methods::new(Class = "BCdat", reads = data.frame(), results_dir = results_dir,
                   label = label,
                   mask = bc_backbone))
  }

  dat <- dat[dat[, 2] > min_reads, ]
  if(dim(dat)[1] != 0) {
    dat <- data.frame(pos = 1:dim(dat)[1], read_count = as.numeric(dat$n), barcode = as.character(unlist(dat[, 1])))
    if(save_it) {
      if(sum(dim(dat)) > 2) {
        utils::write.table(dat[, 2:3], paste(results_dir, label, "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
      } else {
        utils::write.table(dat, paste(results_dir, label, "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
    }

    return(methods::new(Class = "BCdat", reads = dat, results_dir = results_dir,
                            label = label,
                            mask = bc_backbone))
  } else {
    warning(paste("only barocdes with less than", min_reads, "reads detected for backbone:", bc_backbone))
    return(methods::new(Class = "BCdat", reads = data.frame(), results_dir = results_dir,
                          label = label,
                          mask = bc_backbone))
  }
}

#' @title Barcode extraction
#'
#' @description Extracts barcodes according to the given barcode design from a fastq file.
#'
#' @param dat a ShortReadQ object.
#' @param label a character string.
#' @param results_dir a character string which contains the path to the results directory.
#' @param mismatch an positive integer value, default is 0, if greater values are provided they indicate the number of allowed mismatches when identifing the barcode constructe.
#' @param indels under construction.
#' @param bc_backbone a character string or character vector describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param cpus an integer value, indicating the number of available cpus.
#' @param full_output a logical value. If TRUE additional output files will be generated in order to identify errors.
#'
#' @return one or a list of frequency table(s) of barcode sequences.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' bc_backbone <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
#' source_dir <- system.file("extdata", package = "genBaRcode")
#' dat <- ShortRead::readFastq(dirPath = source_dir, pattern = "test_data.fastq")
#'
#' extractBarcodes(dat, label = "test", results_dir = getwd(), mismatch = 0,
#' indels = FALSE, bc_backbone)
#'
#' }

extractBarcodes <- function(dat, label, results_dir = "./", mismatch = 0, indels = FALSE, bc_backbone, full_output = FALSE, cpus = 1) {

  if(mismatch < 0) {
    warning("# mismatch needs to be an integer >= 0")
    mismatch <- 0
  }

  if(cpus < 1) {
    warning("# cpus needs to be >= 1")
    cpus <- 1
  }

   if(length(bc_backbone) > 1) {
      extractBarcodes_multiple(dat, label, results_dir, mismatch, indels, bc_backbone, full_output)
    } else {
     extractBarcodes_single(dat, label, results_dir, mismatch, indels, bc_backbone, full_output, cpus)
    }
}

#' @title Barcode extraction for a single backbone
#'
#' @description Extracts barcodes according to the given barcode design from a fastq file.
#'
#' @param dat a ShortReadQ object.
#' @param label a character string.
#' @param results_dir a character string which contains the path to the results directory.
#' @param mismatch an positive integer value, default is 0, if greater values are provided they indicate the number of allowed mismatches when identifing the barcode constructe.
#' @param indels under construction.
#' @param bc_backbone a character string or character vector describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param cpus an integer value, indicating the number of available cpus.
#' @param full_output a logical value. If TRUE additional output files will be generated in order to identify errors.
#'
#' @return one frequency table of barcode sequences.
#'
#' @importFrom foreach %dopar%

extractBarcodes_single <- function(dat, label, results_dir, mismatch = 0, indels = FALSE, bc_backbone, full_output = FALSE, cpus = 1) {

  read_length <- ShortRead::width(dat)[1]
  mm <- 0

  match_index <- Biostrings::vmatchPattern(bc_backbone, ShortRead::sread(dat),
                             max.mismatch = mm, min.mismatch = mm,
                             with.indels = indels, fixed = FALSE,
                             algorithm = "auto")
  match_matrix <- data.frame(match_index)

  while(mm < mismatch & length(unique(match_matrix$group)) < length(ShortRead::sread(dat))) {
    mm <- mm + 1
    match_index <- Biostrings::vmatchPattern(bc_backbone, ShortRead::sread(dat),
                                             max.mismatch = mm, min.mismatch = mm,
                                             with.indels = indels, fixed = FALSE,
                                             algorithm = "auto")
    match_matrix <- rbind(match_matrix, data.frame(match_index))
  }

  if(sum(duplicated(match_matrix$group)) != 0) {
    warning("Duplicated backbone matches")
    if(full_output) {
      group <- unique(match_matrix$group[duplicated(match_matrix$group)])
      utils::write.table(match_matrix[match_matrix$group %in% group, ], paste(results_dir, label, "_", bc_backbone, "_duplicatedMatcher.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    match_matrix <- match_matrix[!duplicated(match_matrix$group), ]
  }

  if(cpus > 1) {
    cl <- parallel::makeCluster(cpus)
    doParallel::registerDoParallel(cl)

    # due to dopar problems
    s <- NULL
    parts <- seq(1, dim(match_matrix)[1], floor(dim(match_matrix)[1] / cpus))
    tmp <- foreach::foreach(s = 1:(length(parts) - 1)) %dopar% {
      apply(match_matrix[ifelse(s == 1, parts[s], parts[s] + 1):parts[s+1], ], 1, function(x) {
        x <- as.numeric(x)
        substring(as.character(ShortRead::sread(dat)[x[1]]), x[3], x[4])
      })
    }
    parallel::stopCluster(cl)
    dat <- unlist(tmp)
  } else {
    dat <- apply(match_matrix, 1, function(x) {
      x <- as.numeric(x)
      substring(as.character(ShortRead::sread(dat)[x[1]]), x[3], x[4])
    })
  }

  if(!indels & length(table(nchar(dat))) > 1) {

    if(full_output) {
      utils::write.table(dat, paste(results_dir, label, "_", bc_backbone, "_diffLengthBCs.csv", sep=""), sep=";", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    for(m in 1:mismatch) {
      if(sum(match_matrix$start == (1-m)) > 0) {
        dat[match_matrix$start == (1-m)] <- paste0(paste(rep("N", m), collapse = ""), dat[match_matrix$start == (1-m)])
      } else {
        if(sum(match_matrix$end == read_length+m) > 0) {
          dat[match_matrix$start == read_length] <- paste0(dat[match_matrix$start == 0], paste(rep("N", m), collapse = ""))
        } else {
          warning("different barcode lengths")
        }
      }
    }
  }

  if(length(dat) != 0) {
    if(!indels) {
      wobble_pos <- .getWobblePos(bc_backbone)
      dat <- data.frame(unlist(lapply(dat, function(x, wobble_pos) {
        paste(unlist(strsplit(x, split= ""))[wobble_pos], collapse = "")
      }, wobble_pos)))
    } else {
      dat <- as.data.frame(dat)
    }
    return(dplyr::count(dat, dat[, 1], sort = TRUE))
  } else {
      warning(paste("# ", label, ": no backbone patterns detectable!", sep =""))
      return(data.frame())
  }

}

#' @title Barcode extraction for multiple backbones
#'
#' @description Extracts barcodes according to the given barcode design from a fastq file.
#'
#' @param dat a ShortReadQ object.
#' @param label a character string.
#' @param results_dir a character string which contains the path to the results directory.
#' @param mismatch an positive integer value, default is 0, if greater values are provided they indicate the number of allowed mismatches when identifing the barcode constructe.
#' @param indels under construction.
#' @param bc_backbone a character string or character vector describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param full_output a logical value. If TRUE additional output files will be generated in order to identify errors.
#'
#' @return a list of frequency table(s) of barcode sequences.

extractBarcodes_multiple <- function(dat, label, results_dir, mismatch = 0, indels = FALSE, bc_backbone, full_output) {

  dat_list <- list()
  read_length <- ShortRead::width(dat)[1]

  for(i in 1:length(bc_backbone)) {

    mm <- 0

    match_index <- Biostrings::vmatchPattern(bc_backbone[i], ShortRead::sread(dat),
                                                 max.mismatch = mm, min.mismatch = mm,
                                                 with.indels = indels, fixed = FALSE,
                                                 algorithm = "auto")
    match_matrix <- data.frame(match_index)

    while(mm < mismatch & length(unique(match_matrix$group)) < length(ShortRead::sread(dat))) {
      mm <- mm + 1
      match_index <- Biostrings::vmatchPattern(bc_backbone[i], ShortRead::sread(dat),
                                               max.mismatch = mm, min.mismatch = mm,
                                               with.indels = indels, fixed = FALSE,
                                               algorithm = "auto")
      match_matrix <- rbind(match_matrix, data.frame(match_index))
    }

    if(sum(duplicated(match_matrix$group)) != 0) {
      warning("Duplicated backbone matches")
      if(full_output) {
        group <- unique(match_matrix$group[duplicated(match_matrix$group)])
        utils::write.table(match_matrix[match_matrix$group %in% group, ], paste(results_dir, label, "_", bc_backbone, "_duplicatedMatcher.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
      match_matrix <- match_matrix[!duplicated(match_matrix$group), ]
    }

    dat_list[[i]] <- apply(match_matrix, 1, function(x) {
      x <- as.numeric(x)
      substring(as.character(ShortRead::sread(dat)[x[1]]), x[3], x[4])
    })

    if(!indels & length(table(nchar(dat_list[[i]]))) > 1) {

      if(full_output) {
        utils::write.table(dat_list[[i]], paste(results_dir, label, "_", bc_backbone[i], "_diffLengthBCs.csv", sep=""), sep=";", row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
     for(m in 1:mismatch) {
        if(sum(match_matrix$start == (1-m)) > 0) {
          dat_list[[i]][match_matrix$start == (1-m)] <- paste0(paste(rep("N", m), collapse = ""), dat_list[[i]][match_matrix$start == (1-m)])
        } else {
          if(sum(match_matrix$end == read_length+m) > 0) {
            dat_list[[i]][match_matrix$start == read_length] <- paste0(dat_list[[i]][match_matrix$start == 0], paste(rep("N", m), collapse = ""))
          } else {
            warning("different barcode lengths")
          }
        }
     }
    }

    if(length(dat_list[[i]]) != 0) {
        if(!indels) {
          wobble_pos <- .getWobblePos(bc_backbone[i])
          dat_list[[i]] <- data.frame(unlist(lapply(dat_list[[i]], function(x, wobble_pos) {
            paste(unlist(strsplit(x, split= ""))[wobble_pos], collapse = "")
          }, wobble_pos)))
        } else {
          dat_list[[i]] <- as.data.frame(dat_list[[i]])
        }
        dat_list[[i]] <- dplyr::count(dat_list[[i]], dat_list[[i]][, 1], sort = TRUE)
        dat <- dat[-match_matrix$group]
    } else {
      dat_list[[i]] <- data.frame()
    }
  }

  if(sum(unlist(lapply(dat_list, nrow))) > 0) {
    return(dat_list)
  } else {
    warning(paste("# ", label, ": no backbone patterns detectable!", sep =""))
    return(data.frame())
  }
}


#' @title Identifies hybrid barcodes
#'
#' @description Experimental function to identify hybrid barcodes which can occure due to unfinished synthesis of a template
#' in-between PCR cycles.
#'
#' @param dat a character vector containing barcode sequences or a BCdat object.
#' @param min_seq_length a positive integer value indicating the minimal length of the two barcodes which give rise to a hybrid barcode.
#'
#' @return a hybrid-free frequency table of barcode sequences
#' @export
#'
#' @examples
#' data(BC_dat)
#' hybridsIdentification(BC_dat, min_seq_length = 2)

hybridsIdentification <- function(dat, min_seq_length = 2) {

  if(methods::is(dat)[1] == "BCdat") {
    dat <- methods::slot(dat, "reads")[, 2:3]
  } else {
    if(dim(dat)[2] != 2) {
      stop("# Hybrid detection needs a two-column data.frame consisting of read counts and barcode sequences")
    }
  }

  seq_length <- unique(nchar(as.character(dat[, 2])))
  if(length(seq_length) > 1) {
    stop("# Hybrid detection only for sequences of identical length!")
  }

  hybrid_result <- NULL
  hybrid_pos_index <- NULL
  num_of_seqs <- length(dat[, 2])

  for(i in num_of_seqs) {
    BC_1 <- unlist(strsplit(as.character(dat[i, 2]), split = ""))
    for(j in (1:num_of_seqs)[-i]) {
      BC_2 <- unlist(strsplit(as.character(dat[j, 2]), split = ""))
      for(k in min_seq_length:(seq_length-min_seq_length)) {
        hybrid <- paste(c(BC_1[1:k], BC_2[(k+1):seq_length]), collapse = "")
        index <- dat[, 2] == hybrid
        index[c(1:i,j, hybrid_pos_index)] <- FALSE
        if(sum(index) > 0) {
          if(sum(dat[, 2] == hybrid) > 1) { stop(" ## found same BC twice (find.hybrid.after.EC) ## ") }
          hybrid_result <- rbind(hybrid_result, c(i, which(index)), c(j, which(index)))
          hybrid_pos_index <- c(hybrid_pos_index, which(index))
        }
      }
    }
  }

  if(!is.null(hybrid_pos_index)) {
      for(i in unique(hybrid_result[, 2])) {
        ori_BCs <- hybrid_result[hybrid_result[, 2] == i, 1]
        dat[ori_BCs[1], 1] <- dat[ori_BCs[1], 1] + dat[i, 1]
        dat[ori_BCs[2], 1] <- dat[ori_BCs[2], 1] + dat[i, 1]
      }
      if(length(hybrid_result) > 0) {
        dat[, 3] <- dat[-unique(hybrid_result[, 2]), ]
      }

      dat[, 3] <- dat[order(dat[, 1], decreasing = TRUE), ]
      return(dat)
  } else {
    return(dat)
  }
}

#' @title Quality Filtering
#' @description Excludes all sequences of a given fastq file below a certain quality value.
#'
#' @param file_name a character string containing the name of the source file.
#' @param source_dir a character string containing the path to the source directory.
#' @param results_dir a character string containing the path to the directory of the results.
#' @param min_score an integer value representing the minimal average phred score a read has to achieve in order to be accepted.
#'
#' @return a ShortRead object.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' source_dir <- system.file("extdata", package = "genBaRcode")
#' qualityFiltering(file_name = "test_data.fastq", source_dir, results_dir = getwd(), min_score = 30)
#' }

qualityFiltering <- function(file_name, source_dir, results_dir, min_score = 30) {

  dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)
  scores <- methods::as(ShortRead::FastqQuality(Biostrings::quality(attr(dat, "quality"))), "matrix")

  return(dat[(rowSums(scores)/dim(scores)[2]) > min_score])
}
