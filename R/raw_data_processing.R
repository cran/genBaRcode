
#' @title Data processing
#'
#' @description Reads the corresponding fast(a/q) file(s), extracts the defined barcode constructs and counts them. Optionally,
#' a Phred-Score based quality filtering will be conducted and the results will be saved within a csv file.
#' @param file_name a character string or a character vector, containing the file name(s).
#' @param source_dir a character string which contains the path to the source files.
#' @param results_dir a character string which contains the path to the results directory.
#' @param mismatch an positive integer value, default is 0, if greater values are provided they indicate the number of allowed mismtaches when identifying the barcode constructes.
#' @param indels a logical value. If TRUE the chosen number of mismatches will be interpreted as edit distance and allow for insertions and deletions as well (currently under construction).
#' @param label a character string which serves as a label for every kind of created output file.
#' @param bc_backbone a character string describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param bc_backbone_label a character vector, an optional list of barcode backbone names serving as additional identifier within file names and BCdat labels. If not provided ordinary numbers will serve as alternative.
#' @param min_score a positive integer value, all fastq sequence with an average score smaller
#'  then min_score will be excluded, if min_score = 0 there will be no quality score filtering
#' @param min_reads positive integer value, all extracted barcode sequences with a read count smaller than min_reads will be excluded from the results
#' @param save_it a logical value. If TRUE, the raw data will be saved as a csv-file.
#' @param seqLogo a logical value. If TRUE, the sequence logo of the entire NGS file will be generated and saved.
#' @param cpus an integer value, indicating the number of available cpus.
#' @param strategy since the future package is used for parallelisation a strategy has to be stated, the default is "sequential"  (cpus = 1) and "multisession" (cpus > 1). For further information please read future::plan() R-Documentation.
#' @param full_output a logical value. If TRUE, additional output files will be generated.
#' @param bc_extraction a logical value. If TRUE, single reads will be stripped of the backbone and only the "wobble" positions will be left.
#' @param dist_method a character value. If "bc_backbone = 'none'", single reads will be clustered based on a distance measure.
#' Available distance methods are Optimal string aligment ("osa"), Levenshtein ("lv"), Damerau-Levenshtein ("dl"), Hamming ("hamming"), Longest common substring ("lcs"), q-gram ("qgram"), cosine ("cosine"), Jaccard ("jaccard"), Jaro-Winkler ("jw"),
#' distance based on soundex encoding ("soundex"). For more detailed information see stringdist function of the stringdist-package for more information)
#'
#' @return a BCdat object which includes reads, seqs, directories, BC_backbones.
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

processingRawData <- function(file_name, source_dir, results_dir, mismatch = 0, indels = FALSE, label = "", bc_backbone, bc_backbone_label = NULL, min_score = 30, min_reads = 2, save_it = TRUE, seqLogo = FALSE, cpus = 1, strategy = "sequential", full_output = FALSE, bc_extraction = TRUE, dist_method = "hamming") {

  source_dir <- .testDirIdentifier(source_dir)
  results_dir <- .testDirIdentifier(results_dir)

  if(length(bc_backbone_label) != length(bc_backbone)) {
      if(!is.null(bc_backbone_label)) {
        message(paste("#", length(bc_backbone), "BC pattern(s) but", length(bc_backbone_label),"BC pattern label(s)"))
      }
      bc_backbone_label = paste0("backbone_", 1:length(bc_backbone))
  }

  if(length(file_name) > 1) {
    dat <- processingRawData_multiple(file_name, source_dir, results_dir,
                                      mismatch, indels, label, bc_backbone,
                                      bc_backbone_label, min_score, min_reads,
                                      save_it, seqLogo, cpus, strategy, full_output, bc_extraction, dist_method)
  } else {
    if(length(file_name) == 1) {
      dat <- processingRawData_single(file_name, source_dir, results_dir,
                                      mismatch, indels, label, bc_backbone,
                                      bc_backbone_label, min_score, min_reads,
                                      save_it, seqLogo, full_output, cpus, strategy = "sequential", bc_extraction, dist_method)
    } else {
      stop("# file_name needs to be a character string (processingRawData)")
    }
  }

  return(dat)
}

processingRawData_single <- function(file_name, source_dir, results_dir, mismatch, indels, label, bc_backbone, bc_backbone_label, min_score, min_reads, save_it, seqLogo, full_output, cpus, strategy, bc_extraction, dist_method) {

  if(label == "") {
    label <- strsplit(file_name, split = "[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
  }

  ending <- strsplit(file_name, split="[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
  if(min_score > 0 & ending == "fastq" | ending == "FASTQ" | ending == "fq") {
    dat <- qualityFiltering(file_name, source_dir, min_score)
    if(length(dat) == 0) {
      warning("# There are no sequences after quality filtering!")
    }
  } else {
    if(ending == "fasta") {
      if(min_score > 0) {
        warning("# NGS quality filtering only available for FASTQ file formats! Proceeding without filtering!", call. = FALSE, immediate. = TRUE)
      }
      dat <- ShortRead::readFasta(dirPath = source_dir, pattern = file_name)
    } else {
      dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)
    }
  }

  if(seqLogo) {
    l<- nchar(as.character(ShortRead::sread(dat)[1]))
    ggplot2::ggsave(filename = paste0("seqLogo_", label, "_NGS.png"),
                    plot = suppressMessages(plotSeqLogo(as.character(ShortRead::sread(dat))) + ggplot2::scale_x_continuous(breaks = c(1, round(l/2), l))),
                    device = "png",
                    path = results_dir,
                    width = l / 2,
                    height = l / 16)
  }

    dat <- extractBarcodes(dat, label, results_dir,
                           mismatch = mismatch, indels = indels,
                           bc_backbone, full_output, cpus, strategy, bc_extraction, dist_method)

  if(length(bc_backbone) == 1) {
    dat <- prepareDatObject(dat[[1]], results_dir, label, bc_backbone, min_reads, save_it)
  } else {
    for(b in 1:length(bc_backbone)) {
      dat[[b]] <- prepareDatObject(dat[[b]], results_dir, label = paste(label, bc_backbone_label[b], sep = "_"),
                                   bc_backbone[b], min_reads, save_it)
    }
  }
  return(dat)
}

processingRawData_multiple <- function(file_name, source_dir, results_dir, mismatch, indels, label, bc_backbone, bc_backbone_label, min_score, min_reads , save_it, seqLogo, cpus, strategy, full_output, bc_extraction, dist_method) {

  cores <- future::availableCores()
  if(cpus > cores) {
    warning(paste("# available number of CPUs is", cores))
    cpus <- ifelse(cores > 1, cores - 1, 1)
  }

  if(cpus < 1) {
    warning("# cpus needs to be >= 1")
    cpus <- 1
    strategy <- "sequential"
  }

  if(length(label) != length(file_name)) {
    if(label != "") {
      message(paste0("# there are ", length(file_name), " file(s) but ", length(label), " label(s)"))
    }

    label <- unlist(lapply(strsplit(file_name, split = "[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE), function(x) { x[1] }))
  }

  if(cpus > length(file_name)) {
    cpus <- length(file_name)
  }

  if(strategy == "sequential" & cpus > 1) {
    strategy <- "multiprocess"
  }

  if(strategy != "sequential" & cpus == 1) {
    strategy <- "sequential"
  }

  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  suppressWarnings(future::plan(strategy = strategy, workers = cpus))

  tmp <- future.apply::future_lapply(1:length(file_name),
                                     function(i, source_dir, results_dir, mismatch, indels,
                                              label, bc_backbone, bc_backbone_label, min_score,
                                              min_reads, save_it, seqLogo, full_output, bc_extraction, dist_method) {

          processingRawData_single(file_name[i], source_dir, results_dir, mismatch,
                                   indels, label[i], bc_backbone, bc_backbone_label,
                                   min_score, min_reads, save_it, seqLogo, full_output, cpus = 1, strategy = "sequential", bc_extraction, dist_method)

  }, source_dir, results_dir, mismatch, indels,
     label, bc_backbone, bc_backbone_label, min_score,
     min_reads, save_it, seqLogo, full_output, bc_extraction, dist_method)

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
                        BC_backbone = bc_backbone))
  }

  dat <- dat[dat[, 2] > min_reads, ]
  if(dim(dat)[1] != 0) {
    dat <- data.frame(pos = 1:dim(dat)[1], read_count = as.numeric(unlist(dat[, 2])), barcode = as.character(unlist(dat[, 1])))
    if(save_it) {
      if(sum(dim(dat)) > 2) {
        utils::write.table(dat[, 2:3], paste(results_dir, label, "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
      } else {
        utils::write.table(dat, paste(results_dir, label, "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
    }

    return(methods::new(Class = "BCdat", reads = dat, results_dir = results_dir,
                        label = label,
                        BC_backbone = bc_backbone))
  } else {
    warning(paste("only barocdes with less than", min_reads, "reads detected for backbone:", bc_backbone))
    return(methods::new(Class = "BCdat", reads = data.frame(), results_dir = results_dir,
                        label = label,
                        BC_backbone = bc_backbone))
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
#' @param full_output a logical value. If TRUE additional output files will be generated in order to identify errors.
#' @param cpus an integer value, indicating the number of available cpus.
#' @param strategy since the future package is used for parallelisation a strategy has to be stated, the default is "sequential"  (cpus = 1) and "multiprocess" (cpus > 1). For further information please read future::plan() R-Documentation.
#' @param bc_extraction a logical value. If TRUE, single reads will be stripped of the backbone and only the "wobble" positions will be left.
#' @param dist_method a character value. If "bc_backbone = 'none'", single reads will be clustered based on a distance measure.
#' Available distance methods are Optimal string aligment ("osa"), Levenshtein ("lv"), Damerau-Levenshtein ("dl"), Hamming ("hamming"), Longest common substring ("lcs"), q-gram ("qgram"), cosine ("cosine"), Jaccard ("jaccard"), Jaro-Winkler ("jw"),
#' distance based on soundex encoding ("soundex"). For more detailed information see stringdist function of the stringdist-package for more information)
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

extractBarcodes <- function(dat, label, results_dir = "./", mismatch = 0, indels = FALSE, bc_backbone, full_output = FALSE, cpus = 1, strategy = "sequential", bc_extraction, dist_method = "hamming") {

  if(mismatch < 0) {
    warning("# mismatch needs to be an integer >= 0")
    mismatch <- 0
  }

  if(cpus < 1) {
    warning("# cpus needs to be >= 1")
    cpus <- 1
    strategy <- "sequential"
  }

  cores <- future::availableCores()
  if(cpus > cores) {
    warning(paste("# available number of CPUs is", cores))
    cpus <- ifelse(cores > 1, cores - 1, 1)
  }

  if(cpus > 1 && strategy == "sequential") {
    if(length(bc_backbone) == 1) {
      cpus <- 1
    } else {
      strategy <- "multiprocess"
    }
  }
  if(cpus == 1 && strategy == "multiprocess") {
    strategy <- "sequential"
  }

  if(length(bc_backbone) == 1) {
    if(bc_backbone == "none") {
      tmp <- as.character(ShortRead::sread(dat))
      tmp <- data.frame(table(tmp))
      tmp <- asBCdat(dat = tmp[, c(2, 1)], label = label, BC_backbone = "none", resDir = results_dir)

      res <- list()
      if(mismatch > 0) {
        tmp <- errorCorrection_single_variation(BC_dat = tmp, maxDist = mismatch, save_it = FALSE, m = dist_method, EC_analysis = FALSE, only_EC_BCs = FALSE)
      }
      res[[1]] <- methods::slot(tmp, "reads")[, 3:2]

      return(res)
    }
  } else {
    if(sum(bc_backbone == "none") > 0) {
      stop("# At the moment it is only possible to analyse the data either with no backbone at all or with one or multiple actual backbones.")
    }
  }

  tmp <- list(list(), ShortRead::sread(dat))

  for(i in 1:length(bc_backbone)) {
    tmp <- BC_matching(tmp, label[i], results_dir, mismatch, indels, bc_backbone[i], full_output)
  }

  if(bc_extraction) {
  # future limit == 500MB
    if(utils::object.size(tmp) <= 499000000) {

      oplan <- future::plan()
      on.exit(future::plan(oplan), add = TRUE)
      suppressWarnings(future::plan(strategy = strategy, workers = cpus))

      res <- future.apply::future_lapply(1:length(bc_backbone), extractWobble, tmp = tmp, bc_backbone = bc_backbone, indels, label)
    } else {
      res <- lapply(1:length(bc_backbone), extractWobble, tmp = tmp, bc_backbone = bc_backbone, indels, label)
    }
  } else {
    res <- vector(mode = "list", length = length(tmp[[1]]))
    for(i in 1:length(tmp[[1]])) {
      res_tmp <- sort(table(tmp[[1]][[i]]), decreasing = TRUE)
      res[[i]] <- data.frame(seqs = names(res_tmp), n = as.numeric(res_tmp))
    }
  }
  return(res)
}

BC_matching <- function(dat, label, results_dir, mismatch, indels, bc_backbone, full_output) {

  seqs <- dat[[2]]
  read_length <- ShortRead::width(dat[[2]])[1]

  # first, w/o mismatch, because faster
  match_index <- Biostrings::vmatchPattern(bc_backbone, seqs,
                                           max.mismatch = 0, min.mismatch = 0,
                                           with.indels = indels, fixed = FALSE,
                                           algorithm = "auto")
  match_matrix <- as.data.frame(match_index)

  if(sum(duplicated(match_matrix$group)) != 0) {
    warning("# Duplicated backbone matches (mismatch == 0)")
    if(full_output) {
      group <- unique(match_matrix$group[duplicated(match_matrix$group)])
      utils::write.table(match_matrix[match_matrix$group %in% group, ], paste(results_dir, label, "_", bc_backbone, "_duplicatedMatcher_m0.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    match_matrix <- match_matrix[!duplicated(match_matrix$group), ]
  }

  BCs <- substr(as.character(seqs)[match_matrix$group],
                start = match_matrix$start,
                stop = match_matrix$end)
  if(!indels & length(table(nchar(BCs))) > 1) {
    BCs <- fixBorderlineBCs(BCs, mismatch, match_matrix, read_length)
  }
  seqs <- seqs[-match_matrix$group]

  # second, w/ mismatch but smnaller data set
  if(mismatch > 0) {
    match_index <- Biostrings::vmatchPattern(bc_backbone, seqs,
                                             max.mismatch = mismatch, min.mismatch = 1,
                                             with.indels = indels, fixed = FALSE,
                                             algorithm = "auto")
    match_matrix <- data.frame(match_index)
    if(sum(duplicated(match_matrix$group)) != 0) {
      warning("# Duplicated backbone matches (mismatch > 0)")
      if(full_output) {
        group <- unique(match_matrix$group[duplicated(match_matrix$group)])
        utils::write.table(match_matrix[match_matrix$group %in% group, ], paste(results_dir, label, "_", bc_backbone, "_duplicatedMatcher_m, ", mismatch, ", .csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
      match_matrix <- match_matrix[!duplicated(match_matrix$group), ]
    }

    BCsTmp <- substr(as.character(seqs)[match_matrix$group],
                     start = match_matrix$start,
                     stop = match_matrix$end)
    seqs <- seqs[-match_matrix$group]

    if(!indels & length(table(nchar(BCsTmp))) > 1) {
      BCsTmp <- fixBorderlineBCs(BCsTmp, mismatch, match_matrix, read_length)
    }
    BCs <- c(BCs, BCsTmp)
  }

  if(!indels & length(table(nchar(BCs))) > 1) {
    warning("# different barcode lengths (function: BCmatching)")
    if(full_output) {
      print(table(nchar((as.character(unlist(dat))))))
      utils::write.table(BCs, paste(results_dir, label, "_", bc_backbone, "_diffLengthBCs.csv", sep=""), sep=";", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }

  dat[[1]][[length(dat[[1]]) + 1]] <- BCs
  dat[[2]] <- seqs
  return(dat)
}

fixBorderlineBCs <- function(BCs, mismatch, match_matrix, read_length) {

  for(m in 1:mismatch) {
    if(sum(match_matrix$start == (1-m)) > 0) {
      BCs[match_matrix$start == (1-m)] <- paste0(paste(rep("N", m), collapse = ""), BCs[match_matrix$start == (1-m)])
    } else {
      if(sum(match_matrix$end == read_length+m) > 0) {
        BCs[match_matrix$start == read_length] <- paste0(BCs[match_matrix$start == 0], paste(rep("N", m), collapse = ""))
      } else {
        warning("# different barcode lengths (function: fixBorderBCs)")
      }
    }
  }
  return(BCs)
}

extractWobble <- function(i, tmp, bc_backbone, indels, label) {

  BCs <- tmp[[1]][[i]]
  if(length(BCs) != 0) {
    if(!indels) {
      wobble_pos <- .getWobblePos(bc_backbone[i])

      BCs <- strsplit(BCs, split = "")
      tmp <- lapply(BCs, function(x, wobble_pos) {
        paste0(x[wobble_pos], collapse = "")
      }, wobble_pos
      )
      dat <- data.frame(unlist(tmp))

    } else {
      dat <- as.data.frame(dat)
    }
    return(dplyr::count(dat, dat[, 1], sort = TRUE))
  } else {
    #warning(paste("# ", label, ": no backbone patterns detectable!", sep =""))
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
# @examples
# data(BC_dat)
# hybridsIdentification(BC_dat, min_seq_length = 2)

hybridsIdentification <- function(dat, min_seq_length = 10) {

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
          if(sum(dat[, 2] == hybrid) > 1) { stop(" # found same BC twice (find.hybrid.after.EC) ## ") }
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

qualityFiltering <- function(file_name, source_dir, min_score = 30) {

  dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)
  scores <- methods::as(ShortRead::FastqQuality(Biostrings::quality(attr(dat, "quality"))), "matrix")

  return(dat[(rowSums(scores)/dim(scores)[2]) > min_score])
}
