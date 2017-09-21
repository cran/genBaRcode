
#' @title Data processing
#'
#' @description Reads the corresponding fast(q) file(s), extracts the defined barcode constructs and counts them. Optionally,
#' a Phred-Score based quality filtering will be conducted and the results will be saved within a csv file.
#' @param file_name a character string which will serve as file name.
#' @param source_dir a character string which contains the path to the source files.
#' @param results_dir a character string which contains the path to the results directory.
#' @param mismatch an positive integer value, default is 0, if greater values are provided they indicate the number of allowed mismtaches when identifying the barcode constructes.
#' @param label a character string which serves as a label for every kind of created output file.
#' @param bc_pattern a character string describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param quality_filtering a logical value. If TRUE a quality filtering will be applied before extracting the barcode sequences
#' @param min_score a positive integer value, only relevant if quality_filtering is TRUE, all fastq sequence with an average score smaller
#'  then min_score will be excluded
#' @param min_reads positive integer value, all extracted barcode sequences with a read count smaller than min_reads will be excluded from the results
#' @param unix under construction
#' @param save_it a logical value. If TRUE, the raw data will be saved as a csv-file.
#' @param cpus a positive integer identifying the number of usable CPUs.
#'
#' @return a BCdat object which includes reads, seqs, directories, masks.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' bc_pattern <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
#'
#' source_dir <- system.file("extdata", package = "genBaRcode")
#'
#' processingRawData(file_name = "test_data.fastq", source_dir, results_dir = getwd(), mismatch = 0,
#' label = "test", bc_pattern, quality_filtering = FALSE, min_score = 30,
#' min_reads = 2, unix = FALSE, save_it = TRUE)
#'
#' }

processingRawData <- function(file_name, source_dir, results_dir, mismatch = 0, label = "", bc_pattern, quality_filtering = FALSE, min_score = 30, min_reads = 2, unix = FALSE, save_it = TRUE, cpus = 1) {

  if(.Platform$OS.type != "unix" & unix) {
    warning("# no unix based system recognized, set parameter 'unix' to FALSE")
    unix <- FALSE
  }

  source_dir <- .testDirIdentifier(source_dir)
  results_dir <- .testDirIdentifier(results_dir)

  if(label == "") {
    label <- strsplit(file_name, split = "[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
  }

  ending <- strsplit(file_name, split="[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
  if(quality_filtering) {
    if(ending == "fastq" | ending == "FASTQ" | ending == "fq") {
      dat <- qualityFiltering(file_name, source_dir, min_score)
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

  if(unix) {
      stop("under construction...")
      # extractSeqs_unix(file_name, label, results_dir, source_dir, mismatch, bc_pattern)
      # dat <- utils::read.table(paste(results_dir, label, "_matching_seqs.tmp", sep=""))
      # dat <- dat[dat[, 1] > min_reads, ]
      # dat <- data.frame(pos = 1:dim(dat)[1], read_count = dat[order(dat[, 1], decreasing = TRUE), 1], barcode = dat[order(dat[, 1], decreasing = TRUE), 2])
  } else {
      dat <- extractBarcodes(dat, label, results_dir, mismatch = mismatch, indels = FALSE, bc_pattern, cpus)
      if(length(dat) != 0) {
        dat <- dat[order(dat, decreasing = TRUE)]
        dat <- dat[dat > min_reads]
        if(length(dat) != 0) {
          dat <- data.frame(pos = 1:length(dat), read_count = as.numeric(dat), barcode = names(dat))
        }
      }
  }

  if(save_it) {
    if(sum(dim(dat)) > 2) {
      utils::write.table(dat[, 2:3], paste(results_dir, label, "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
    } else {
      utils::write.table(dat, paste(results_dir, label, "_BCs.csv", sep=""), sep=";", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  }

  if(!is.null(dim(dat))) {
    dat <- methods::new(Class = "BCdat", reads = dat, results_dir = results_dir,
                        label = label,
                        mask = bc_pattern)
  } else {
    warning(paste("# ", label, ": No barcode sequences saved!", sep = ""))
    dat <- methods::new(Class = "BCdat", reads = as.data.frame(NULL), results_dir = results_dir,
                        label = label,
                        mask = bc_pattern)
  }

  return(dat)
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
#' @param bc_pattern a character string describing the barcode design, variable positions have to be marked with the letter 'N'.
#' @param cpus a positive integer identifying the number of usable CPUs.
#'
#' @return a frequency table of barcode sequences.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' bc_pattern <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
#' source_dir <- system.file("extdata", package = "genBaRcode")
#' dat <- ShortRead::readFastq(dirPath = source_dir, pattern = "test_data.fastq")
#'
#' extractBarcodes(dat, label = "test", results_dir = getwd(), mismatch = 0,
#' indels = FALSE, bc_pattern)
#'
#' }

extractBarcodes <- function(dat, label, results_dir, mismatch = 0, indels = FALSE, bc_pattern, cpus = 1) {

  if(mismatch < 0) {
    warning("# mismatch needs to be an integer >= 0")
    mismatch <- 0
  }

  if(cpus < 1) {
    warning("# cpus needs to be >= 1")
    cpus <- 1
  }

  if(cpus > getOption("mc.cores", 2L)) {
    warning(paste("# available number of CPUs is", getOption("mc.cores", 2L)))
    cpus <- getOption("mc.cores", 2L)
  }

  if(indels) {

    match_index <- unlist(lapply(ShortRead::sread(dat), function(x) {
                                Biostrings::matchPattern(bc_pattern, x,
                                      max.mismatch = mismatch, min.mismatch = 0,
                                      with.indels = TRUE, fixed = FALSE,
                                      algorithm = "auto")
                    }))
  } else {

    wobble_pos <- .getWobblePos(bc_pattern)
    match_index <- Biostrings::vmatchPattern(bc_pattern, ShortRead::sread(dat),
                            max.mismatch = mismatch, min.mismatch = 0,
                            with.indels = FALSE, fixed = FALSE,
                            algorithm = "auto")
  }

  dat <- ShortRead::sread(dat)[match_index]
  if(sum(ShortRead::width(dat) != 0) != 0) {
      dat <- dat[ShortRead::width(dat) != 0]

      if(indels) {
              return(table(as.character(dat)))
      } else {
              return(table(unlist(parallel::mclapply(as.character(dat), function(x, wobble_pos) {
                paste(unlist(strsplit(x, split= ""))[wobble_pos], collapse = "")
              }, wobble_pos, mc.cores = cpus))))
      }
  } else {
      warning(paste("# ", label, ": no backbone patterns detectable!", sep =""))
      return(NULL)
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
#' source_dir <- system.file("extdata", package = "genBaRcode")
#' qualityFiltering(file_name = "test_data.fastq", source_dir, results_dir = getwd(), min_score = 30)

qualityFiltering <- function(file_name, source_dir, results_dir, min_score = 30) {

  dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)
  scores <- methods::as(ShortRead::FastqQuality(Biostrings::quality(attr(dat, "quality"))), "matrix")

  return(dat[(rowSums(scores)/dim(scores)[2]) > min_score])
}
