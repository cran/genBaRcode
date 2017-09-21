
#' @title Shiny App
#' @description Launches the corresponding shiny app.
#'
#' @param dat_dir a character string, identifying the path to one or more fast(q) files which shall be analysed, default is the
#' path to the package inherent example fastq file
#'
#' @export

genBaRcode_app <- function(dat_dir = system.file("extdata", package = "genBaRcode")) {

  options("genBaRcode-shinyDir" = paste0(dat_dir, .Platform$file.sep))
  options("genBaRcode-info" = "")

  appDir <- system.file("shiny_app", package = "genBaRcode")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing the `genBaRcode` package.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")

}

#' @title Data Type Conversion
#'
#' @description Converts a data.frame into a BCdat object.
#'
#' @param dat a data.frame object with two columns containing read counts and barcode sequences.
#' @param label a optional character string used as label.
#' @param mask a optional character string, describing the barcode backbone structure.
#' @param resDir a optional character string, identifying the path to the results directory, default is current working directory.
#'
#' @return a BCdat object.
#' @export

as.BCdat <- function(dat, label = "without_label", mask = "", resDir = getwd()) {

  if(dim(dat) != 2) {
    stop("# Data obect needs two columns!")
  }

  return(methods::new(Class = "BCdat", reads = dat, results_dir = resDir,
                      label = label,
                      mask = mask))
}

#' Converts hex colors into gephi usable rgb colors
#'
#' @param colrs a character vector containing a list of hex colors
#'
#' @return
#'

.hex2rgbColor <- function(colrs) {

  if(is.vector(colrs)) {
    unlist(lapply(colrs, function(col) {
      tmp <- c(strtoi(substr(col, 2, 3), 16),
               strtoi(substr(col, 4, 5), 16),
               strtoi(substr(col, 6, 7), 16))
      paste("'", paste(tmp, collapse = ","), "'", sep = "")
    }))
  } else {
    warning("# hex2rgbColor needs to get a vector of character values")
    return(NULL)
  }
}

#' @title Data Input
#'
#' @description Reads in a data table and returns a BCdat objects.
#'
#' @param path a character string containing the path to a saved read count table (two columns containing read counts and barcode
#' sequences).
#' @param label a character string containing a label of the data set.
#' @param mask a character string containing the barcode structure information.
#' @param file_name a character string containing the name of the file to read in.
#' @param s a character value, identifying the column separating char.
#'
#' @return a BCdat object.
#' @export

readBCdat <- function(path = "./", label = "", mask = "", file_name, s = ";") {

  dat <- utils::read.table(paste(path, file_name, sep = ""), header = TRUE, sep = s)
  dat <- data.frame(pos = 1:dim(dat)[1], read_count = as.numeric(dat$read_count), barcode = as.character(dat$barcode))

  if(label == "") {
    label <- unlist(strsplit(file_name, split = "[.]"))[1]
  }
  if(path == "./") {
    path <- getwd()
  }

  return(methods::new(Class = "BCdat", reads = dat, results_dir = path,
               label = label,
               mask = mask))

}

#' @title Index Generation
#'
#' @description Generates a matrix index to create a square triangular matrix.
#'
#' @param n an integer indicating the size of the resulting index matrix.
#'
#' @return a locigal matrix of size \code{n} x \code{n}

.getDiagonalIndex <- function(n) {

  index <- matrix(TRUE, nrow = n, ncol = n)
  return(
    matrix(
        unlist(
          lapply(1:(n), function(i) {
                                            c(rep(TRUE, i-1), rep(FALSE, length(i:n)))
                                    })
          )
        , ncol = n, byrow=TRUE
      )
    )
}

#' @title getWobblePos
#'
#' @description Extracts barcode positions.
#'
#' @param bc_pattern a character vector.

.getWobblePos <- function(bc_pattern = "") {

  bc_pattern <- unlist(strsplit(bc_pattern, ""))
  wobble_pos <- which(bc_pattern == "N")

  return(wobble_pos)
}

#' @title Internal function
#'
#' @description Checks directory paths for correctness and if nessesary corrects them.
#'
#' @param s a character string.

.testDirIdentifier <- function(s) {

  dir_sep <- ifelse(.Platform$OS.type == "unix", "/", "\\")
  tmp <- unlist(strsplit(s, split = ""))
  if(tmp[length(tmp)] != dir_sep) {
    s <- paste(s, dir_sep, collapse = "", sep = "")
  }
  return(s)

}

#' @title Internal function
#'
#' @description Identifies the barcode positions within the mask and generates a awk command.
#'
#' @param wooble_pos a character string.

.getBarcodeFilter <- function(wooble_pos) {

  seq_length <- NULL
  counter <- 1
  for(i in 1:(length(wooble_pos)-1)) {
    if(wooble_pos[i+1]-wooble_pos[i] == 1) {
              counter <- counter +1
    } else {
              seq_length <- c(seq_length, counter)
              counter <- 1
    }
  }
  seq_length <- c(seq_length, counter)

  awk_cmd <- "| awk '{print "
  wooble_pos_tmp <- wooble_pos
  for(i in 1:length(seq_length)) {
    awk_cmd <- paste(awk_cmd, "substr($0, ", wooble_pos_tmp[1], ", ", seq_length[i], ")", sep="")
    if(i < length(seq_length)) {
        wooble_pos_tmp <- wooble_pos_tmp[-(1:seq_length[i])]
    }
  }
  BC_filter <- paste(awk_cmd, "}'", sep="")

  return(BC_filter)
}

#' @title Internal function
#'
#' @description Creates a search file for a command line grep search.
#'
#' @param bc_pattern a character string (barcode pattern).
#' @param patterns_file a character string (file name)

.createPatternFile <- function(bc_pattern, patterns_file) {

  bc_pattern <- unlist(strsplit(bc_pattern, ""))
  nucs <- c("A", "T", "C" ,"G")
  counter <- 2
  results <- list()
  results[[1]] <- bc_pattern

  for(i in 1:length(bc_pattern)) {
    if(bc_pattern[i] != ".") {
      tmp <- nucs[-which(nucs == bc_pattern[i])]
      for(k in 1:3) {
        results[[counter]] <- bc_pattern
        results[[counter]][i] <- tmp[k]
        counter <- counter + 1
      }
    }
  }
  utils::write.table(lapply(results, function(x) { paste(x, collapse = "") }), patterns_file, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}

#' @title Hamming distance calculation
#' @description Calculates the minimum hamming distance to a set of predefined barcodes for a given list of barcode.
#' @param BC_dat a BCdat object
#' @param ori_BCs a character vector containing barcodes to which the minimal hamming distance will be calculated.

.getMinHammDist <- function(BC_dat, ori_BCs) {

  HD_values <- unlist(lapply(methods::slot(BC_dat, "reads")$"barcode", function(x) {
        min(stringdist::stringdist(ori_BCs, x, "h"))
  }))

  return(HD_values)
}

#' @title Color list generation
#'
#' @description Generates a collection of colors for a list of barcodes based on their identified minimum hamming distances.
#'
#' @param minHD a numeric vector of all the minimum hamming distances
#' @param type a character string. Possible Values are "rainbow", "heat.colors", "topo.colors" (see package "grDevices")

.generateColors <- function(minHD, type = "rainbow") {

  BC_dist <- sort(unique(minHD))

  if(type == "rainbow") {
    color_list <- grDevices::rainbow(length(unique(minHD)))
  }
  if(type == "heat") {
    color_list <- grDevices::heat.colors(length(unique(minHD)))
  }
  if(type == "topo.colors") {
    color_list <- grDevices::topo.colors(length(unique(minHD)))
  }

  dist_col <- unlist(lapply(minHD, function(x) {
              color_list[BC_dist == x]
          }))

  names(color_list) <- BC_dist

  return(list(dist_col, color_list))
}

#' @title DNA string manipulation
#'
#' @description Converts a vector of character strings (DNA sequences) into its reverse complement.
#'
#' @param seq_dat a character vector containing DNA sequences

.revComp <- function(seq_dat) {

  # if(!is.character(seq_dat)) {
  #   stop("# a character vector is required")
  # }

  tmp <- unique(unlist(strsplit(as.character(seq_dat), split = "")))
  if(sum(tmp %in% c("A", "C", "T", "G", "N")) != length(tmp)) {
    stop("# only the letter 'A', 'T', 'C', 'G' and 'N' are allow!")
  }

  word_length <- unique(nchar(as.character(seq_dat)))
  if(length(word_length) == 1) {
      .revComp_EqLength(seq_dat, word_length)
  } else {
      .revComp_UneqLength(seq_dat)
  }

}

#' @title DNA string manipulation for equal string sizes
#' @description Converts a vector of equally long character strings into its reverse complement.
#'
#' @param seq_dat a character vector.
#' @param word_length an integer giving the word length.

.revComp_EqLength <- function(seq_dat, word_length) {

  seq_dat <- unlist(strsplit(seq_dat, ""))
  seq_dat[seq_dat == "A"] <- "a"
  seq_dat[seq_dat == "T"] <- "t"
  seq_dat[seq_dat == "G"] <- "g"
  seq_dat[seq_dat == "C"] <- "c"

  seq_dat[seq_dat == "a"] <- "T"
  seq_dat[seq_dat == "t"] <- "A"
  seq_dat[seq_dat == "g"] <- "C"
  seq_dat[seq_dat == "c"] <- "G"

  return(apply(matrix(seq_dat, ncol = word_length, byrow = TRUE), 1, function(x) paste(x[length(x):1], collapse = "")))

}

#' @title DNA string manipulation for unequal string sizes
#' @description Converts a vector of unequally long character strings into the reverse complement.
#'
#' @param seq_dat A character vector.

.revComp_UneqLength <- function(seq_dat) {

  seq_dat <- lapply(seq_dat, function(x) {
    x <- unlist(strsplit(x, ""))
    x[x == "A"] <- "a"
    x[x == "T"] <- "t"
    x[x == "G"] <- "g"
    x[x == "C"] <- "c"

    x[x == "a"] <- "T"
    x[x == "t"] <- "A"
    x[x == "g"] <- "C"
    x[x == "c"] <- "G"
    paste(x[length(x):1], collapse="")
  })

  return(seq_dat)
}


