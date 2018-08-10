
#' @title Shiny App
#' @description Launches the corresponding shiny app.
#'
#' @param dat_dir a character string, identifying the path to one or more fast(q) files which shall be analysed, default is the
#' path to the package inherent example fastq file
#'
#' @export

genBaRcode_app <- function(dat_dir = system.file("extdata", package = "genBaRcode")) {

  tmp <- unlist(strsplit(dat_dir, split = ""))
  if(tmp[length(tmp)] != .Platform$file.sep) {
    dat_dir <- paste0(dat_dir, .Platform$file.sep)
  }

  options("genBaRcode-shinyDir" = dat_dir)
  options("genBaRcode-info" = "")

  appDir <- system.file("shiny_app", package = "genBaRcode")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing the `genBaRcode` package.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", quiet = TRUE)

}

#' Predefined Barcode Backbone Sequences
#'
#' allows one to choose between predefined backbone sequences. Excecution of the function without any parameter value will display all
#' available backbone sequences. The id parameter will accept the name of the backbone or the rownumber of the shown selection.
#'
#' @param id an integer or character value in order to choose a specific backbone.
#'
#' @return a character string.
#' @export

getBackbone <- function(id = NULL) {

  BB_names <- c("BC32-GFP",
                "BC32-Venus",
                "BC32-eBFP",
                "BC32-T-Sapphire",
                "BC16-GFP",
                "BC16-Venus",
                "BC16-mCherry",
                "BC16-Cerulean")
  BB_seqs <- c("ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN",
               "CGANNAGANNCTTNNCGANNCTANNGGANNCTTNNCGANNAGANNCTTNNCGANNCTANNGGANNCTTNNCGANNAGANN",
               "CTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNN",
               "CAGNNATCNNCTTNNCGANNGGANNCTANNCTTNNCAGNNATCNNCTTNNCGANNGGANNCTANNCTTNNCAGNNATCNN",
               "ATCNNTAGNNTCCNNAAGNNTCGNNAAGNNTCGNNAGTNNTAG",
               "CTANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNGAT",
               "CTANNCAGNNATCNNCTTNNCGANNGGANNCTANNCTTNNGAT",
               "CTANNCACNNAGANNCTTNNCGANNCTANNGGANNCTTNNGAT")

  if(is.null(id)) {
    return(print(format(data.frame(name = BB_names, sequences = BB_seqs), justify = "left", width = 20), right = FALSE))
  } else {
    if(is.numeric(id)) {
      BB_seqs[id]
    } else {
      if(is.character(id))
        if(sum(BB_names == id) == 0) {
            message("# No backbone with such name known.")
            return()
        } else {
          return(BB_seqs[BB_names == id])
        }
    }
  }
}

#' @importFrom foreach %dopar%
makeParallel <- function(dat, func, cpus, split = FALSE, ...) {

  if(split) {
    psize <- floor(length(dat) / cpus)
    dat <- lapply(1:cpus, function(x) {
      if(x == 1) {
        dat[1:psize]
      } else {
        if(x == cpus) {
          dat[((psize*(x-1))+1):length(dat)]
        } else {
          dat[((psize*(x-1))+1):(psize*x)]
        }
      }
    })
  }

  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)

  # due to dopar problems
  i <- NULL
  tmp <- foreach::foreach(i = 1:length(dat)) %dopar% {
    func(dat[[i]], ...)
  }

  parallel::stopCluster(cl)

  return(tmp)
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
#' @importFrom methods new

asBCdat <- function(dat, label = "without_label", mask = "", resDir = getwd()) {

  if(dim(dat)[2] < 2 | dim(dat)[2] > 3) {
    stop("# Data obect needs two or three columns!")
  }

  if(dim(dat)[2] == 2) {
    dat <- data.frame(pos = 1:dim(dat)[1], read_count = as.numeric(dat[, 1]), barcode = as.character(dat[, 2]))
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
#' @importFrom methods new

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
#' @param bc_backbone a character vector.

.getWobblePos <- function(bc_backbone = "") {

  bc_backbone <- unlist(strsplit(bc_backbone, ""))
  wobble_pos <- which(bc_backbone == "N")

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
#' @param wobble_pos a character string.

.getBarcodeFilter <- function(wobble_pos) {

  seq_length <- NULL
  counter <- 1
  for(i in 1:(length(wobble_pos)-1)) {
    if(wobble_pos[i+1]-wobble_pos[i] == 1) {
              counter <- counter +1
    } else {
              seq_length <- c(seq_length, counter)
              counter <- 1
    }
  }
  seq_length <- c(seq_length, counter)

  awk_cmd <- "| awk '{print "
  wobble_pos_tmp <- wobble_pos
  for(i in 1:length(seq_length)) {
    awk_cmd <- paste(awk_cmd, "substr($0, ", wobble_pos_tmp[1], ", ", seq_length[i], ")", sep="")
    if(i < length(seq_length)) {
        wobble_pos_tmp <- wobble_pos_tmp[-(1:seq_length[i])]
    }
  }
  BC_filter <- paste(awk_cmd, "}'", sep="")

  return(BC_filter)
}

#' @title Internal function
#'
#' @description Creates a search file for a command line grep search.
#'
#' @param bc_backbone a character string (barcode pattern).
#' @param patterns_file a character string (file name)

.createPatternFile <- function(bc_backbone, patterns_file) {

  bc_backbone <- unlist(strsplit(bc_backbone, ""))
  nucs <- c("A", "T", "C" ,"G")
  counter <- 2
  results <- list()
  results[[1]] <- bc_backbone

  for(i in 1:length(bc_backbone)) {
    if(bc_backbone[i] != ".") {
      tmp <- nucs[-which(nucs == bc_backbone[i])]
      for(k in 1:3) {
        results[[counter]] <- bc_backbone
        results[[counter]][i] <- tmp[k]
        counter <- counter + 1
      }
    }
  }
  utils::write.table(lapply(results, function(x) { paste(x, collapse = "") }), patterns_file, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}

#' @title Distance calculation
#' @description Calculates the minimum distance to a set of predefined barcodes for a given list of barcode.
#' @param BC_dat a BCdat object
#' @param ori_BCs a character vector containing barcodes to which the minimal hamming distance will be calculated.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).

.getMinDist <- function(BC_dat, ori_BCs, m = "hamming") {

  HD_values <- unlist(lapply(methods::slot(BC_dat, "reads")$"barcode", function(x) {
        min(stringdist::stringdist(ori_BCs, x, m))
  }))

  return(HD_values)
}

#' @title Color list generation
#'
#' @description Generates a collection of colors for a list of barcodes based on their identified minimum hamming distances.
#'
#' @param minHD a numeric vector of all the minimum hamming distances.
#' @param type a character string. Possible Values are "rainbow", "heat.colors", "topo.colors" (see package "grDevices").
#' @param alpha a numeric value between 0 and 1, modifies colour transparency.

.generateColors <- function(minHD, type = "rainbow", alpha = 1) {

  BC_dist <- sort(unique(minHD))

  if(type == "rainbow") {
    color_list <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"), alpha = alpha)(length(unique(minHD)))
  }
  if(type == "heat") {
    color_list <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), alpha = alpha)(length(unique(minHD)))
  }
  if(type == "topo.colors") {
    color_list <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"), alpha = alpha)(length(unique(minHD)))
  }
  if(type == "greens") {
    color_list <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Greens"), alpha = alpha)(length(unique(minHD)))
  }
  if(type == "wild") {
    color_list <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"), alpha = alpha)(length(unique(minHD)))
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


