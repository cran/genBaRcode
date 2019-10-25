
checkDir <- function(dir, errors = NULL) {

  if (!dir.exists(dir)) {
    errors <- c(errors, "# Unfortunatly, the specified folder does not exists.")
  } else {
    if (length(dir) > 1) {
      errors <- c(errors, "# Only one results path required.")
    }
  }

  return(errors)
}

checkReads <- function(reads, errors = NULL) {

  if (!is.data.frame(reads)) {
    errors <- c(errors, paste0("# Data is of type ", methods::is(reads)[1], ".  Should be a data.frame."))
  }

  if (dim(reads)[2] > 2) {
    errors <- c(errors, paste0("# Data consists of ", dim(reads)[2], " column(s).  Should be 2."))
  } else {
    if (dim(reads)[2] == 2) {
      if (sum(names(reads) == c("read_count", "barcode")) != 2) {
        errors <- c(errors, "# Data needs to consist of the columns 'read_count' and 'barcode', in that order.")
      }
    }
  }
  return(errors)
}

#' @importFrom dplyr %>%
checkBackbone <- function(backbone, errors = NULL) {

  if (backbone != "not defined") {
    if (backbone != "none") {
      if (!is.character(backbone)) {
        errors <- c(errors, paste0("# The barcode backbone is of type ", methods::is(backbone)[1], ".  Should be a character string."))
      }
      if (length(backbone) > 1) {
        errors <- c(errors, "# Only one backbone is supported.")
      } else {
        elements <- strsplit(backbone, split = "") %>% unlist %>% table %>% as.data.frame()
        IUPAC_nucCode <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", ".", "-")
        if (sum(!as.character(elements[, 1]) %in% IUPAC_nucCode) != 0) {
          errors <- c(errors, "# Backbones are only valid if consisting of IUPAC-nucleotide-code symbols")
        }
      }
    }
  }
  return(errors)
}

checkLabel <- function(label, errors = NULL) {

  if (length(label) != 1) {
    errors <- c(errors, "# Only one label needed.")
  }
  if (!is.character(label)) {
    errors <- c(errors, paste0("# Label is of type ", methods::is(label)[1], ".  Should be a character string."))
  }
  return(errors)
}

checkBarcodeData <- function(object) {

  errors <- character()
  errors <- checkReads(object@reads, errors)
  errors <- checkBackbone(object@BC_backbone, errors)
  errors <- checkDir(object@results_dir, errors)
  errors <- checkLabel(object@label, errors)

  if (length(errors) == 0) TRUE else errors
}

#' BCdat class.
#'
#' @slot reads data.frame containing barcode sequences and their corresponding read counts.
#' @slot results_dir character string of the working directory path.
#' @slot label character string identifying the particular experiment (will be part of the names of any file created).
#' @slot BC_backbone character string of the used barcode design (also called barcode backbone).
#'
#' @exportClass BCdat
BCdat <- methods::setClass("BCdat",
                            slots = list(
                                       reads = "data.frame",
                                       results_dir = "character",
                                       label = "character",
                                       BC_backbone = "character"
                            ),
                            prototype = list(
                              reads = data.frame(read_count = NULL, barcode = NULL),
                              results_dir = NA_character_,
                              label = NA_character_,
                              BC_backbone = NA_character_
                            ),
                            validity = checkBarcodeData
)


#' @importFrom dplyr %>%
setMethod("show", signature = c("BCdat"),
          definition = function(object){

            len <- 45

            cat(" class: BCdat\n\n")
            cat("  number of barcode sequences:", dim(object@reads)[1], "\n")
            if (sum(dim(object@reads)) > 1) {
              cat("  read count distribution: min", min(object@reads$read_count),
                  " mean", round(mean(object@reads$read_count), digits = 2),
                  " median", stats::median(object@reads$read_count),
                  " max", max(object@reads$read_count), "\n")
              l <- sort(unique(nchar(as.character(object@reads$barcode))))
              if (length(l) == 1) {
                cat("  barcode sequence length:", l[1], "\n")
              }
              if (length(l) == 2) {
                cat("  barcode sequence lengths: ", l[1], ", ", l[2], "\n")
              }
              if (length(l) == 3) {
                cat("  barcode sequence lengths: ", l[1], ", ", l[2], ", ", l[3], "\n")
              }
              if (length(l) > 3) {
                cat("  barcode sequence lengths: ", l[1], ", ", l[2], ", ", l[3], ", ...", l[length(l)], "\n")
              }
              cat("\n")
              cat(" barcode read counts:\n")
              if (dim(object@reads)[1] > 10) {
                if (sum(l > (len + 3)) > 0) { # >35 because of three additional letters ("...") inbetween the particular word
                  seqs <- lapply(as.character(object@reads[1:10, 2]), function(x) {
                    ind <- x %>% nchar %>% -len
                    tmp_seq <- (strsplit(x, split = "") %>% unlist)[c(1:round(len/2), (round(len/2) + ind):nchar(x))]
                    paste0(paste0(tmp_seq[1:round(len/2)], collapse = ""),
                           "...",
                           paste0(tmp_seq[(round(len/2)+1):length(tmp_seq)], collapse = ""),
                           collapse = "")
                  }) %>% unlist
                  print(data.frame(read_count = object@reads[1:10, 1],
                                   barcode = seqs), row.names = FALSE)
                } else {
                  print(object@reads[1:10, ], row.names = FALSE)
                }
                cat("                       ...")
              } else {
                if (sum(l > (len + 3)) > 0) { # >35 because of three additional letters ("...") inbetween the particular word
                  seqs <- lapply(as.character(object@reads[, 2]), function(x) {
                    ind <- x %>% nchar %>% -len
                    tmp_seq <- (strsplit(x, split = "") %>% unlist)[c(1:round(len/2), (round(len/2) + ind):nchar(x))]
                    paste0(paste0(tmp_seq[1:round(len/2)], collapse = ""),
                           "...",
                           paste0(tmp_seq[(round(len/2)+1):length(tmp_seq)], collapse = ""), collapse = "")
                  }) %>% unlist
                  print(data.frame(read_count = object@reads[, 1],
                                   barcode = seqs), row.names = FALSE)
                } else {
                  print(object@reads[, ], row.names = FALSE)
                }
              }
            }
            cat("\n")
            cat(" results dir: \n")
            cat("      ", object@results_dir, "\n")
            cat(" barcode backbone: \n")
            cat("      ", object@BC_backbone, "\n")
            cat(" label: \n")
            cat("      ", object@label, "\n")
          }
)


# methods::setGeneric(
#   name = "reads<-",
#   def = function(object, value) {
#     standardGeneric("reads<-")
#   }
# )
#
#   methods::setMethod(
#     f = "reads<-",
#     signature = "BCdat",
#     definition = function(object, value) {
#       errors <- checkReads(value, NULL)
#       if (length(errors)) {
#           object@reads <- value
#       } else {
#           stop(errors)
#       }
#     }
#   )
#
# methods::setGeneric(
#   name = "resDir",
#   def = function(object) {
#     standardGeneric("resDir")
#   }
# )
#
#   methods::setMethod(
#     f = "resDir",
#     signature = character(),
#     definition = function(object) {
#       return(object@results_dir)
#     }
#   )
#
#
# methods::setGeneric(
#   name = "resDir<-",
#   def = function(object, value) {
#     standardGeneric("resDir<-")
#   }
# )
#
#   methods::setMethod(
#     f = "resDir<-",
#     signature = "BCdat",
#     definition = function(object, value) {
#       errors <- checkReads(value, NULL)
#       if (length(errors)) {
#         object@reads <- value
#       } else {
#         stop(errors)
#       }
#       return(object)
#     }
#   )
#
#
# methods::setGeneric(
#   name = "fileLabel",
#   def = function(object) {
#     standardGeneric("fileLabel")
#   }
# )
#
# methods::setMethod(
#   f = "fileLabel",
#   signature = "BCdat",
#   definition = function(object) {
#     return(object@label)
#   }
# )
#
# methods::setGeneric(
#   name = "fileLabel<-",
#   def = function(object, value) {
#     standardGeneric("fileLabel<-")
#   }
# )
#
#
# methods::setGeneric(
#   name = "backbone",
#   def = function(object) {
#     standardGeneric("backbone")
#   }
# )
#
#
# methods::setMethod(
#   f = "backbone",
#   signature = "BCdat",
#   definition = function(object) {
#     return(object@BC_backbone)
#   }
# )
#
# methods::setGeneric(
#   name = "backbone<-",
#   def = function(object, value) {
#     standardGeneric("backbone<-")
#   }
# )
#
# methods::setReplaceMethod(
#   f = "backbone",
#   signature = "BCdat",
#   definition = function(object, value) {
#     if(methods::is(value, "character")) {
#         object@BC_backbone <- value
#     } else {
#         stop("# the BC_backbone slot demands a character")
#     }
#     return(object)
#   }
# )

# #' @include BCdata-class.R
# NULL

# #' Title READS
# #'
# #' @param object asdf
# #'
# #' @return stuff
# #' @export
# #' @docType BCdat
# #' @rdname reads-BCdat
# methods::setGeneric("reads",
#                     function(object) {
#                       standardGeneric("reads")
#                     }
# )

# #' @rdname reads-BCdat
# #' @aliases reads,BCdat
# #' @exportMethod
# methods::setMethod("reads", signature(object = "BCdat"),
#                    function(object) {
#                      return(object@reads)
#                    }
# )

###################
