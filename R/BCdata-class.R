

#' @title checkBarcodeData
#'
#' @description Checks data slots of BCdat object for correctness.
#'
#' @param object a BCdat object.

checkBarcodeData <- function(object) {

  errors <- character()

    if (!is.data.frame(object@reads)) {
      msg <- paste("Data is of type ", methods::is(object@reads)[1], ".  Should be data.frame.", sep = "")
      errors <- c(errors, msg)
    }

    if (dim(object@reads)[2] > 3) {
      msg <- paste("Data consists of ", dim(object@reads)[2], " column(s).  Should be 2 or 3.", sep = "")
      errors <- c(errors, msg)
    } else {
      if(dim(object@reads)[2] == 3) {
          if(sum(names(object@reads) == c("pos", "read_count", "barcode")) != 3) {
            msg <- paste("Data needs to consist of the columns 'pos', 'read_count' and 'barcode'", sep = "")
            errors <- c(errors, msg)
          }
      }
      if(dim(object@reads)[2] == 2) {
        if(sum(names(object@reads) == c("read_count", "barcode")) != 2) {
          msg <- paste("Data needs to consist of the columns 'read_count' and 'barcode'", sep = "")
          errors <- c(errors, msg)
        }
      }
  }


  if (length(object@results_dir) != 1) {
    msg <- paste("Only one results-path needed.")
    errors <- c(errors, msg)
  }

  if (length(object@label) != 1) {
    msg <- paste("Only one label needed.")
    errors <- c(errors, msg)
  }

  if (length(object@mask) != 1) {
    msg <- paste("Only one mask needed.")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) TRUE else errors
}


#' @title BCdat class
#'
#' @description S4 data class containing every relevant information.
#'
#' @slot reads data.frame containing barcode sequences and their corresponding read counts.
#' @slot results_dir character string of the working directory path.
#' @slot label character string identifying the particular experiment (will be part of the names of any file created).
#' @slot mask character string of the used barcode design.
#'
#' @return a BCdat object.
#' @export

BCdat <- methods::setClass("BCdat", slots = list(
                                       reads = "data.frame",
                                       results_dir = "character",
                                       label = "character",
                                       mask = "character"
                                   ), validity = checkBarcodeData
                          )


methods::setGeneric (
  name = "reads",
  def = function(object) {
    standardGeneric("reads")
  }
)

methods::setMethod (
  f = "reads",
  signature = "BCdat",
  definition = function(object) {
    return(object@reads)
  }
)


methods::setGeneric (
  name = "reads<-",
  def = function(object, value) {
    standardGeneric("reads<-")
  }
)

methods::setReplaceMethod (
  f = "reads",
  signature = "BCdat",
  definition = function(object, value) {
    if(methods::is(value, "data.frame")) {
      if(sum(names(value) %in% c("read_count", "barcode")) != 2) {
        stop("# BCdat demands at least 2 columns named 'read_count' and 'barcode'")
      } else {
        object@reads <- data.frame(pos = 1:dim(value)[1], value)
      }
    } else { stop("# class BCdat demands a data.frame") }
    return(object)
  }
)

methods::setGeneric (
  name = "resDir",
  def = function(object) {
    standardGeneric("resDir")
  }
)

methods::setMethod (
  f = "resDir",
  signature = "BCdat",
  definition = function(object) {
    return(object@results_dir)
  }
)


methods::setGeneric (
  name = "resDir<-",
  def = function(object, value) {
    standardGeneric("resDir<-")
  }
)

methods::setReplaceMethod (
  f = "resDir",
  signature = "BCdat",
  definition = function(object, value) {
    if(methods::is(value, "character")) {
        object@results_dir <- value
    } else {
        stop("# class BCdat@results_dir demands a character")
    }
    return(object)
  }
)

methods::setGeneric (
  name = "fileLabel",
  def = function(object) {
    standardGeneric("fileLabel")
  }
)

methods::setMethod (
  f = "fileLabel",
  signature = "BCdat",
  definition = function(object) {
    return(object@label)
  }
)

methods::setGeneric (
  name = "fileLabel<-",
  def = function(object, value) {
    standardGeneric("fileLabel<-")
  }
)


methods::setReplaceMethod (
  f = "fileLabel",
  signature = "BCdat",
  definition = function(object, value) {
    if(methods::is(value, "character")) {
        object@label <- value
    } else {
        stop("# the label slot demands a character")
    }
    return(object)
  }
)

methods::setGeneric (
  name = "maskBC",
  def = function(object) {
    standardGeneric("maskBC")
  }
)


methods::setMethod (
  f = "maskBC",
  signature = "BCdat",
  definition = function(object) {
    return(object@mask)
  }
)

methods::setGeneric (
  name = "maskBC<-",
  def = function(object, value) {
    standardGeneric("maskBC<-")
  }
)

methods::setReplaceMethod (
  f = "maskBC",
  signature = "BCdat",
  definition = function(object, value) {
    if(methods::is(value, "character")) {
        object@mask <- value
    } else {
        stop("# the mask slot demands a character")
    }
    return(object)
  }
)

setMethod("show", signature = "BCdat",
          definition = function(object){
            cat(" class: BCdat\n\n")
              cat("  number of barcode sequences:", dim(object@reads)[1], "\n")
              if(sum(dim(object@reads)) > 1) {
                  cat("  read count distribution: min", min(object@reads$read_count),
                            " mean", round(mean(object@reads$read_count), digits = 2),
                            " median", stats::median(object@reads$read_count),
                            " max", max(object@reads$read_count), "\n")
                  l <- unique(nchar(as.character(object@reads$barcode)))
                  if(length(l) == 1) {
                    cat("  barcode sequence length:", l[1], "\n")
                  }
                  if(length(l) == 2) {
                    cat("  barcode sequence lengths: ", l[1], ", ", l[2], "\n")
                  }
                  if(length(l) == 3) {
                    cat("  barcode sequence lengths: ", l[1], ", ", l[2], ", ", l[3], "\n")
                  }
                  if(length(l) > 3) {
                    cat("  barcode sequence lengths: ", l[1], ", ", l[2], ", ", l[3], ", ...", "\n")
                  }
                  cat("\n")
                  cat(" barcode read counts:\n")
                  if(dim(object@reads)[1] > 10) {
                      print(object@reads[1:10, 2:3], row.names = FALSE)
                      cat("                   ...")
                  } else {
                      print(object@reads[, 2:3], row.names = FALSE)
                  }
              }
              cat("\n")
              cat(" results dir: \n")
              cat("      ", object@results_dir, "\n")
              cat(" barcode layout: \n")
              cat("      ", object@mask, "\n")
              cat(" label: \n")
              cat("      ", object@label, "\n")
          }
)
