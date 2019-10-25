
#' Accessing the Read-Count slot of a BCdat objects.
#'
#' @param object a BCdat object.
#'
#' @return A data.frame containing the read count table of the object paramter.
#' @examples
#' data(BC_dat)
#' getReads(BC_dat)
#'
#' @export
getReads <- function(object) {
  if (methods::is(object) == "BCdat") {
    return(methods::slot(object, "reads"))
  } else {
    stop("No BCdat object.")
  }
}

#' Replacing the Read-Count slot of a BCdat objects.
#'
#' @param object a BCdat object.
#' @param value a data.frame caontaining two columns called "read_count" and "barcode".
#'
#' @return a BCdat object.
#' @examples
#' data(BC_dat)
#' require("dplyr")
#'
#' bcs <- unlist(lapply(1:20, function(x) {
#'            c("A", "C", "T", "G") %>% sample(replace = TRUE, size = 32) %>% paste0(collapse = "")
#'        }))
#' new_read_count_table <- data.frame(read_count = sample(1:1000, size = 20), barcode = bcs)
#' BC_dat_alt <- setReads(BC_dat, new_read_count_table)
#'
#' @export
setReads <- function(object, value) {

  errors <- checkReads(value, NULL)
  if (length(errors) == 0 & methods::is(object) == "BCdat") {
    methods::slot(object, "reads") <- value
  } else {
    if (methods::is(object) != "BCdat") {
      errors <- c(errors, "Object is not of BCdat class.")
    }
    stop(errors)
  }
  return(object)
}

#' Accessing the Results Directory slot of a BCdat objects.
#'
#' @param object a BCdat object.
#'
#' @return A character string.
#' @examples
#' data(BC_dat)
#' getResultsDir(BC_dat)
#'
#' @export
getResultsDir <- function(object) {

  if (methods::is(object) == "BCdat") {
    return(methods::slot(object, "results_dir"))
  } else {
    stop("This function does only work with BCdat-class objects.")
  }
}

#' Replacing the Results Directory slot of a BCdat objects.
#'
#' @param object a BCdat object.
#' @param value a character string of an existing path.
#'
#' @return a BCdat object.
#' @examples
#' data(BC_dat)
#' new_path <- getwd()
#' BC_dat_alt <- setResultsDir(BC_dat, new_path)
#'
#' @export
setResultsDir <- function(object, value) {

  res <- checkDir(value, NULL)

  if (is.null(res)) {
    if (methods::is(object) == "BCdat") {
      methods::slot(object, "results_dir") <- value
    } else {
      stop("This function does only work with BCdat-class objects.")
    }
  } else {
    stop(res)
  }

  return(object)
}

#' Accessing the Barcode Backbone slot of a BCdat objects.
#'
#' @param object a BCdat object.
#'
#' @return A character string.
#' @examples
#' data(BC_dat)
#' getBackbone(BC_dat)
#'
#' @export
getBackbone <- function(object) {

  if (methods::is(object) == "BCdat") {
    return(methods::slot(object, "BC_backbone"))
  } else {
    stop("This function does only work with BCdat-class objects.")
  }
}

#' Replacing the Barcode Backbone slot of a BCdat objects.
#'
#' @param object a BCdat object.
#' @param value a character string consisting of exclusively IUPAC-nucleotide-code conform letters.
#'
#' @return a BCdat object.
#' @examples
#' data(BC_dat)
#' new_backbone <- getBackboneSelection("BC32-T-Sapphire")
#' BC_dat_alt <- setBackbone(BC_dat, new_backbone)
#'
#' @export
setBackbone <- function(object, value) {

  res <- checkBackbone(value)

  if (is.null(res)) {
    if (methods::is(object) == "BCdat") {
      methods::slot(object, "BC_backbone") <- value
    } else {
      stop("This function does only work with BCdat-class objects.")
    }
  } else {
    stop(res)
  }

  return(object)
}

#' Accessing the Label slot of a BCdat objects.
#'
#' @param object a BCdat object.
#'
#' @return A character string.
#' @examples
#' data(BC_dat)
#' getLabel(BC_dat)
#'
#' @export
getLabel <- function(object) {

  if (methods::is(object) == "BCdat") {
    return(methods::slot(object, "label"))
  } else {
    stop("This function does only work with BCdat-class objects.")
  }
}

#' Replacing the Label slot of a BCdat objects.
#'
#' @param object a BCdat object.
#' @param value a character string.
#'
#' @return a BCdat object.
#' @examples
#' data(BC_dat)
#' new_label <- "foo-bar"
#' BC_dat_alt <- setLabel(BC_dat, new_label)
#'
#' @export
setLabel <- function(object, value) {

  res <- checkLabel(value)

  if (is.null(res)) {
      if (methods::is(object) == "BCdat") {
            methods::slot(object, "label") <- value
      } else {
            stop("This function does only work with BCdat-class objects.")
      }
  } else {
    stop(res)
  }

  return(object)
}

