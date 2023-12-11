
#' Clustered HD Plot
#'
#' This function will create a jitter plot displaying the maximal distances within each of the barcode sequence clusters.
#'
#' @param datEC a BC_dat object, returned by the error_correction function with the EC_analysis parameter set to TRUE.
#' @param size a numeric value, specifying the dot size.
#'
#' @return a ggplot2 object.
#' @export
#'
# @examples
# data(BC_dat)
# BC_dat_EC <- errorCorrection(BC_dat, maxDist = 8, EC_analysis = TRUE)
# error_correction_clustered_HDs(BC_dat_EC, size = 0.75)

error_correction_clustered_HDs <- function(datEC, size = 0.75) {

  tmp <- reshape2::melt(
    lapply(datEC$EC_seqs, function(x) {
      if(length(x) > 1) {
          max(as.matrix(stringdist::stringdistmatrix(x, method = "hamming")))
      } else {
        NA
      }
    }))

  tmp <- tmp[!is.na(tmp$value), ]

  no <- sum(unlist(lapply(datEC$EC_seqs, function(x) {
    if(length(x) > 1) {
      1
    }
  })))

  if(no == 0) {
    warning("# There were no sequences clustered during error correction.")
    p <- NULL
  } else {

    p <- ggplot2::ggplot(tmp, ggplot2::aes_string(x = factor("EC-BCs"), y = "value")) +
          ggplot2::geom_jitter(width = 0.2, height = 0.2, size = size) +
          ggplot2::ylab("max. HDs") +
          ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
          ggplot2::xlab(paste0("clustered BC distances (", no," clusters)")) +
          ggplot2::scale_y_continuous(limits = c(0, max(tmp$value, na.rm = TRUE) + 1))
  }

  return(p)
}

#' Circle Plot
#'
#' creates a circle plot based on the additional data gathered by the error_correction function (EC_analysis needs to be set to TRUE). This function is intended to visualize the error correction procedure.
#'
#' @param edges a data frame containing edge definitions by two columns calles "from" and "to". Such data frame will be returned by the error_correction function with the EC_analysis parameter set to TRUE.
#' @param vertices a data frame with at least one column containing a list of nodes (also returned by the error_correction function with the EC_analysis parameter set to TRUE)
#'
#' @return a ggplot2 object.
#' @export
#'
# @examples
# data(BC_dat)
# BC_dat_EC <- errorCorrection(BC_dat, maxDist = 8, EC_analysis = TRUE)
# error_correction_circlePlot(edges = BC_dat_EC$edges, vertices = BC_dat_EC$vertices)

error_correction_circlePlot <- function(edges, vertices) {

  fromTo <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

  p <- ggraph::ggraph(fromTo, 'circlepack') +
    ggraph::geom_node_circle(ggplot2::aes_string(fill = "read_counts"), n = 50) +
    ggplot2::coord_fixed() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   line = ggplot2::element_blank(), rect = ggplot2::element_blank()) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Read Counts"))

  return(p)
}

#' Tree Plot
#'
#' creates a Tree Plot visualising of the barcode clustering as part of the error correction process.
#'
#' @param edges a data frame containing edge definitions by two columns calles "from" and "to". Such data frame will be returned by the error_correction function with the EC_analysis parameter set to TRUE.
#' @param vertices a data frame with at least one column containing a list of nodes (also returned by the error_correction function with the EC_analysis parameter set to TRUE)
#'
#' @return a ggplot2 object.
#' @export
#'
# @examples
# data(BC_dat)
# BC_dat_EC <- errorCorrection(BC_dat, maxDist = 8, EC_analysis = TRUE)
# error_correction_treePlot(edges = BC_dat_EC$edges, vertices = BC_dat_EC$vertices)

error_correction_treePlot <- function(edges, vertices) {

  fromTo <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

  p <- ggraph::ggraph(fromTo, layout = 'partition', circular = TRUE) +
    ggraph::geom_edge_diagonal(lineend = 'round') +
    ggraph::scale_edge_width(range = c(0.2, 1.5)) +
    ggraph::geom_node_point(ggplot2::aes_string(colour = "depth", size = "read_counts")) +
    ggplot2::coord_fixed() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   line = ggplot2::element_blank(),
                   rect = ggplot2::element_blank()) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::guides(colour = FALSE,
                    size = FALSE)

  return(p)
}


#' Plotting a VennDiagram
#'
#' @description plotVennDiagramm will create a Venn Diagram ans is based on the VennDiagram package.
#' It accepts a list of BCdat objects and will return a ggplot2 output object.
#'
#' @param BC_dat a list of BCdat objects.
#' @param alpha_value color transparency value [0-1].
#' @param colrs a character vector containing the desired colors, if NA the colors will be chosen automatically.
#' @param border_color a character value specifying the desired border color, if NA no border will be drawn.
#' @param plot_title a character value.
#' @param legend_sort a character or factor vector in case the order of legend items needs to be changed.
#' @param annotationSize an integer value specifying the venn diagramm internal text size.
#'
#' @return ggplot2 object.
#' @export

plotVennDiagram <- function(BC_dat, alpha_value = 0.4, colrs = NA, border_color = NA, plot_title = "", legend_sort = NULL, annotationSize = 5) {

  futile.logger::flog.threshold('FATAL', name = "VennDiagramLogger")

  if(length(BC_dat) > 1) {
    v <- VennDiagram::venn.diagram(lapply(BC_dat, function(x) methods::slot(x, "reads")$"barcode"),
                    filename = NULL,
                    category.names = unlist(lapply(BC_dat, function(x) methods::slot(x, "label"))))
  } else {
    if(methods::is(BC_dat) == "BCdat") {
      v <- VennDiagram::venn.diagram(list(methods::slot(BC_dat, "reads")$"barcode"),
                      filename = NULL,
                      category.names = methods::slot(BC_dat, "label"))
    } else {
      warning("plotVennDiagramm needs a list of or a single BCdat object!")
      return(1)
    }
  }

  l <- sum(unlist(lapply(unclass(v), methods::is)) == "polygon")

  coord_x <- coord_y <- cate <- NULL
  for(i in 1:ifelse(l == 2, l, (l/2))) {
    coord_x <- c(coord_x, unclass(unclass(v)[[i]])$x)
    coord_y <- c(coord_y, unclass(unclass(v)[[i]])$y)
    if(length(BC_dat) > 1) {
      cate <- c(cate, rep(methods::slot(BC_dat[[i]], "label"), length(unclass(unclass(v)[[i]])$x)))
    } else {
      cate <- c(cate, rep(methods::slot(BC_dat, "label"), length(unclass(unclass(v)[[i]])$x)))
    }
  }

  dat <- data.frame(x = coord_x, y = coord_y, category = cate, stringsAsFactors = TRUE)

  if(length(legend_sort) != 0) {
    if(length(unique(dat$category)) != length(legend_sort)) {
      warning("Number of legend_sort entries does not match the number of BC_dat objects!")
    } else {
      dat$category <- factor(dat$category, levels = legend_sort)
    }
  }

  vennPlot <-
    ggplot2::ggplot(dat, ggplot2::aes_string(x = "x", y = "y", fill = "category")) +
    ggplot2::geom_polygon(alpha = alpha_value, linetype = "solid", colour = border_color) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.001),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::ggtitle(label = plot_title)

  if(sum(is.na(colrs)) == 0) {
    if(length(colrs) == length(BC_dat)) {
      vennPlot <- vennPlot + ggplot2::scale_fill_manual(values = colrs)
    } else {
      if(length(colrs) == 1) {
        vennPlot <- vennPlot + ggplot2::scale_fill_manual(values = rep(colrs, length(BC_dat)))
      } else {
        warning("Number of colors has to match the length of the data list!")
      }
    }
  }

  l <- which(unlist(lapply(unclass(v), methods::is)) == "text")
  l <- l[1:(length(l) - length(BC_dat))]

  for(i in l) {
    vennPlot <- vennPlot + ggplot2::annotate(geom = "text",
                                    x = as.numeric(unclass(unclass(v)[[i]])$x),
                                    y = as.numeric(unclass(unclass(v)[[i]])$y),
                                    label = unclass(unclass(v)[[i]])$label,
                                    color = "black", size = annotationSize)
  }

  return(vennPlot)
}


#' @title Plotting a Distance Network
#' @description plotDistanceVisNetwork will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the ggplot2 and the ggnetwork packages. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a distance between two barcodes
#' of maximal \code{minDist}.
#' If \code{ori_BCs} is provided the node color also refelects the distance of a particular barcode to one of the given
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param minDist an integer value representing the maximal distance value for which the graph will
#' contain edges.
#' @param loga a logical value indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#'
#' @return a visNetwork object.
#' @export
#'
# @examples
# data(BC_dat)
# plotDistanceVisNetwork(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL,
# complete = FALSE, col_type = "rainbow")

plotDistanceVisNetwork <- function(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL, complete = FALSE, col_type = "rainbow", m = "hamming") {

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(visNetwork::visNetwork(data.frame(), data.frame()))
  }

  net <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m)
  if(complete) {
    net <- t(apply(net, 1, function(x) {
      x[x > min(x[x > 0])] <- 0
      return(x)
    }))
    weight <- net
    net[net > 0] <- 1
  } else {
    net[net > minDist] <- 0
  }

  net[.getDiagonalIndex(dim(net)[1])] <- 0
  net <- igraph::graph.adjacency(net)

  if (loga) {
    v <- log(methods::slot(BC_dat, "reads")$"read_count")
  } else {
    v <- methods::slot(BC_dat, "reads")$"read_count"
  }

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(visNetwork::visNetwork(data.frame(), data.frame()))
  }

  if (length(ori_BCs) > 0) {

    minDists <- .getMinDist(BC_dat, ori_BCs, m)
    colrs <- .generateColors(minDists, type = col_type)

    nodes <- data.frame(id = 1:length(igraph::V(net)),
                        label = paste0(methods::slot(BC_dat, "reads")$"barcode", " - ",
                                       methods::slot(BC_dat, "reads")$"read_count"),
                        title = paste0(methods::slot(BC_dat, "reads")$"barcode", " - ",
                                       methods::slot(BC_dat, "reads")$"read_count"),
                        value = v,
                        shadow = TRUE,
                        color.background = as.character(colrs[[1]]),
                        color.border = as.character(colrs[[1]]),
                        color.highlight.border = "gray", stringsAsFactors = TRUE)
  } else {
    nodes <- data.frame(id = 1:length(igraph::V(net)),
                        label = paste0(methods::slot(BC_dat, "reads")$"barcode", " - ",
                                       methods::slot(BC_dat, "reads")$"read_count"),
                        title = paste0(methods::slot(BC_dat, "reads")$"barcode", " - ",
                                       methods::slot(BC_dat, "reads")$"read_count"),
                        value = v,
                        shadow = TRUE,
                        color.background = "black",
                        color.border = "darkgray",
                        color.highlight.border = "gray", stringsAsFactors = TRUE)
  }

  edges <- data.frame(from = igraph::get.edgelist(net)[, 1], to = igraph::get.edgelist(net)[, 2], stringsAsFactors = TRUE)

  graph_visNetwork <- visNetwork::visNetwork(nodes, edges, width = "100%")
  graph_visNetwork <- visNetwork::visOptions(graph_visNetwork,
                                             highlightNearest = list(enabled = TRUE),
                                             nodesIdSelection = TRUE) #list(enabled = TRUE, useLabels = TRUE)
  graph_visNetwork <- visNetwork::visIgraphLayout(graph_visNetwork, layout = "layout_nicely", smooth = TRUE)
  graph_visNetwork <- visNetwork::visHierarchicalLayout(graph_visNetwork, graph_visNetwork)

  return(graph_visNetwork)
}

#' @title Plotting a Distance Network (error correction)
#' @description plotDistanceVisNetwork will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the ggplot2 and the ggnetwork packages. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a distance between two barcodes
#' of maximal \code{minDist}.
#' If \code{ori_BCs} is provided the effects of the error correction function will be color-coded only for those sequences.
#'
#' @param BC_dat a BCdat object.
#' @param BC_dat_EC the corresponding error corrected BCdat object (EC_analysis has to be TRUE)
#' @param minDist an integer value representing the maximal distance value for which the graph will
#' contain edges.
#' @param loga a logical value indicating the use or non-use of logarithmic read count values.
#' @param equal_node_sizes a logical value. If TRUE, every node will have the sames size.
#' @param BC_threshold an integer value representing the number of barcodes for which the color-coding should be applied (starting with the barcodes with the most read counts).
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#'
#' @return a visNetwork object.
#' @export

plotDistanceVisNetwork_EC <- function(BC_dat, BC_dat_EC, minDist = 1, loga = TRUE, equal_node_sizes = TRUE, BC_threshold = NULL, ori_BCs = NULL, complete = FALSE, col_type = "rainbow", m = "hamming") {

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(visNetwork::visNetwork(data.frame(), data.frame()))
  }

  colrsEC <- rep("#000000",  length(methods::slot(BC_dat, "reads")$barcode))
  if(is.null(BC_threshold)) {
    th <- length(BC_dat_EC$EC_seqs):1
  } else {
    th <- BC_threshold:1
  }
  if(!is.null(ori_BCs)) {
    th <- which(ori_BCs %in% methods::slot(BC_dat_EC$BC_dat, "reads")$barcode)
  }

  index <- NULL
  for(i in th) {
    index <- which(methods::slot(BC_dat, "reads")$barcode %in% rev(BC_dat_EC$EC_seqs)[[i]])
    colrsEC[index] <- RColorBrewer::brewer.pal(12, "Set3")[round(stats::runif(1, min = 1, max = 12))]
  }

  net <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m)
  if(complete) {
    net <- t(apply(net, 1, function(x) {
      x[x > min(x[x > 0])] <- 0
      return(x)
    }))
    weight <- net
    net[net > 0] <- 1
  } else {
    net[net > minDist] <- 0
  }

  net[.getDiagonalIndex(dim(net)[1])] <- 0

  net <- igraph::graph.adjacency(net)

  if(equal_node_sizes) {
    v <- 3
  } else {
    if(loga) {
      v <- log(methods::slot(BC_dat, "reads")$"read_count")
    } else {
      v <- methods::slot(BC_dat, "reads")$"read_count"
    }
  }

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(visNetwork::visNetwork(data.frame(), data.frame()))
  }

  nodes <- data.frame(id = 1:length(igraph::V(net)),
                        label = paste0(methods::slot(BC_dat, "reads")$"barcode", " - ",
                                       methods::slot(BC_dat, "reads")$"read_count"),
                        title = paste0(methods::slot(BC_dat, "reads")$"barcode", " - ",
                                       methods::slot(BC_dat, "reads")$"read_count"),
                        value = v,
                        shadow = TRUE,
                        color.background = colrsEC,
                        color.border = colrsEC,
                        color.highlight.border = "gray", stringsAsFactors = TRUE)

  edges <- data.frame(from = igraph::get.edgelist(net)[, 1], to = igraph::get.edgelist(net)[, 2], stringsAsFactors = TRUE)

  graph_visNetwork <- visNetwork::visNetwork(nodes, edges, width = "100%")
  graph_visNetwork <- visNetwork::visOptions(graph_visNetwork, highlightNearest = TRUE)
  graph_visNetwork <- visNetwork::visIgraphLayout(graph_visNetwork, "layout_nicely")
  graph_visNetwork <- visNetwork::visHierarchicalLayout(graph_visNetwork, graph_visNetwork)

  return(graph_visNetwork)
}


#' @title Plotting a Distance Network (error correction)
#'
#' @description ggplotDistanceGraph will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the ggplot2 and the ggnetwork packages. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a distance between two barcodes
#' of maximal \code{minDist}.
#' If \code{ori_BCs} is provided the node color also refelects the distance of a particular barcode to one of the initial
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param BC_dat_EC the error corrected BCdat object (the EC_analysis parameter needs to be set to TRUE).
#' @param minDist an integer value representing the maximal distance for which the graph will
#' contain edges.
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param equal_node_sizes a logical value. If TRUE, every node will have the same size.
#' @param BC_threshold a nnumeric value, limiting the number of barcodes for which their error correction "history" will be colored (if BC_threshold = 5 then the five biggest barcodes will be evaluated)
#' @param ori_BCs a vector of character strings containing barcode sequences (without the fixed positions of the barcode construct). Similar to BC_threshold but allowing for barcode identification via sequence.
#' @param lay a character string, identifying the prefered layout algorithm (see ggnetwork layout option).
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#' @param scale_nodes a numeric value, scaling the node size.
#' @param scale_edges a numeric value, scaling the edge size.
#'
#' @return a ggplot2 object
#' @export

ggplotDistanceGraph_EC <- function(BC_dat, BC_dat_EC, minDist = 1, loga = TRUE, equal_node_sizes = TRUE, BC_threshold = NULL, ori_BCs = NULL, lay = "fruchtermanreingold", complete = FALSE, col_type = "rainbow", m = "hamming", scale_nodes = 1, scale_edges = 1) {

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(ggplot2::ggplot())
  }

  colrsEC <- rep("#000000",  length(methods::slot(BC_dat, "reads")$barcode))
  if (is.null(BC_threshold)) {
    th <- length(BC_dat_EC$EC_seqs):1
  } else {
    th <- BC_threshold:1
  }
  if (!is.null(ori_BCs)) {
    th <- which(ori_BCs %in% methods::slot(BC_dat_EC$BC_dat, "reads")$barcode)
  }

  index <- NULL
  for (i in th) {
    index <- which(methods::slot(BC_dat, "reads")$barcode %in% rev(BC_dat_EC$EC_seqs)[[i]])
    colrsEC[index] <- c(RColorBrewer::brewer.pal(12, "Set3"),
                        RColorBrewer::brewer.pal(8, "Dark2"),
                        RColorBrewer::brewer.pal(9, "Pastel1"))[round(stats::runif(1, min = 1, max = 29))]
  }

  net <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m)
  if (complete) {
    net <- t(apply(net, 1, function(x) {
      x[x > min(x[x > 0])] <- 0
      return(x)
    }))
    net[net > 0] <- 1
  } else {
    net[net > minDist] <- 0
  }

  if (sum(rowSums(net)) <= 2) {
    warning("# Only 1 or 0 edges left, due to a ggnetwork problem (version 0.5.1), also graphs with only one edge can't be displayed")
    return(ggplot2::ggplot() +
             ggplot2::theme_minimal() +
             ggplot2::ggtitle("no edges left") +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }

  net <- network::network(net, directed = FALSE)

  if (equal_node_sizes) {
    v_size <- 3
  } else {
      if(loga) {
        v_size <- log(methods::slot(BC_dat, "reads")$"read_count")
      } else {
        v_size <- methods::slot(BC_dat, "reads")$"read_count"
      }
  }
  network::set.vertex.attribute(net, "size", v_size * scale_nodes)
  network::set.vertex.attribute(net, "names", as.character(methods::slot(BC_dat, "reads")$"barcode"))

  if (!is.null(ori_BCs)) {
    minDists <- .getMinDist(BC_dat, ori_BCs, m)
    network::set.vertex.attribute(net, "minDists", minDists)
  } else {
    network::set.vertex.attribute(net, "minDists", NA)
  }

  network::set.vertex.attribute(net, "colors", colrsEC)

  net <- ggnetwork::ggnetwork(net, layout = lay, cell.jitter = 0.75)
  net$minDists <- as.factor(net$minDists)

  p <- ggplot2::ggplot(net, ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
    ggnetwork::geom_edges(color = "grey50", size = 0.3 * scale_edges) +
    ggnetwork::theme_blank()

  p <- p + ggnetwork::geom_nodes(size = v_size * scale_nodes, color = colrsEC)

  return(p)
}


#' @title Plotting a Distance Network
#'
#' @description ggplotDistanceGraph will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the ggplot2 and the ggnetwork packages. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a distance between two barcodes
#' of maximal \code{minDist}.
#' If \code{ori_BCs} is provided the node color also refelects the distance of a particular barcode to one of the initial
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param minDist an integer value representing the maximal distance for which the graph will
#' contain edges.
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param lay a character string, identifying the prefered layout algorithm (see ggnetwork layout option, "?gplot.layout"). Default value is "fruchtermanreingold", but possible are also "circle", "eigen", "kamadakawai", "spring" and many more. Or the user provides a two-column matrix with as many rows as there are nodes in the network, in which case the matrix is used as nodes coordinates.
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes ("rainbow", "heat.colors", "topo.colors", "greens", "wild" - see package "grDevices").
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#' @param scale_nodes a numeric value, scaling the node size.
#' @param scale_edges a numeric value, scaling the edge size.
#' @param legend_size a numeric value, scaling the legend symbol size, if legend_size equals 0, the legend will be dismissed.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' data(BC_dat)
#' ggplotDistanceGraph(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL, lay = "fruchtermanreingold",
#' complete = FALSE, col_type = "rainbow")
#'
#' }
#'

ggplotDistanceGraph <- function(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL, lay = "fruchtermanreingold", complete = FALSE, col_type = "rainbow", m = "hamming", scale_nodes = 1, scale_edges = 1, legend_size = 4) {

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(ggplot2::ggplot())
  }

  net <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m)
  if (complete) {
    net <- t(apply(net, 1, function(x) {
                x[x > min(x[x > 0])] <- 0
                return(x)
              }))
    net[net > 0] <- 1
  } else {
    net[net > minDist] <- 0
  }

  if (sum(rowSums(net)) <= 2) {
    warning("# Only 1 or 0 edges left, due to a ggnetwork problem (version 0.5.1), also graphs with only one edge can't be displayed")
    return(ggplot2::ggplot() +
             ggplot2::theme_minimal() +
             ggplot2::ggtitle("no edges left") +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }

  net <- network::network(net, directed = FALSE)

  if (loga) {
    v_size <- log(methods::slot(BC_dat, "reads")$"read_count")
  } else {
    v_size <- methods::slot(BC_dat, "reads")$"read_count"
  }

  network::set.vertex.attribute(net, "size", v_size * scale_nodes)
  network::set.vertex.attribute(net, "names", as.character(methods::slot(BC_dat, "reads")$"barcode"))

  if (!is.null(ori_BCs)) {
    minDists <- .getMinDist(BC_dat, ori_BCs, m)
    colrs <-  .generateColors(minDists, type = col_type)

    network::set.vertex.attribute(net, "minDist", minDists)
    network::set.vertex.attribute(net, "colors", colrs[[1]])
    tmp <- TRUE
  } else {
    network::set.vertex.attribute(net, "minDist", NA)
    tmp <- FALSE
  }

  net <- ggnetwork::ggnetwork(net, layout = lay, cell.jitter = 0.75)
  net$minDist <- as.factor(net$minDist)

  p <- ggplot2::ggplot(net, ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend", text = "names")) +
        ggnetwork::geom_edges(color = "grey50", size = 0.3 * scale_edges) +
        ggnetwork::theme_blank()

  if (tmp) {
    p <- p + ggnetwork::geom_nodes(size = v_size * scale_nodes, ggplot2::aes_string(color = "minDist")) +
              ggplot2::scale_colour_manual(name = ifelse(m == "hamming", "HDs", "min. Distances"), values = colrs[[2]]) +
              if(legend_size == 0) {
                ggplot2::guides(colour = "none")
              } else {
                ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legend_size)))
              }
  } else {
    p <-  p + ggnetwork::geom_nodes(size = v_size * scale_nodes)
  }

  return(p)
}

#' @title Plotting a Distance Network
#'
#' @description plotDistanceIgraph will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the igraph package. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a distance between two barcodes
#' of maximal \code{minD}.
#' If \code{ori_BCs} is provided the node color also refelects the distance of a particular barcode to one of the initial
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param minDist an integer value representing the maximal distance value for which the graph will
#' contain edges.
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param threeD a logical value to chose between 2D and 3D visualisation.
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#' @param leg_pos a character string, containing the position of the legend (e.g. topleft), if NULL no legend will be plotted
#' @param inset a numeric value, specifying the distance from the margins as a fraction of the plot region
#' @param title a character string, containing the legend title
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#'
#' @return an igraph object.
#' @export

plotDistanceIgraph <- function(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL, threeD = FALSE, complete = FALSE, col_type = "rainbow", leg_pos = "left", inset = -0.125, title = "Distance", m = "hamming") {

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(igraph::plot.igraph(igraph::graph.adjacency(0), xlim = c(0,0), ylim = c(0,0)))
  }

  adj_dat <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m)

  if(complete) {
    adj_dat <- t(apply(adj_dat, 1, function(x) {
      x[x > min(x[x > 0])] <- 0
      return(x)
    }))
    adj_dat[adj_dat > 0] <- 1
  } else {
    adj_dat[adj_dat > minDist] <- 0
  }

  colnames(adj_dat) <- methods::slot(BC_dat, "reads")$"barcode"

  g <- igraph::graph.adjacency(adj_dat, mode = c("undirected"), diag = FALSE, add.colnames = "names")

  igraph::V(g)$size <- log(methods::slot(BC_dat, "reads")$"read_count")
  igraph::V(g)$label <- ""

  if(length(ori_BCs) > 0) {

    minDists <- .getMinDist(BC_dat, ori_BCs, m)
    colrs <- .generateColors(minDists, type = col_type)

    igraph::V(g)$color <- colrs[[1]]
  } else {
    igraph::V(g)$color <- "gray"
  }

  if(!threeD) {
    igraph::plot.igraph(g, layout = igraph::layout_nicely) # layout_with_fr
    if(length(ori_BCs) > 0 & !is.null(leg_pos)) {
      graphics::legend(leg_pos, legend = names(colrs[[2]]), fill = colrs[[2]], bty = "n", title = title, inset = inset)
    }
  } else {
    igraph::rglplot(g, layout = igraph::layout_nicely(graph = g, dim = 3))
  }
}

#' @title Plotting Time Series Data
#'
#' @description Uses the result of the generateTimeSeriesData function as inout and generates a visualisation of the clonal
#' contributions over a number of given time points (similar to a stacked barplot).
#'
#' @param ov_dat a numeric matrix consisting of all time points as columns and all barcode sequences as rows and the corresponding read counts as numerical values (see function \code{generateTimeSeriesData()}).
#' @param colr a vector of character strings identifying a certain color palette.
#' @param tp a numeric vector containing the time points of measurement (in case of unequally distributed time points).
# @param bar_width a numeric value, specifiying the bar width at every time pime point.
#' @param x_label a character string providing the x-axis label.
#' @param y_label a character string providing the y-axis label.
#' @return a ggplot2 object.
#'
#' @export
#' @examples
#' ov_dat <- matrix(round(runif(1:100, min = 0, max = 1000)), ncol = 5)
#' rownames(ov_dat) <- paste("barcode", 1:20)
#' plotTimeSeries(ov_dat)

plotTimeSeries <- function(ov_dat, colr = NULL, tp = NULL, x_label = "time", y_label = "contribution") {

  if(!is.matrix(ov_dat)) {
    stop("# plotTimeSeries requires a numeric matrix! Try generateTimeSeriesData() first.")
  }

  if(sum(ov_dat[, 1]) != 1) {
    message("# normalization of read count data")
    tmp <- colSums(ov_dat)
    tmp[tmp == 0] <- 1
    ov_dat <- t(t(ov_dat) / tmp)
  }

  if(is.null(colr)) {
    colr <- c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))

  } else {
    if(length(colr) != dim(ov_dat)[1]) {
      stop("# wrong number of colors defined")
    }
  }

  if(!is.null(tp)) {
    if(length(tp) != dim(ov_dat)[2]) {
      stop("# tp contains less values than the ov_dat data matrix contains columns")
    }
    colnames(ov_dat) <- tp
  } else {
    tp <- 1:dim(ov_dat)[2]
  }

  ov_dat <- reshape2::melt(ov_dat)
  colnames(ov_dat) <- c("barcodes", "time", "frequency")
  ov_dat$barcodes <- as.factor(ov_dat$barcodes)
  ov_dat$time <- as.numeric(ov_dat$time)
  ov_dat$frequency <- as.numeric(ov_dat$frequency)

  if(length(colr) == dim(ov_dat)[1]) {
    coolr <- colr
  } else {
    coolr <- sample(colr, dim(ov_dat)[1], replace = TRUE)
  }

  ggOV <- ggplot2::ggplot(ov_dat) +
    ggplot2::theme_light() +
    ggplot2::geom_area(ggplot2::aes_string(group = "barcodes", fill = "barcodes", x = "time", y = "frequency")) +
    ggplot2::scale_fill_manual(values = coolr) +
    ggplot2::scale_color_manual(values = coolr) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab(y_label) + ggplot2::xlab(x_label) + #ggplot2::xlim(tp[1], tp[length(tp)]) +
    #ggplot2::scale_x_discrete(breaks=tp, labels=tp, limits=tp)
    ggplot2::scale_x_continuous(breaks = 1:length(tp), labels = tp)

  return(ggOV)
}

#' @title Creating a gdf File
#' @description createGDF creates a data file usable with the free graph visualisation tool gephi. The nodes
#' represent barcodes and its respective size reflects the corresponding read counts. Edges between nodes indicate
#' a distance between two barcodes of maximal \code{minD}.
#' If \code{ori_BCs} is provided the node color refelects the distance of a particular barcode to one
#' of the provided barcode sequences.
#'
#' @param BC_dat a BCdat object.
#' @param minDist an integer value representing the maximal distance value for which the graph will
#' contain edges.
#' @param loga a logical value indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param col_type character sting, choosing one of the available color palettes.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#'
#' @return NULL
#'
#' @export
#' @examples
#'
#' \dontrun{
#'
#' data(BC_dat)
#' createGDFFile(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL, col_type = "rainbow")
#'
#' }

createGDF <- function(BC_dat, minDist = 1, loga = TRUE, ori_BCs = NULL, col_type = "rainbow", m = "hamming") {

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(NULL)
  }

  file_name <- paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), ".gdf", sep = "")

  adj_dat <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m)
  adj_dat[adj_dat > minDist] <- 0

  colnames(adj_dat) <- methods::slot(BC_dat, "reads")$"barcode"
  g <- igraph::graph.adjacency(adj_dat, mode=c("undirected"), diag=FALSE, add.colnames = "names")

  if(loga) {
    tmp_reads <- log(as.numeric(methods::slot(BC_dat, "reads")$"read_count"))
  } else {
    tmp_reads <- as.numeric(methods::slot(BC_dat, "reads")$"read_count")
  }

  final <- rbind(c("nodedef>name VARCHAR","label VARCHAR","size DOUBLE"),
                 cbind( paste("node_",1:dim(methods::slot(BC_dat, "reads"))[1],sep = ""),
                        as.character(methods::slot(BC_dat, "reads")$"barcode"),
                        tmp_reads))

  if(length(ori_BCs) > 0) {
    minDists <- .getMinDist(BC_dat, ori_BCs, m)
    colrs <- .hex2rgbColor(.generateColors(minDists, type = col_type)[[1]])
    final <- cbind(final, c("color VARCHAR", colrs))
  }

  utils::write.table(final, file = file_name, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)


  if(length(unclass(g)[[3]]) > 0) {
    graph_data <- cbind(paste("node_",(unclass(g)[[4]]+1), sep = ""), paste("node_",(unclass(g)[[3]]+1), sep = ""), rep("FALSE",length(unclass(g)[[3]])))
  } else {
    graph_data <- NULL
  }

  utils::write.table(rbind(c("edgedef>node1 VARCHAR","node2 VARCHAR", "directed BOOLEAN"), graph_data), file = file_name, append = TRUE, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
}

#' @title Plotting a Cluster Tree
#' @description Generates a tree plot based on a herachical clustering of the complete distance matrix.
#'
#' @param BC_dat a BCdat object.
#' @param tree_est a character string, indicating the particular cluster algorithm, possible algorithms are "Neighbor-Joining" ("NJ") and "Unweighted Pair Group Method" ("UPGMA").
#' @param type a character string, the graph layout style ("unrooted", "phylogram", "cladogram", "fan", "radial").
#' @param tipLabel a logical value, indicating the use of labeled tree leaves.
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#'
#' @export

plotClusterTree <- function(BC_dat, tree_est = "NJ", type = "unrooted", tipLabel = FALSE, m = "hamming") {

  if(dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(NULL)
  }

  if(dim(methods::slot(BC_dat, "reads"))[1] <= 2) {
    return(graphics::plot(0,0, axes = FALSE, xlim = c(1,1), ylim =c(1,1), xlab = "", ylab = "", main = "at least 2 barcode sequences needed"))
  }

  hamming_dist <- stats::as.dist(stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m))
  hamming_phylo <- ape::as.phylo(stats::hclust(hamming_dist))

  if(tree_est == "UPGMA") {
    tree_data  <- phangorn::upgma(hamming_dist)
  }
  if(tree_est == "NJ") {
    tree_data  <- phangorn::NJ(hamming_dist)
  }

  tree_data$tip.label <- as.character(methods::slot(BC_dat, "reads")$"barcode")
  graphics::plot(tree_data, type, show.tip.label = tipLabel)
}

#' @title Plotting a Cluster ggTree
#' @description Generates a tree plot based on a herachical clustering of the complete distance matrix.
#'
#' @param BC_dat a BCdat object.
#' @param tree_est a character string, indicating the particular cluster algorithm, possible algorithms are "Neighbor-Joining" ("NJ") and "Unweighted Pair Group Method" ("UPGMA").
#' @param type a character string, the graph layout style ('rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle' or 'daylight').
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information).
#'
#' @return a ggtree object.
#' @export
#' @examples
#' \dontrun{
#' data(BC_dat)
#' plotClusterGgTree(BC_dat, tree_est = "UPGMA", type = "circular")
#' }

plotClusterGgTree <- function(BC_dat, tree_est = "NJ", type = "rectangular", m = "hamming") {

  if(dim(methods::slot(BC_dat, "reads"))[1] <= 2) {
    return(ggplot2::ggplot() +
             ggplot2::theme_minimal() +
             ggplot2::ggtitle("at least 2 barcode sequences needed") +
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }

  hamming_dist <- stats::as.dist(stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = m))
  hamming_phylo <- ape::as.phylo(stats::hclust(hamming_dist))

  if(tree_est == "UPGMA") {
    tree_data  <- phangorn::upgma(hamming_dist)
  }
  if(tree_est == "NJ") {
    tree_data  <- phangorn::NJ(hamming_dist)
  }

  tree_data$tip.label <- as.character(methods::slot(BC_dat, "reads")$"barcode")
  ggtree::ggtree(tree_data, layout = type)

}

#' @title Plotting a Kirchenplot
#' @description Generates a barplot based on read counts. If \code{ori_BCs} is provided the bar color reflects the
#' distance between a particular barcode to one of the provided barcode sequences.
#'
#' @param BC_dat a BCdat object.
#' @param ori_BCs a vector of character strings containing known barcode sequences (without the fixed positions of the barcode construct).
#' @param ori_BCs2 a vector of character strings containing a 2nd set of known barcode sequences (also without the fixed positions).
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param col_type character string, choosing one of the availabe color palettes ("rainbow", "heat.colors", "topo.colors", "greens", "wild" - see package "grDevices")
#' @param m a character string, Method for distance calculation, default value is Hamming distance. Possible values
#' are "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex" (see stringdist function
#' of the stringdist-package for more information). If neither `ori_BCs` nor `ori_BCs2` are provided with input the choice of `m` does not matter.
#' @param setLabels a character vector, containing three strings serving as plot labels.
#'
#' @return a ggplot2 object
#' @export
generateKirchenplot <- function(BC_dat, ori_BCs = NULL, ori_BCs2 = NULL, loga = TRUE, col_type = NULL, m = "hamming", setLabels = c("BC-Set 1", "Rest", "BC-Set 2")) {

  if (length(m) != 1) {
    warning("Only one method for distance calculation can be utilized. Therefore only the first will be used")
    m <- m[1]
  } else {
    if (sum(m %in% c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")) != 1) {
      stop("Unknown method for distance calculation chosen.")
    }
  }

  if (!is.null(ori_BCs2)) {
    generateKirchenplot_separate(BC_dat, ori_BCs1 = ori_BCs, ori_BCs2 = ori_BCs2, loga, col_type, setLabels, m) + ggplot2::theme_light()
  } else {
    generateKirchenplot_single(BC_dat, ori_BCs, loga, col_type, m) + ggplot2::theme_light()
  }
}

generateKirchenplot_single <- function(BC_dat, ori_BCs = NULL, loga = TRUE, col_type = NULL, method) {

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(ggplot2::ggplot())
  }

  if (is.null(col_type)) {
    col_type <- "darkgray"
  }

  data_dim <- dim(methods::slot(BC_dat, "reads"))[1]

  ggbar <- ggplot2::ggplot(data.frame(pos = 1:data_dim, methods::slot(BC_dat, "reads"), stringsAsFactors = TRUE), ggplot2::aes_string(x = "pos",
                                                                        y = "read_count")) +
    ggplot2::geom_bar(stat = "identity", fill = "darkgray") +
    ggplot2::xlab("barcodes") + ggplot2::ylab("barcode reads\n") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  if (loga) {
    ggbar <- ggbar + ggplot2::scale_y_continuous(trans = "log2") + ggplot2::ylab("log2(barcode reads)\n")
  }

  if (!is.null(ori_BCs)) {

    if (col_type == "darkgray") {
      col_type <- "rainbow"
    }

    if (method == "hamming") {
      minDist <- .getMinDist(BC_dat, ori_BCs, method)
      dist_col <- .generateColors(minDist, type = col_type)

      ggbar <- ggbar + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = factor(minDist))) +
        ggplot2::scale_fill_manual("HDs", values = dist_col[[2]]) +
        ggplot2::labs(fill = "HDs")
    } else {
      distance <- .getMinDist(BC_dat, ori_BCs, method)
      dist_col <- .generateColors(distance, type = col_type)

      ggbar <- ggbar + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = factor(distance))) +
        ggplot2::scale_fill_manual("known BCs", values = dist_col[[2]]) +
        ggplot2::labs(fill = method)
    }
  }
  return(ggbar)
}

generateKirchenplot_separate <- function(BC_dat, ori_BCs1 = NULL, ori_BCs2 = NULL, loga = TRUE,
                                         col_type = NULL, setLabels = c("BC-Set 1", "Rest", "BC-Set 2"),
                                         method = "hamming") {

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(ggplot2::ggplot())
  }

  if (!is.null(ori_BCs1)) {

    ind1 <- methods::slot(BC_dat, "reads")$barcode %in% ori_BCs1
    ind2 <- methods::slot(BC_dat, "reads")$barcode %in% ori_BCs2

    type <- rep(NA, length(ind1))
    type[ind1] <- setLabels[1]
    type[ind2] <- setLabels[3]
    type[!(ind1 | ind2)] <- setLabels[2]

    if (method == "hamming") {
      minDist <- .getMinDist(BC_dat, c(ori_BCs1, ori_BCs2), method)
      dist_col <- .generateColors(minDist, type = ifelse(is.null(col_type), "rainbow", col_type))

      names(dist_col[[2]]) <- unique(sort(as.numeric(minDist)))

      data_dim <- dim(methods::slot(BC_dat, "reads"))[1]

      res <- data.frame(pos = 1:data_dim, methods::slot(BC_dat, "reads"), type = factor(type, levels = c(setLabels[1], setLabels[2], setLabels[3])), minDists = factor(minDist, levels = unique(sort(as.numeric(names(dist_col[[2]]))))), colrs = dist_col[[1]], stringsAsFactors = TRUE)

      ggbar <- ggplot2::ggplot(res, ggplot2::aes_string(x = "pos", y = "read_count", fill = "minDists")) +
        ggplot2::geom_bar(stat = "identity")  +
        ggplot2::facet_grid(type ~ .) +
        ggplot2::xlab("barcodes") + ggplot2::ylab("barcode reads\n") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank()) +
        #ggplot2::scale_y_continuous(trans = "log2") +
        #ggplot2::ylab("log2(barcode reads)\n") +
        ggplot2::scale_fill_manual(values = dist_col[[2]], name = paste0("HDs"))

    } else {
        distance <- .getMinDist(BC_dat, ori_BCs1, method)
        dist_col <- .generateColors(distance, type = ifelse(is.null(col_type), "rainbow", col_type))

        data_dim <- dim(methods::slot(BC_dat, "reads"))[1]

        res <- data.frame(pos = 1:data_dim, methods::slot(BC_dat, "reads"), type = factor(type, levels = c(setLabels[1], "contaminations", setLabels[2])), distances = factor(distance, levels = unique(sort(as.numeric(names(dist_col[[2]]))))), colrs = dist_col[[1]], stringsAsFactors = TRUE)

        ggbar <- ggplot2::ggplot(res, ggplot2::aes_string(x = "pos", y = "read_count", fill = "distances")) +
          ggplot2::geom_bar(stat = "identity")  +
          ggplot2::facet_grid(type ~ .) +
          ggplot2::xlab("barcodes") + ggplot2::ylab("barcode reads\n") +
          ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank()) +
          ggplot2::scale_fill_manual(values = dist_col[[2]], name = paste0("distance"))
    }
  }

  if (loga) {
    ggbar <- ggbar + ggplot2::scale_y_continuous(trans = "log2") + ggplot2::ylab("log2(barcode reads)\n")
  }

    return(ggbar)
  }


#' @title Plotting a Barplot
#'
#' @description Generates a barplot visualising the abundances of unique read count frequencies.
#'
#' @param BC_dat a BCdat object.
#' @param b an integer value, defining the number of bins. Overridden by bw. Defaults to 30. (see `?ggplot2::geom_histogram`)
#' @param bw an integer value, defining the width of the bins.
#' @param show_it a logical vaue. If TRUE, the respective values are printed on the console?
#' @param log a logical vaue. If TRUE, the y-axis will be on a log scale.
#' @param dens a logical vaue. If TRUE, the density of the read frequencies will be plotted.
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data(BC_dat)
#' plotReadFrequencies <- function(BC_dat, b = 10, show_it = TRUE)
plotReadFrequencies <- function(BC_dat, b = 30, bw = NULL, show_it = FALSE, log = FALSE, dens = FALSE) {

  if (dens & log) {
    log <- FALSE
    warning("Its either density = TRUE or log = TRUE.")
  }

  if (dim(methods::slot(BC_dat, "reads"))[1] == 0) {
    return(ggplot2::ggplot())
  }

  if (show_it) {
    tmp <- table(methods::slot(BC_dat, "reads")$read_count)
    print(data.frame('read count' = names(tmp), 'number of BCs' = as.numeric(tmp)), stringsAsFactors = TRUE)
  }

  plot_data <- methods::slot(BC_dat, "reads")
  if (log) {
    plot_data[, 1] <- plot_data[, 1] + 1
  }

  if (dens) {
    ..density.. <- NULL
    ggbar <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "read_count")) + #, colour = "black"
      ggplot2::geom_histogram(ggplot2::aes(y = ..density..), binwidth = 1, fill = "gray") +
      ggplot2::geom_density(alpha = .2, fill = "#66CCFF") +
      ggplot2::xlab("read count") +
      ggplot2::ylab("barcode read frequencies\n")
  }  else {
    ggbar <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "read_count")) +
      ggplot2::xlab("read count") + # , colour = "black"
      ggplot2::geom_histogram(bins = b, binwidth = bw, fill = "gray") +
      ggplot2::ylab(ifelse(log, "log(number of barcodes)\n", "number of barcodes\n"))
    if (log) {
      ggbar <- ggbar + ggplot2::scale_y_log10()
    }
  }

  return(ggbar + ggplot2::theme_light())
}

#' @title Plotting Quality Score Distribution
#' @description Creates a plot of the quality values accommodated by the fastq file.
#'
#' @param source_dir a character string of the path to the source directory.
#' @param file_name a character string of the file name.
#' @param type a character string, possible values are "mean" and "median".
#' @param rel a logical value. If TRUE the y-axis will show relative frequency instead of the absolut counts.
#'
#' @return  a ggplot2 object.
#' @export
#' @examples
#'
#' \dontrun{
#'
#' source_dir <- system.file("extdata", package = "genBaRcode")
#'
#' plotQualityScoreDis(source_dir, file_name = "test_data.fastq", type = "mean")
#'
#' }
plotQualityScoreDis <- function(source_dir, file_name, type = "median", rel = FALSE) {

  dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)

  w <- ShortRead::width(dat)[1]
  scores <- methods::as(ShortRead::FastqQuality(Biostrings::quality(attr(dat, "quality"))), "matrix")
  scores <- data.frame(mean = rowSums(scores)/w, median = apply(scores, 1, stats::median), stringsAsFactors = TRUE)

  gghist <- ggplot2::ggplot(scores, ggplot2::aes_string(x = type)) +
    ggplot2::xlab(paste(type, "score per sequence")) +
    ggplot2::theme_light()

  if (rel) {
    ..count.. <- NULL
    gghist +
      ggplot2::geom_histogram(binwidth = 1, ggplot2::aes(y = (..count../sum(..count..)))) +
      ggplot2::ylab("frequencies")
  } else {
    gghist + ggplot2::geom_histogram(binwidth = 1)
  }
}

#' @title Plotting Nucleotide Frequency
#'
#' @description Creates a plot visualising the nucleotide frequency within the entire fastq file.
#'
#' @param source_dir a character string containing the path to the sequencing file.
#' @param file_name a character string containng the name of the sequencing file.
#'
#' @return a ggplot2 object.
#' @export

plotNucFrequency <- function(source_dir, file_name) {

  qa_summary <- ShortRead::qa(dirPath = paste(source_dir, "/", file_name, sep = ""))

  dat <- data.frame(pos = 1:5,
                    label = names(qa_summary[["baseCalls"]]),
                    values = as.numeric(qa_summary[["baseCalls"]] / sum(qa_summary[["baseCalls"]])),
                    stringsAsFactors = TRUE)

  ggplot2::ggplot(dat, ggplot2::aes_string(x = "pos", y = "values")) + ggplot2::geom_bar(stat = "identity") +
    ggplot2::ylab("frequency") + ggplot2::xlab("") +
    ggplot2::scale_x_continuous(breaks = 1:5, labels = names(qa_summary[["baseCalls"]])) + ggplot2::theme_light()

}

#' @title Plotting Quality Score per Cycle
#'
#' @description Visualises the mean, median, 25% and 75% quantile of the quality score per sequencing cycle.
#'
#' @param source_dir a character string containing the path to the sequencing file.
#' @param file_name a character string containng the name of the sequencing file.
#'
#' @return a ggplot2 object.
#' @export

plotQualityScorePerCycle <- function(source_dir, file_name) {

  qa_summary <- ShortRead::qa(dirPath = paste(source_dir, "/", file_name, sep = ""))
  dat <- qa_summary[["perCycle"]]$quality

  mean_dat <- quant_dat1 <- quant_dat2 <- quant_dat3 <- NULL
  for(i in unique(dat$Cycle)) {
    index <- dat$Cycle == i
    mean_dat <- c(mean_dat, sum(dat$Score[index] * dat$Count[index]) / sum(dat$Count[index]))
    quant_dat1 <- c(quant_dat1, stats::quantile(S4Vectors::Rle(dat$Score[index], dat$Count[index]), 0.25))
    quant_dat2 <- c(quant_dat2, stats::quantile(S4Vectors::Rle(dat$Score[index], dat$Count[index]), 0.5))
    quant_dat3 <- c(quant_dat3, stats::quantile(S4Vectors::Rle(dat$Score[index], dat$Count[index]), 0.75))
  }

  tmp <- data.frame(Cycle = unique(dat$Cycle),
                    Quality1 = mean_dat,
                    Quality2 = quant_dat1,
                    Quality3 = quant_dat2,
                    Quality4 = quant_dat3, stringsAsFactors = TRUE)

  tmp <- reshape2::melt(tmp, id = "Cycle")

  colrs <- RColorBrewer::brewer.pal(11, "RdBu")[c(2, 8, 10, 9)]

  ggplot2::ggplot(tmp, ggplot2::aes_string(x = "Cycle", y = "value", linetype = "variable", color = "variable")) +
        ggplot2::geom_line() +
        ggplot2::ylab("quality score") + ggplot2::xlab("cycle") +
        ggplot2::scale_color_manual(name  = "",
                         breaks = c("Quality1", "Quality2", "Quality3", "Quality4"),
                         labels = c("mean", "25% quantile", "median", "75% quantile"),
                         values = colrs) +
        ggplot2::scale_linetype_manual(name = "",
                                       breaks = c("Quality1", "Quality2", "Quality3", "Quality4"),
                                       labels = c("mean", "25% quantile", "median", "75% quantile"),
                                       values = c("solid", "dotted", "solid", "dotdash")) + ggplot2::theme_light()

}

#' Plots a sequence logo
#'
#' @param BC_dat a chatacter vector or BCdat object containing the respective sequences
#' @param colrs a character vector containing the desired colors for the nucleotides A, T, C, G and N (in that order)
#'
#' @return a ggplot2 object
#' @export

plotSeqLogo <- function(BC_dat, colrs = NULL) {

  if (sum(methods::is(BC_dat) == "BCdat") == 1) {
    BC_dat <- as.character(methods::slot(BC_dat, "reads")$barcode)
  }
  if (length(BC_dat) == 0) {
    return(ggplot2::ggplot())
  }

  if (length(colrs) != 5) {
    if (!is.null(colrs)) {
      warning(" # 5 colors needed (A, T, C, G, N)")
    }
    colrs <- c(RColorBrewer::brewer.pal(5, "Set3")[-2], "darkgray")
  }

  freqs <- unlist(lapply(BC_dat, function(x) {
    unlist(strsplit(as.character(x), split = ""))
  }))
  freqs <- matrix(freqs, byrow = TRUE, ncol = nchar(as.character(BC_dat[1])))
  freqs <- apply(freqs, 2, function(x) {
              table(factor(x, levels = c("A", "T", "C", "G", "N")))
            })

  freqdf <- as.data.frame(t(freqs))
  freqdf <- freqdf/length(BC_dat)
  freqdf$pos = as.numeric(as.character(rownames(freqdf)))

  lmf <- reshape2::melt(freqdf, id.var = "pos")
  colrs <- colrs[as.numeric(lmf$variable)]

  ggplot2::ggplot(data = lmf, ggplot2::aes_string(x = "pos", y = "value", fill = "variable")) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.75, color = "white") +
    ggplot2::scale_fill_manual(values = colrs) +
    ggplot2::theme_bw() +
    ggplot2::xlab("position") + ggplot2::ylab("percent") +
    ggplot2::theme(legend.title = ggplot2::element_blank())

}

# plotSeqLogo <- function(BC_dat, colrs = NULL) {
#
#   if (sum(methods::is(BC_dat) == "BCdat") == 1) {
#     BC_dat <- as.character(methods::slot(BC_dat, "reads")$barcode)
#   }
#   if (length(BC_dat) == 0) {
#     return(ggplot2::ggplot())
#   }
#
#   if (length(colrs) != 5) {
#     if (!is.null(colrs)) {
#       warning(" # 5 colors needed (A, T, C, G, N)")
#     }
#     colrs <- c(RColorBrewer::brewer.pal(5, "Set3")[-2], "black")
#   }
#
#   cs <- .make_col_scheme(chars = c('A', 'T', 'C', 'G', 'N'), cols = colrs)
#
#   ggplot2::ggplot() +
#     .geom_logo(data = BC_dat, namespace = 'AGTCN', method = "p", col_scheme = cs) +
#     ggplot2::theme_light() +
#     ggplot2::theme(#axis.ticks.y = ggplot2::element_line(size = 0),
#                    panel.grid.major = ggplot2::element_blank(),
#                    panel.grid.minor = ggplot2::element_blank()) +
#     ggplot2::ylab("probability")
# }

