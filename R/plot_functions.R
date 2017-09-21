
#' @title Plotting a Hamming Distance Network
#' @description plotHammDistVisNetwork will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the ggplot2 and the ggnetwork packages. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a hamming distance between two barcodes
#' of maximal \code{minHD}.
#' If \code{ori_BCs} is provided the node color also refelects the hamming distance of a particular barcode to one of the given
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param minHD an integer value representing the maximal hamming distance value for which the graph will
#' contain edges.
#' @param loga a logical value indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#'
#' @return a visNetwork object.
#' @export
#'
#' @examples
#' data(BC_dat)
#' plotHammDistVisNetwork(BC_dat, minHD = 1, loga = TRUE, ori_BCs = NULL,
#' complete = TRUE, col_type = "rainbow")

plotHammDistVisNetwork <- function(BC_dat, minHD = 1, loga = TRUE, ori_BCs = NULL, complete = TRUE, col_type = "rainbow") {

  net <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = "hamming")
  if(complete) {
    net <- t(apply(net, 1, function(x) {
      x[x > min(x[x > 0])] <- 0
      return(x)
    }))
    weight <- net
    net[net > 0] <- 1
  } else {
    net[net > minHD] <- 0
  }

  net[.getDiagonalIndex(dim(net)[1])] <- 0

  net <- igraph::graph.adjacency(net)

  if(loga) {
    v <- log(methods::slot(BC_dat, "reads")$"read_count")
  } else {
    v <- methods::slot(BC_dat, "reads")$"read_count"
  }

  if(length(ori_BCs) > 0) {

    minHDs <- .getMinHammDist(BC_dat, ori_BCs)
    colrs <- .generateColors(minHDs, type = col_type)

    nodes <- data.frame(id = 1:length(igraph::V(net)),
                        label = methods::slot(BC_dat, "reads")$"barcode",
                        title = methods::slot(BC_dat, "reads")$"barcode",
                        value = v,
                        shadow = TRUE,
                        color.background = as.character(colrs[[1]]),
                        color.border = "red",
                        color.highlight.border <- "gray")
  } else {
    nodes <- data.frame(id = 1:length(igraph::V(net)),
                        label = methods::slot(BC_dat, "reads")$"barcode",
                        title = methods::slot(BC_dat, "reads")$"barcode",
                        value = v,
                        shadow = TRUE,
                        color.background = "black",
                        color.border = "darkgray",
                        color.highlight.border <- "gray")
  }


  edges <- data.frame(from = igraph::get.edgelist(net)[, 1], to = igraph::get.edgelist(net)[, 2])

  graph_visNetwork <- visNetwork::visNetwork(nodes, edges, width = "100%")
  graph_visNetwork <- visNetwork::visOptions(graph_visNetwork, highlightNearest = TRUE)
  graph_visNetwork <- visNetwork::visIgraphLayout(graph_visNetwork, "layout_nicely")
  graph_visNetwork <- visNetwork::visHierarchicalLayout(graph_visNetwork, graph_visNetwork)

  return(graph_visNetwork)
}

#' @title Plotting a Hamming Distance Network
#'
#' @description ggplotHammDistGraph will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the ggplot2 and the ggnetwork packages. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a hamming distance between two barcodes
#' of maximal \code{minHD}.
#' If \code{ori_BCs} is provided the node color also refelects the hamming distance of a particular barcode to one of the initial
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param minHD an integer value representing the maximal hamming distance for which the graph will
#' contain edges.
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param lay a character string, identifying the prefered layout algorithm (see ggnetwork layout option).
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#' @param outline an integer value which adjusts the thickness of the black outline of each node.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' data(BC_dat)
#' ggplotHammDistGraph(BC_dat, minHD = 1, loga = TRUE, ori_BCs = NULL, lay = "fruchtermanreingold",
#' complete = FALSE, col_type = "rainbow")
#'
#' }
#'

ggplotHammDistGraph <- function(BC_dat, minHD = 1, loga = TRUE, ori_BCs = NULL, lay = "fruchtermanreingold", complete = TRUE, col_type = "rainbow", outline = 0.5) {

  net <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = "hamming")
  if(complete) {
    net <- t(apply(net, 1, function(x) {
                x[x > min(x[x > 0])] <- 0
                return(x)
              }))
    weight <- net
    net[net > 0] <- 1
  } else {
    net[net > minHD] <- 0
  }

  if(sum(rowSums(net)) == 0) {
    return(graphics::plot(-1, -1, graphics::text(1, 1, "No edges left", cex = 2), xlim = c(0, 2), ylim = c(0, 2), frame.plot = FALSE, axes = FALSE, xlab = "", ylab = ""))
  }

  net <- network::network(net, directed = FALSE)

  if(loga) {
    v_size <- log(methods::slot(BC_dat, "reads")$"read_count")
  } else {
    v_size <- methods::slot(BC_dat, "reads")$"read_count"
  }

  network::set.vertex.attribute(net, "size", v_size)
  network::set.vertex.attribute(net, "names", as.character(methods::slot(BC_dat, "reads")$"barcode"))
  if(!is.null(ori_BCs)) {
    minHDs <- .getMinHammDist(BC_dat, ori_BCs)
    colrs <-  .generateColors(minHDs, type = col_type)

    network::set.vertex.attribute(net, "minHDs", minHDs)
    network::set.vertex.attribute(net, "colors", colrs[[1]])
    tmp <- TRUE
  } else {
    network::set.vertex.attribute(net, "minHDs", NA)
    tmp <- FALSE
  }

  net <- ggnetwork::ggnetwork(net, layout = lay, cell.jitter = 0.75)

  p <- ggplot2::ggplot(net, ggplot2::aes(x = net$x, y = net$y, xend = net$xend, yend = net$yend)) +
          ggnetwork::geom_edges(color = "grey50", size = 0.3) +
          ggnetwork::geom_nodes(ggplot2::aes(size = net$size)) +
          ggnetwork::theme_blank() +
          ggplot2::scale_size_area(ifelse(loga, "log(read count)", "reads")) +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4)))
  if(tmp) {
    p <- p + ggplot2::scale_colour_manual(name = "min HDs", values = colrs[[2]]) +
          ggnetwork::geom_nodes(ggplot2::aes(color = as.factor(minHDs), size = (net$size-(net$size*outline))))
  }

  return(p)
}

#' @title Plotting a Hamming Distance Network
#'
#' @description plotHammDistIgraph will create a graph-like visualisation (ripple plot) of the corresponding barcode sequences
#' and their similarity based on the igraph package. The nodes represent the barcode sequences and their
#' respective size reflects the corresponding read counts. Edges between nodes indicate a hamming distance between two barcodes
#' of maximal \code{minHD}.
#' If \code{ori_BCs} is provided the node color also refelects the hamming distance of a particular barcode to one of the initial
#' barcodes.
#'
#' @param BC_dat a BCdat object.
#' @param minHD an integer value representing the maximal hamming distance value for which the graph will
#' contain edges.
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param threeD a logical value to chose between 2D and 3D visualisation.
#' @param complete a logical value. If TRUE, every node will have at least one edge.
#' @param col_type a character sting, choosing one of the available color palettes.
#'
#' @return an igraph object.
#' @export
#'
#' @examples
#' data(BC_dat)
#' plotHammDistIgraph(BC_dat, minHD = 1, loga = TRUE, ori_BCs, threeD = FALSE,
#' complete = TRUE, col_type = "rainbow")


plotHammDistIgraph <- function(BC_dat, minHD = 1, loga = TRUE, ori_BCs, threeD = FALSE, complete = TRUE, col_type = "rainbow") {

  adj_dat <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = "hamming")

  if(complete) {
    adj_dat <- t(apply(adj_dat, 1, function(x) {
      x[x > min(x[x > 0])] <- 0
      return(x)
    }))
    adj_dat[adj_dat > 0] <- 1
  } else {
    adj_dat[adj_dat > minHD] <- 0
  }

  colnames(adj_dat) <- methods::slot(BC_dat, "reads")$"barcode"

  g <- igraph::graph.adjacency(adj_dat, mode = c("undirected"), diag = FALSE, add.colnames = "names")

  igraph::V(g)$size <- log(methods::slot(BC_dat, "reads")$"read_count")
  igraph::V(g)$label <- ""

  if(length(ori_BCs) > 0) {

    minHDs <- .getMinHammDist(BC_dat, ori_BCs)
    colrs <- .generateColors(minHDs, type = col_type)

    igraph::V(g)$color <- colrs[[1]]
  } else {
    igraph::V(g)$color <- "gray"
  }

  if(!threeD) {
    #igraph::plot.igraph(g, layout = igraph::layout_nicely, edge.width=1+(1/igraph::E(g)$weight)) # layout_with_fr
    igraph::plot.igraph(g, layout = igraph::layout_nicely) # layout_with_fr
    graphics::legend("topleft", legend = names(colrs[[2]]), fill = colrs[[2]], bty = "n", title = "HDs")
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
#' @param bar_width a numeric value specifying the visual space between two plotted measurements.
#' @param labs a character vector containing sample labels.
#' @param x_label a character string providing the x-axis label.
#' @param y_label a character string providing the y-axis label.
#' @return a ggplot2 object.
#'
#' @export
#' @examples
#' ov_dat <- matrix(round(runif(1:100, min = 0, max = 1000)), ncol = 5)
#' rownames(ov_dat) <- paste("barcode", 1:20)
#' plotTimeSeries(ov_dat)

plotTimeSeries <- function(ov_dat, colr = NULL, tp = NULL, bar_width = 0.05, labs = NULL, x_label = "time", y_label = "contribution") {

  if(!is.matrix(ov_dat)) {
    stop("# plotTimeSeries requires a numeric matrix")
  }

  if(sum(ov_dat[, 1]) != 1) {
    message("# normalization of read count data")
    ov_dat <- t(t(ov_dat) / colSums(ov_dat))
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
    bar_pos <- tp
  } else {
    bar_pos <- 1:dim(ov_dat)[2]
  }

  x <- y <- NULL
  for(i in 1:(length(bar_pos) - 1) ) {
    for(j in 0:(dim(ov_dat)[1] - 1)) {
      if(j==0){
        x <- c(x,
                  bar_pos[i+1]-bar_width/2,
                  bar_pos[i]+bar_width/2,
                  bar_pos[i]+bar_width/2,
                  bar_pos[i+1]-bar_width/2
                )
        y <- c(y, 0.0,
                  0.0,
                  sum(ov_dat[1:(j+1),i]),
                  sum(ov_dat[1:(j+1),(i+1)])
                )

      } else {
        x <- c(x, bar_pos[i+1]-bar_width/2,
                  bar_pos[i]+bar_width/2,
                  bar_pos[i]+bar_width/2,
                  bar_pos[i+1]-bar_width/2
              )
        y <- c(y, sum(ov_dat[1:j,(i+1)]),
                  sum(ov_dat[1:j,i]),
                  sum(ov_dat[1:(j+1),i]),
                  sum(ov_dat[1:(j+1),(i+1)])
                )
      }
    }
  }

  plot_dat <- data.frame(x = x,
                         y = y,
                         n = rep(1:(dim(ov_dat)[1] * (dim(ov_dat)[2] - 1)), each = 4),
                         id = rep(unlist(lapply(1:(dim(ov_dat)[1]), function(m) rep(m, times = 4))), dim(ov_dat)[2] - 1))

  coolr <- rep(colr, dim(ov_dat)[1] - 1)

  ggOV <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = plot_dat$x, y = plot_dat$y)) +
    ggplot2::geom_polygon(ggplot2::aes(group = plot_dat$n, col = as.factor(plot_dat$id), fill = as.factor(plot_dat$id))) +
    ggplot2::scale_fill_manual(values = coolr) +
    ggplot2::scale_color_manual(values = coolr) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab(y_label)

  if(is.null(labs)) {
    ggOV <- ggOV + ggplot2::scale_x_discrete(name = x_label, limits = bar_pos)
  } else {
    ggOV <- ggOV + ggplot2::scale_x_discrete(name = x_label, limits = bar_pos, labels = labs)
  }

  return(ggOV)
}


#' @title Creating a Gephi File
#' @description createGephiFile creates a data file usable with the free graph visualisation tool gephi. The nodes
#' represent barcodes and its respective size reflects the corresponding read counts. Edges between nodes indicate
#' a hamming distance between two barcodes of maximal \code{minHD}.
#' If \code{ori_BCs} is provided the node color refelects the hamming distance of a particular barcode to one
#' of the provided barcode sequences.
#'
#' @param BC_dat a BCdat object.
#' @param minHD an integer value representing the maximal hamming distance value for which the graph will
#' contain edges.
#' @param loga a logical value indicating the use or non-use of logarithmic read count values.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param col_type character sting, choosing one of the available color palettes.
#'
#' @return NULL
#'
#' @export
#' @examples
#'
#' \dontrun{
#'
#' data(BC_dat)
#' createGephiFile(BC_dat, minHD = 1, loga = TRUE, ori_BCs = NULL, col_type = "rainbow")
#'
#' }

createGephiFile <- function(BC_dat, minHD = 1, loga = TRUE, ori_BCs = NULL, col_type = "rainbow") {

  file_name <- paste(methods::slot(BC_dat, "results_dir"), methods::slot(BC_dat, "label"), ".gdf", sep = "")

  adj_dat <- stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = "hamming")
  adj_dat[adj_dat > minHD] <- 0

  colnames(adj_dat) <- methods::slot(BC_dat, "reads")$"barcode"
  g <- igraph::graph.adjacency(adj_dat, mode=c("undirected"), diag=FALSE, add.colnames = "names")

  if(loga) {
    tmp_reads <- log(as.numeric(methods::slot(BC_dat, "reads")$"read_count"))
  } else {
    tmp_reads <- as.numeric(methods::slot(BC_dat, "reads")$"read_count")
  }

  final <- rbind(c("nodedef>name VARCHAR","label VARCHAR","size DOUBLE"),
                 cbind( paste("node_",1:dim(methods::slot(BC_dat, "reads"))[1],sep=""),
                        as.character(methods::slot(BC_dat, "reads")[, 3]),
                        tmp_reads))

  if(length(ori_BCs) > 0) {
    minHDs <- .getMinHammDist(BC_dat, ori_BCs)
    colrs <- .hex2rgbColor(.generateColors(minHDs, type = col_type)[[1]])
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
#' @description Generates a tree plot based on a herachical clustering of the complete hamming distance matrix.
#'
#' @param BC_dat a BCdat object.
#' @param tree_est a character string, indicating the particular cluster algorithm, possible algorithms are "Neighbor-Joining" ("NJ") and "Unweighted Pair Group Method" ("UPGMA").
#' @param type a character string, indication the graph layout style.
#' @param tipLabel a logical value, indication the use of labeled tree leaves.
#'
#' @return a ggtree object.
#' @export
#' @examples
#' data(BC_dat)
#' plotClusterTree(BC_dat, tree_est = "UPGMA", type = "unrooted", tipLabel = FALSE)

plotClusterTree <- function(BC_dat, tree_est = c("NJ", "UPGMA"), type = c("unrooted", "phylogram", "cladogram", "fan", "radial"), tipLabel = FALSE) {

  hamming_dist <- stats::as.dist(stringdist::stringdistmatrix(methods::slot(BC_dat, "reads")$"barcode", methods::slot(BC_dat, "reads")$"barcode", method = "hamming"))
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

#' @title Plotting a Kirchenplot
#' @description Generates a barplot based on read counts. If \code{ori_BCs} is provided the bar color reflects the hamming
#' distance of a particular barcode to one of the provided barcode sequences.
#'
#' @param BC_dat a BCdat object.
#' @param ori_BCs a vector of character strings containing the barcode sequences (without the fixed positions of the barcode construct).
#' @param loga a logical value, indicating the use or non-use of logarithmic read count values.
#' @param col_type character sting, choosing one of the availabe color palettes, e.g. rainbow.
#'
#' @return a ggplot2 object
#' @export
#' @examples
#' data(BC_dat)
#' generateKirchenplotBCdis(BC_dat, ori_BCs, loga = TRUE, col_type = NULL)

generateKirchenplotBCdis <- function(BC_dat, ori_BCs = NULL, loga = TRUE, col_type = NULL) {

  if(is.null(col_type)) {
    col_type <- "darkgray"
  }

  ggbar <- ggplot2::ggplot(methods::slot(BC_dat, "reads"), ggplot2::aes(x = methods::slot(BC_dat, "reads")$pos,
                                                                        y = methods::slot(BC_dat, "reads")$read_count)) +
                                  ggplot2::geom_bar(stat = "identity", fill = "darkgray") +
                                  ggplot2::xlab("barcodes") + ggplot2::ylab("barcode reads\n") +
                                  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                        axis.ticks.x = ggplot2::element_blank())

  if(loga) {
    ggbar <- ggbar + ggplot2::scale_y_continuous(trans = "log2") + ggplot2::ylab("log2(barcode reads)\n")
  }

  if(!is.null(ori_BCs)) {

    if(col_type == "darkgray") {
              col_type <- "rainbow"
    }

    minHD <- .getMinHammDist(BC_dat, ori_BCs)
    dist_col <- .generateColors(minHD, type = col_type)

    ggbar <- ggbar + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = factor(minHD))) +
            ggplot2::scale_fill_manual("HDs", values = dist_col[[2]]) +
            ggplot2::labs(fill = "HDs")
  }

  return(ggbar)
}


#' @title Plotting a Kirchenplot
#'
#' @description Generates a barplot visualising the read count distribution.
#'
#' @param BC_dat a BCdat object.
#' @param b an integer value, defining the number of bins.
#' @param show_it a logical vaue. If TRUE, the respective values are printed on the console?
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data(BC_dat)
#' generateKirchenplotRCdis <- function(BC_dat, b = 10, show_it = TRUE)

generateKirchenplotRCdis <- function(BC_dat, b = 30, show_it = FALSE) {

  if(show_it) {
    tmp <- table(methods::slot(BC_dat, "reads")$read_count)
    print(data.frame('read count' = names(tmp), 'number of BCs' = as.numeric(tmp)))
  }

  ggbar <- ggplot2::ggplot(methods::slot(BC_dat, "reads"), ggplot2::aes(methods::slot(BC_dat, "reads")$read_count)) +
    ggplot2::geom_histogram(bins = b) +
    ggplot2::xlab("read count") + ggplot2::ylab("number of barcodes\n")

  return(ggbar)
}


#' @title Plotting Quality Score Distribution
#' @description Creates a plot of the quality values accommodated by the fastq file.
#'
#' @param source_dir a character string of the path to the source directory.
#' @param file_name a character string of the file name.
#' @param type a character string, possible values are "mean" and "median".
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
#'

plotQualityScoreDis <- function(source_dir, file_name, type) {

  dat <- ShortRead::readFastq(dirPath = source_dir, pattern = file_name)

  w <- ShortRead::width(dat)[1]
  scores <- methods::as(ShortRead::FastqQuality(Biostrings::quality(attr(dat, "quality"))), "matrix")
  scores <- data.frame(average = rowSums(scores)/w, median = apply(scores, 1, stats::median))

  if(type == "median") {
    gghist <- ggplot2::ggplot(scores, ggplot2::aes(x = scores$median)) +
      ggplot2::xlab("median score per sequence")
  }
  if(type == "mean") {
    gghist <- ggplot2::ggplot(scores, ggplot2::aes(x = scores$average)) +
      ggplot2::xlab("mean score per sequence")
  }

  gghist + ggplot2::geom_histogram(binwidth = 1)

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
                    values = as.numeric(qa_summary[["baseCalls"]] / sum(qa_summary[["baseCalls"]])))

  ggplot2::ggplot(dat, ggplot2::aes(x = dat$pos, y = dat$values)) + ggplot2::geom_bar(stat = "identity") +
    ggplot2::ylab("frequency") + ggplot2::xlab("") +
    ggplot2::scale_x_continuous(breaks = 1:5, labels = names(qa_summary[["baseCalls"]]))

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
#'

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
                    Quality4 = quant_dat3)

  tmp <- reshape2::melt(tmp, id = "Cycle")

  colrs <- RColorBrewer::brewer.pal(11, "RdBu")[c(2, 8, 10, 9)]

  ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$Cycle, y = tmp$value, linetype = tmp$variable, color = tmp$variable)) +
        ggplot2::geom_line() +
        ggplot2::ylab("Quality Score") +
        ggplot2::scale_color_manual(name  = "",
                         breaks = c("Quality1", "Quality2", "Quality3", "Quality4"),
                         labels = c("mean", "25% quantile", "median", "75% quantile"),
                         values = colrs) +
        ggplot2::scale_linetype_manual(name = "",
                                       breaks = c("Quality1", "Quality2", "Quality3", "Quality4"),
                                       labels = c("mean", "25% quantile", "median", "75% quantile"),
                                       values = c("solid", "dotted", "solid", "dotdash"))

}
