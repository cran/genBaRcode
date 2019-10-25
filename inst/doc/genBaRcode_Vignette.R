## ----setup, include = FALSE----------------------------------------------

#rmarkdown::html_vignette

  knitr::opts_knit$set(
              self.contained = TRUE)

  knitr::opts_chunk$set(
    #collapse = TRUE,
    dpi = 55,
    fig.retina = 1,
    comment = "#>"
  )

  require("genBaRcode")
  require("ggplot2")
  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    if (!requireNamespace("BiocManager", quietly = TRUE)) {
#        install.packages("BiocManager")
#    }
#  
#    BiocManager::install(c("Biostrings", "ShortRead", "S4Vectors", "ggtree"))
#  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    require("genBaRcode")
#  
#    bb <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
#    source_dir <- system.file("extdata", package = "genBaRcode")
#  
#    BC_data <- processingRawData(file_name = "test_data.fastq.gz",
#                                source_dir = source_dir,
#                                results_dir = "/my/results/directory/",
#                                mismatch = 0,
#                                label = "test",
#                                bc_backbone = bb,
#                                bc_backbone_label = "BC_1",
#                                min_score = 30,
#                                min_reads = 2,
#                                save_it = FALSE,
#                                seqLogo = FALSE,
#                                cpus = 1,
#                                strategy = "sequential",
#                                full_output = FALSE,
#                                wobble_extraction = TRUE,
#                                dist_measure = "hamming")
#  

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  getBackboneSelection()
  
  bb <- getBackboneSelection(1)
  show(bb)
  
  bb <- getBackboneSelection("BC32-eBFP")
  show(bb)
  

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  bb <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
  source_dir <- system.file("extdata", package = "genBaRcode")

  # if no results_dir is provided the source_dir automatically also becomes the results_dir
  BC_data <- processingRawData(file_name = "test_data.fastq.gz", 
                              source_dir = source_dir, 
                              mismatch = 0, 
                              label = "test", 
                              bc_backbone = bb, 
                              bc_backbone_label = "BC_1", 
                              min_score = 30, 
                              min_reads = 2, 
                              save_it = FALSE, 
                              seqLogo = FALSE, 
                              cpus = 1, 
                              strategy = "sequential", 
                              full_output = FALSE, 
                              wobble_extraction = TRUE, 
                              dist_measure = "hamming")
  

## ----echo = FALSE, eval=TRUE, collapse=TRUE------------------------------

  methods::slot(BC_data, "results_dir") <- "/my/results/dir/"


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  # if no results_dir is provided the source_dir automatically also becomes the results_dir
  BC_data_multiple <- processingRawData(file_name = "test_data.fastq.gz", 
                              source_dir = source_dir,
                              mismatch = 0, 
                              label = "test", 
                              bc_backbone = getBackboneSelection(1:2), 
                              bc_backbone_label = c("BC_1", "BC_2"), 
                              min_score = 30, 
                              min_reads = 2, 
                              save_it = FALSE, 
                              seqLogo = FALSE, 
                              cpus = 1, 
                              strategy = "sequential", 
                              full_output = FALSE, 
                              wobble_extraction = FALSE, 
                              dist_measure = "hamming")
  

## ----echo = FALSE, eval=TRUE, collapse=TRUE------------------------------

  methods::slot(BC_data_multiple[[1]], "results_dir") <- "/my/results/dir/"
  methods::slot(BC_data_multiple[[2]], "results_dir") <- "/my/results/dir/"


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_multiple)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  # if no results_dir is provided the source_dir automatically also becomes the results_dir
  BC_data_2 <- processingRawData(file_name = "test_data.fastq.gz",
                              source_dir = source_dir,
                              mismatch = 4,
                              label = "test",
                              bc_backbone = "none",
                              min_score = 30,
                              min_reads = 2,
                              save_it = FALSE,
                              seqLogo = FALSE,
                              cpus = 1,
                              strategy = "sequential",
                              full_output = FALSE,
                              wobble_extraction = FALSE,
                              dist_measure = "hamming")
  

## ----echo = FALSE, eval=TRUE, collapse=TRUE------------------------------

  methods::slot(BC_data_2, "results_dir") <- "/my/results/dir/"


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_2)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  head(getReads(BC_data))

  show(getResultsDir(BC_data))
  
  show(getBackbone(BC_data))
  
  show(getLabel(BC_data))


## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    BC_data <- setReads(BC_data, data.frame(read_count = c(1:5), barcode = letters[1:5]))
#    BC_data <- setResultsDir(BC_data, "/my/test/folder/")
#    BC_data <- setBackbone(BC_data, "AAANNNNGGG")
#    BC_data <- setLabel(BC_data, "new label")
#  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    BC_data <- readBCdat(path = "/my/test/folder/",
#                         label = "test",
#                         BC_backbone = "AAANNNNCCCC",
#                         file_name = "test.csv",
#                         s = ";")
#  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    BC_data_EC <- errorCorrection(BC_dat = BC_data,
#                                 maxDist = 4,
#                                 save_it = FALSE,
#                                 cpus = 1,
#                                 strategy = "sequential",
#                                 m = "hamming",
#                                 type = "standard",
#                                 only_EC_BCs = TRUE,
#                                 EC_analysis = FALSE,
#                                 start_small = TRUE)
#  

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  BC_data_EC <- errorCorrection(BC_dat = BC_data, 
                               maxDist = 4,
                               save_it = FALSE, 
                               cpus = 1,
                               strategy = "sequential", 
                               m = "hamming", 
                               type = "standard",
                               only_EC_BCs = TRUE, 
                               EC_analysis = FALSE, 
                               start_small = TRUE)

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_EC)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  BC_data_EC <- errorCorrection(BC_dat = BC_data, 
                               maxDist = 4,
                               save_it = FALSE, 
                               cpus = 1,
                               strategy = "sequential", 
                               m = "hamming", 
                               type = "standard",
                               only_EC_BCs = TRUE, 
                               EC_analysis = FALSE, 
                               start_small = FALSE)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_EC)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  BC_data_EC <- errorCorrection(BC_dat = BC_data, 
                               maxDist = 4,
                               save_it = FALSE, 
                               cpus = 1,
                               strategy = "sequential", 
                               m = "hamming", 
                               type = "graph based",
                               only_EC_BCs = TRUE, 
                               EC_analysis = FALSE, 
                               start_small = FALSE)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_EC)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  BC_data_EC <- errorCorrection(BC_dat = BC_data, 
                               maxDist = 4,
                               save_it = FALSE, 
                               cpus = 1,
                               strategy = "sequential", 
                               m = "hamming", 
                               type = "connectivity based",
                               only_EC_BCs = TRUE, 
                               EC_analysis = FALSE, 
                               start_small = FALSE)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_EC)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  BC_data_EC <- errorCorrection(BC_dat = BC_data, 
                               maxDist = 4,
                               save_it = FALSE, 
                               cpus = 1,
                               strategy = "sequential", 
                               m = "hamming", 
                               type = "clustering",
                               only_EC_BCs = TRUE, 
                               EC_analysis = FALSE, 
                               start_small = FALSE)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(BC_data_EC)


## ----eval=TRUE, fig.width=2.5, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  s_dir <- system.file("extdata", package = "genBaRcode")

  plotNucFrequency(source_dir = s_dir, file_name = "test_data.fastq.gz")


## ----eval=TRUE, fig.height=1.5, fig.width=5, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  plotQualityScoreDis(source_dir = s_dir, file_name = "test_data.fastq.gz", type = "mean")


## ----eval=TRUE, fig.height=1.5, fig.width=5, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  plotQualityScoreDis(source_dir = s_dir, file_name = "test_data.fastq.gz", type = "median")


## ----eval=TRUE, fig.width=6, fig.height=4, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  plotQualityScorePerCycle(source_dir = s_dir, file_name = "test_data.fastq.gz")


## ----eval=TRUE, fig.width=6.5, fig.height=1.5, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  show(BC_data)

  plotSeqLogo(BC_dat = BC_data, colrs = NULL)


## ----eval=TRUE, fig.width=6.5, fig.height=1.5, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  # color order correlates to the following nucleotide order A, T, C, G, N
  col_vec <- c("#000000", 
               "#000000", 
               RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)], 
               "#000000")
  show(col_vec)

  plotSeqLogo(BC_dat = BC_data, colrs = col_vec)


## ----eval=TRUE, fig.width=6, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----
  
  show(BC_data)
  
  generateKirchenplot(BC_dat = BC_data)


## ----eval=FALSE, fig.width=6, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----
#  
#    known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT",
#                   "CACGATCCGCTTCTATCGCGTGCACTACATGT",
#                   "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")
#  
#    generateKirchenplot(BC_dat = BC_data, ori_BCs = known_BCs)
#  

## ----echo=FALSE, eval=TRUE, fig.width=6, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT", 
                 "CACGATCCGCTTCTATCGCGTGCACTACATGT", 
                 "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")
  
  generateKirchenplot(BC_dat = BC_data, ori_BCs = known_BCs) + ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                                                                              legend.key.size = ggplot2::unit(4, "mm"),
                                                                              legend.title = ggplot2::element_text(size = 7))


## ----eval=TRUE, fig.width=7.2, fig.height=4, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT", 
                 "CACGATCCGCTTCTATCGCGTGCACTACATGT", 
                 "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")
  contaminations <- c("CACGATCCGCTTCTATCGCGTGCACTACATGC", 
                      "ATTGGGTCCGTCTGAGGGCGTCTCTGCGCCTT", 
                      "CACGATCCGCTTCTATCGCGTGCGCTACATGT", 
                      "TACGATCCGCTTCTATCGCGTGCACTACATGT")
    
  generateKirchenplot(BC_dat = BC_data, ori_BCs = known_BCs, ori_BCs2 = contaminations)


## ----eval=FALSE, fig.width=7.2, fig.height=4, fig.align="center",  fig.cap="Figure 5.4: Extracted barcodes and their abundancies.", collapse=TRUE----
#  
#    known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT",
#                   "CACGATCCGCTTCTATCGCGTGCACTACATGT",
#                   "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")
#    contaminations <- c("CACGATCCGCTTCTATCGCGTGCACTACATGC",
#                        "ATTGGGTCCGTCTGAGGGCGTCTCTGCGCCTT",
#                        "CACGATCCGCTTCTATCGCGTGCGCTACATGT",
#                        "TACGATCCGCTTCTATCGCGTGCACTACATGT")
#  
#    generateKirchenplot(BC_dat = BC_data,
#                        ori_BCs = known_BCs, ori_BCs2 = contaminations,
#                        setLabels = c("known BCs", "stuff", "contaminations"),
#                        loga = TRUE, col_type = "wild", m = "lv")
#  

## ----eval=TRUE, fig.width=2.5, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  plotReadFrequencies(BC_dat = BC_data)


## ----eval=FALSE, fig.width=2.5, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----
#  
#    plotReadFrequencies(BC_dat = BC_data, log = TRUE)
#    plotReadFrequencies(BC_dat = BC_data, dens = TRUE)
#  

## ----eval=FALSE, fig.width=2.5, fig.height=2, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----
#  
#    plotReadFrequencies(BC_dat = BC_data, bw = 30)
#    plotReadFrequencies(BC_dat = BC_data, b = 30)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#    plotDistanceVisNetwork(BC_dat = BC_data, minDist = 1, loga = TRUE, m = "hamming")
#    plotDistanceIgraph(BC_dat = BC_data, minDist = 1, loga = TRUE, m = "hamming")
#  

## ----eval=TRUE, fig.width=3, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  ggplotDistanceGraph(BC_dat = BC_data, minDist = 1, loga = TRUE, m = "hamming")


## ----eval=TRUE, fig.width=4.5, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='asis', collapse=TRUE----

  known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT", 
                 "CACGATCCGCTTCTATCGCGTGCACTACATGT", 
                 "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")

  ggplotDistanceGraph(BC_dat = BC_data, 
                      minDist = 1, loga = TRUE, m = "hamming", 
                      ori_BCs = known_BCs, lay = "circle", complete = FALSE, 
                      col_type = "topo.colors", legend_size = 2)


## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    createGDF(BC_dat = BC_data, minDist = 1, loga = TRUE, m = "hamming")
#  

## ----eval=TRUE, fig.width=4, fig.height=4, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  plotClusterTree(BC_dat = BC_data, tree_est = "UPGMA", 
                  type = "fan", tipLabel = FALSE, m = "hamming")

## ----eval=TRUE, fig.width=3, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  plotClusterGgTree(BC_dat = BC_data, tree_est = "NJ", 
                    type = "rectangular", m = "hamming")


## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    BC_data_EC <- errorCorrection(BC_dat = BC_data,
#                                   maxDist = 4,
#                                   save_it = FALSE,
#                                   cpus = 1,
#                                   strategy = "sequential",
#                                   m = "hamming",
#                                   type = "standard",
#                                   only_EC_BCs = FALSE,
#                                   EC_analysis = TRUE,
#                                   start_small = FALSE)
#  
#    error_correction_clustered_HDs(datEC = BC_data_EC, size = 0.75)
#  

## ----echo=FALSE, fig.width=2, fig.height=2.5, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  BC_data_EC <- errorCorrection(BC_dat = BC_data, 
                                 maxDist = 4,
                                 save_it = FALSE, 
                                 cpus = 1,
                                 strategy = "sequential", 
                                 m = "hamming", 
                                 type = "standard",
                                 only_EC_BCs = FALSE, 
                                 EC_analysis = TRUE, 
                                 start_small = FALSE)
  
  error_correction_clustered_HDs(datEC = BC_data_EC, size = 0.75) + ggplot2::theme(axis.title = ggplot2::element_text(size = 8))
  

## ----eval=TRUE, fig.width=4, fig.height=4, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  error_correction_circlePlot(edges = BC_data_EC$edges, vertices = BC_data_EC$vertices)


## ----eval=TRUE, fig.width=3, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  error_correction_treePlot(edges = BC_data_EC$edges, vertices = BC_data_EC$vertices)
    

## ----eval=TRUE, fig.width=4, fig.height=4, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  ggplotDistanceGraph_EC(BC_dat = BC_data, BC_dat_EC = BC_data_EC, 
                         minDist = 1, loga = TRUE, m = "hamming")
  

## ----eval=FALSE, fig.width=3, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----
#  
#    plotDistanceVisNetwork_EC(BC_dat = BC_data, BC_dat_EC = BC_data_EC,
#                              minDist = 1, loga = TRUE, m = "hamming")
#  

## ----eval=TRUE, fig.width=3, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT", 
                 "CACGATCCGCTTCTATCGCGTGCACTACATGT", 
                 "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")

  ggplotDistanceGraph_EC(BC_dat = BC_data, BC_dat_EC = BC_data_EC, 
                         minDist = 1, loga = TRUE, m = "hamming", ori_BCs = known_BCs)
  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    plotDistanceVisNetwork_EC(BC_dat = BC_data, BC_dat_EC = BC_data_EC,
#                              minDist = 1, loga = TRUE, m = "hamming", ori_BCs = known_BCs)
#  

## ----eval=TRUE, fig.width=3, fig.height=3, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  known_BCs <- c("GGTCGAAGCTTCTTTCGGGCCGCACGGCTGCT", 
                 "CACGATCCGCTTCTATCGCGTGCACTACATGT", 
                 "ATTGGGTCCGTCTGAGGGCGTTTCTGCGCCTT")

  ggplotDistanceGraph_EC(BC_dat = BC_data, BC_dat_EC = BC_data_EC, 
                         minDist = 1, loga = TRUE, m = "hamming", BC_threshold = 2)
  

## ----eval=TRUE, , fig.width=3, fig.height=2.5, fig.pos = 'H', fig.align='center', fig.show='hold', collapse=TRUE----

  # path to the package internal data file
  source_dir <- system.file("extdata", package = "genBaRcode")

  BC_data_tp1 <- processingRawData(file_name = "test_data.fastq.gz", 
                              source_dir, 
                              mismatch = 10, 
                              label = "tp1", 
                              bc_backbone = getBackboneSelection(1), 
                              bc_backbone_label = "BC_1", 
                              min_score = 10, 
                              save_it = FALSE)
  BC_data_tp1 <- errorCorrection(BC_data_tp1, maxDist = 2)

  BC_data_tp2 <- processingRawData(file_name = "test_data.fastq.gz", 
                              source_dir, 
                              mismatch = 1, 
                              label = "tp2", 
                              bc_backbone = getBackboneSelection(1), 
                              bc_backbone_label = "BC_1", 
                              min_score = 30,
                              min_reads = 1000,
                              save_it = FALSE)
  BC_data_tp2 <- errorCorrection(BC_data_tp2, maxDist = 4, type = "clustering")

  BC_data_tp3 <- processingRawData(file_name = "test_data.fastq.gz", 
                              source_dir, 
                              mismatch = 0, 
                              label = "tp3", 
                              bc_backbone = getBackboneSelection(1), 
                              bc_backbone_label = "BC_1", 
                              min_score = 37, 
                              save_it = FALSE)
  BC_data_tp3 <- errorCorrection(BC_data_tp3, maxDist = 8, type = "graph based")
  
  BC_list <- list(BC_data_tp1, BC_data_tp2, BC_data_tp3)
  BC_matrix <- generateTimeSeriesData(BC_dat_list = BC_list)
  plotTimeSeries(ov_dat = BC_matrix)
  plotVennDiagram(BC_dat = BC_list)


## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    # choose colors
#    test_colors <- RColorBrewer::brewer.pal(12, "Set3")
#  
#    plotTimeSeries(ov_dat = BC_matrix[1:12, ],
#                   colr = test_colors, tp = c(1,3,4),
#                   x_label = "test data", y_label = "test freqs")
#  
#    plotVennDiagram(BC_dat = BC_list, alpha_value = 0.25,
#                    colrs = c("green", "red", "blue"), border_color = "orange",
#                    plot_title = "this is the title",
#                    legend_sort = c("tp2_EC", "tp3_EC", "tp1_EC"),
#                    annotationSize = 2.5)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#    # start Shiny app with the package internal test data file
#    genBaRcode_app()
#  
#    # start Shiny app with access to a predefined directory
#    genBaRcode_app(dat_dir = "/my/test/directory/")
#  

## ----eval=TRUE, out.width = 40, collapse=TRUE----------------------------

  getBackboneSelection()
  
  bb <- getBackboneSelection(1)
  show(bb)
  
  bb <- getBackboneSelection("BC32-eBFP")
  show(bb)
  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    BC_data <- readBCdat(path = "/my/test/firectory", label = "test_label", s = ";",
#                         BC_backbone = "ACTNNGGCNNTGANN", file_name = "test_file.csv")
#  

## ----eval=FALSE, collapse=TRUE-------------------------------------------
#  
#    test_data_frame <- data.frame(read_count = seq(100, 400, 100),
#                                  barcode = c("AAAAAAAA", "GGGGGGGG",
#                                              "TTTTTTTT", "CCCCCCCC"))
#  
#    BC_data <- asBCdat(dat = test_data_frame,
#                        label = "test_label",
#                        BC_backbone = "CCCNNAAANNTTTNNGGGNN",
#                        resDir = "/my/results/directory/")
#  

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  test_data_frame <- data.frame(read_count = seq(100, 400, 100), 
                                barcode = c("AAAAAAAA", "GGGGGGGG", 
                                            "TTTTTTTT", "CCCCCCCC"))

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(test_data_frame)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  BC_data_1 <- asBCdat(dat = test_data_frame, 
                      label = "test_label_1", 
                      BC_backbone = "CCCNNAAANNTTTNNGGGNN", 
                      resDir = getwd())

  test_data_frame <- data.frame(read_count = c(300, 99, 150, 400), 
                                barcode = c("TTTTTTTT", "AATTTAAA", 
                                            "GGGGGGGG", "CCCCCCCC"))

## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(test_data_frame)                              


## ----eval=TRUE, collapse=TRUE--------------------------------------------
  BC_data_2 <- asBCdat(dat = test_data_frame, 
                      label = "test_label_2", 
                      BC_backbone = "CCCNNAAANNTTTNNGGGNN", 
                      resDir = getwd())

  test <- genBaRcode:::com_pair(BC_dat1 = BC_data_1, BC_dat2 = BC_data_2)


## ----eval=TRUE, collapse=TRUE--------------------------------------------

  show(test)
  

