
makeReactiveTrigger <- function() {
  rv <- reactiveValues(a = 0)
  list(
    depend = function() {
      rv$a
    },
    trigger = function() {
      rv$a <- isolate(rv$a + 1)
    },
    reset = function() {
      rv$a <- 0
    }
  )
}

makeHistory <- function() {
  hist <- reactiveValues(old = NULL)

  list(
    old = function() {
      hist$old
    },
    update = function(x) {
      hist$old <- x
    }
  )
}

makeGeneralHistory <- function() {
  hist <- reactiveValues(show = NULL, show_sCode = NULL)

  list(
    show = function() {
      hist$show
    },
    update = function(x) {
      hist$show <- x
    },
    show_sCode = function() {
      hist$show_sCode
    },
    update_sCode = function(y) {
      hist$show_sCode <- y
    }
  )
}

shinyServer(
  function(input, output, session) {

    givenDir <- genBaRcode:::.testDirIdentifier(getOption("genBaRcode-shinyDir"))
    flag <- givenDir != paste0(system.file("extdata", package = "genBaRcode"), .Platform$file.sep)

    fileTrigger <- makeReactiveTrigger()
    history <- makeHistory()
    genHist <- makeGeneralHistory()

    ############################
    ### UI functional linkage
    ############################

    output$selection <- renderUI({
        if (fileTrigger$depend() == 0) {
          processing_UI()
        } else {
          options("genBaRcode-info" = "")
          plot_UI()
        }
    })

    output$G_and_T <- renderUI({
      if (fileTrigger$depend() == 0) {
        NULL
      } else {
        G_and_T_UI()
      }
    })

    output$caption <- renderText({
      capText()
    })

    ##################
    ### caption
    ##################

    capText <- reactive({
      if(fileTrigger$depend() == 0 | is.null(input$plot_click)) {
        return(getOption("genBaRcode-info"))
      } else {

          ind <- NULL

          if(input$plot_type == "Kirchenplot") {
              ind <- which(((ggplot2::ggplot_build(finalPlot())$data[[1]]$xmin <= input$plot_click$x) + (ggplot2::ggplot_build(finalPlot())$data[[1]]$xmax >= input$plot_click$x)) == 2)
          } else {
            if(input$plot_type == "HD Graph") {
              ind_x <- abs(ggplot2::ggplot_build(finalPlot())$data[[2]]$x - input$plot_click$x)
              ind_y <- abs(ggplot2::ggplot_build(finalPlot())$data[[2]]$y - input$plot_click$y)
              ind <- which.min(ind_x + ind_y)
            }
          }
          if(!is.null(ind)) {
            cap <- paste("- ", slot(dat(), "reads")$barcode[ind], ": ", slot(dat(), "reads")$read_count[ind], " -", sep = "")
            isolate(genHist$update(rbind(c(as.character(slot(dat(), "reads")$barcode[ind]),
                                              as.character(slot(dat(), "reads")$read_count[ind])), genHist$show())))
            return(cap)
          } else {
            return("")
          }

      }
    })

    #############################
    ### plot creation
    #############################

    dat <- reactive({
      if(length(input$f_name) == 1) {
        if(input$error_corr & input$plot_type != "HD Graph" & input$plot_type != "interactive HD Graph") {
          isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "ec", gDir = givenDir))))
          tmp <- BC()$datEC
          tmp@reads <- slot(BC()$datEC, "reads")[1:input$maxBCs_EC, ]
        } else {
          tmp <- BC()$dat
          tmp@reads <- slot(BC()$dat, "reads")[1:input$maxBCs, ]
        }
      } else {
        if(input$error_corr) {
          tmp <- BC()$datEC[1:input$maxBCs_EC, ]
        } else {
          tmp <- BC()$dat[1:input$maxBCs, ]
        }

      }
        return(tmp)
    })

    finalPlot2 <- reactive({
      if(input$plot_type == "interactive HD Graph" & fileTrigger$depend() != 0) {
        return(plotHammDistVisNetwork(dat(), minHD = 1, loga = input$loga, oriBCs(), complete = input$compl, col_type = input$palette))
      }
    })

    finalPlot <-  reactive({

      if(fileTrigger$depend() == 0) {
        return(NULL)
      }

      if(input$plot_type == "Kirchenplot") {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "kirchenplot", gDir = givenDir))))
        return(generateKirchenplotBCdis(dat(), ori_BCs = oriBCs(), loga = input$loga, col_type = input$palette))
      }

      if(input$plot_type == "Read Dist.") {
        return(generateKirchenplotRCdis(dat(), b = input$bins, show_it = FALSE))
      }

      # if(input$plot_type == "SeqLogo") {
      #   return(plotSeqLogo(as.character(slot(dat(), "reads")$"barcode")))
      # }

      if(input$plot_type == "Tree") {
        if(input$tree_est == "Neighbor-Joining") { tEst <- "NJ" }
        if(input$tree_est == "Unweighted Pair Group Method (UPGMA)") { tEst <- "UPGMA" }

        return(plotClusterTree(dat(), type = input$tree_style, tree_est = tEst, tipLabel = input$show_tip_label))
      }

      if(input$plot_type == "HD Graph") {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "HD_Graph", gDir = givenDir))))
        return(ggplotHammDistGraph(dat(), minHD = 1, loga = input$loga, oriBCs(), lay = input$graph_layout, complete = input$compl, col_type = input$palette))
      }

      if(length(input$f_name) > 1) {
          return(plotTimeSeries(dat()))
      }
    })

    ####################
    ### table creation
    ####################

    output$table_overview <- renderDataTable({
      if(length(input$f_name) == 1) {
        d <- data.frame(feature = c("number of barcodes",
                                    "number of barcodes (EC)",
                                    "read count min",
                                    "read count median",
                                    "read count mean",
                                    "read count max",
                                    "read count min (EC)",
                                    "read count median (EC)",
                                    "read count mean (EC)",
                                    "read count max (EC)"),
                        data = c(dim(slot(BC()$dat, "reads"))[1],
                                 dim(slot(BC()$datEC, "reads"))[1],
                                 summary(slot(BC()$dat, "reads")[, 1])[c(1, 3, 4, 6)],
                                 summary(slot(BC()$datEC, "reads")[, 1])[c(1, 3, 4, 6)])
        )
      }
      if(length(input$f_name) > 1) {
        d <- data.frame(feature = c("number of barcodes",
                                    "number of barcodes (EC)",
                                    paste("time point", 1:dim(BC()$dat)[2], "(reads)")),
                        data = c(dim(BC()$dat)[1],
                                 dim(BC()$datEC)[1],
                                 colSums(BC()$dat))
        )
      }
      if((length(input$f_name) == 0) | is.null(input$f_name)) {
        d <- NULL
      }

      return(d)
    })

    output$table_dat <- renderDataTable({

      if(length(input$f_name) == 1) {
        return(slot(BC()$dat, "reads")[, 2:3])
      }
      if(length(input$f_name) > 1) {
        tmp <- data.frame(BC = row.names(BC()$dat))
        tmp <- cbind(tmp, as.data.frame(BC()$dat))
        colnames(tmp) <- c("BC", paste("tp", 1:dim(BC()$dat)[2]))
        return(tmp)
      }
    })

    output$table_EC <- renderDataTable({

      if(length(input$f_name) == 1) {
        return(slot(BC()$datEC, "reads")[, 2:3])
      }
      if(length(input$f_name) > 1) {
        tmp <- data.frame(BC = row.names(BC()$datEC))
        tmp <- cbind(tmp, as.data.frame(BC()$datEC))
        colnames(tmp) <- c("BC", paste("tp", 1:dim(BC()$datEC)[2]))
        return(tmp)
      }
    })

    output$table_clickHist <- renderDataTable({

      if(length(input$f_name) == 1 & !is.null(genHist$show())) {
        tmp <-  genHist$show()
        colnames(tmp) <- c("barcode", "read counts")
        return(tmp)
      }
      if(length(input$f_name) > 1) {
        return("")
      }

    })

    output$table_sCodeHist <- renderDataTable({
        tmp <- cbind(genHist$show_sCode(), "")
        colnames(tmp) <- c("", "")
        return(tmp)
    })

    ####################
    ### UI design
    ####################

    processing_UI <- function() {
      span(
        br(),
        selectInput("f_name", "choose file:", choices = c("", list.files(givenDir, pattern = ".fast")), multiple = TRUE),
        fileInput("file1", "choose known BCs file", accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv")),

        selectInput("backbone", "choose backbone:", choices = c("", "BC16",
                                                                    "BC32-GFP",
                                                                    "BC32-Venus",
                                                                    "BC32-eBFP",
                                                                    "BC32-T-Sapphire",
                                                                    "own Design")),

        conditionalPanel(
          condition = "input.backbone == 'own Design'",
          textInput("ownDesign", label = "backbone", value = "", width = NULL, placeholder = "N indicates variable positions e.g. ACTNNGCANN")
        ),

        fluidRow(
          column(8,
                 numericInput("mm", "mismatches:", 0, min = 0)),
          column(8,
                 numericInput("minReads", "min. reads:", 3, min = 0)),
          column(8,
                 numericInput("maxHD", "maxHD EC", 1, min = 0))
        ),
        checkboxInput("quaFil", "quality filtering", TRUE),

        actionButton(inputId = "go", label = "Go"),
        actionButton(inputId = "exit", label = "Exit")
      )
    }

    plot_UI <- function() {

      if(length(input$f_name) == 1) {
        cs <- c("", "Kirchenplot", "Read Dist.", "HD Graph", "interactive HD Graph", "Tree") #"SeqLogo",
      } else {
        cs <- "Time Series"
      }

      span(
        fluidRow(
          column(8,
                 selectInput("plot_type", "plot type:",
                             choices = cs)
          ),
          column(4,
                 conditionalPanel(
                   condition = "(!(input.error_corr)) | input.plot_type == 'HD Graph' | input.plot_type == 'interactive HD Graph'",
                   numericInput("maxBCs", "max. BCs:",
                                value = ifelse(length(input$f_name) == 1, dim(slot(BC()$dat, "reads"))[1], dim(BC()$dat)[1]),
                                min = 0,
                                max = ifelse(length(input$f_name) == 1, dim(slot(BC()$dat, "reads"))[1], dim(BC()$dat)[1])
                   )
                 ),
                 conditionalPanel(
                   condition = "input.error_corr & input.plot_type != 'HD Graph' & input.plot_type != 'interactive HD Graph'",
                   numericInput("maxBCs_EC", "max. BCs:",
                                value = ifelse(length(input$f_name) == 1, dim(slot(BC()$datEC, "reads"))[1], dim(BC()$datEC)[1]),
                                min = 0,
                                max = ifelse(length(input$f_name) == 1, dim(slot(BC()$datEC, "reads"))[1], dim(BC()$datEC)[1])
                   )
                 )
           )
        ),
        conditionalPanel(
          condition = "input.plot_type == 'Read Dist.'",
          numericInput("bins", "bins:", value = 30, min = 0, max = 100)
        ),

        conditionalPanel(
          condition = "input.plot_type == 'HD Graph'",
          selectInput("graph_layout", "graph layout:",
                      choices = c("fruchtermanreingold", "kamadakawai", "adj"))
        ),

        conditionalPanel(
          condition = "input.plot_type == 'Tree'",
          selectInput("tree_est", "tree estimation alg:",
                      choices = c("Neighbor-Joining", "Unweighted Pair Group Method (UPGMA)")),
          selectInput("tree_style", "tree layout:",
                      choices = c("unrooted", "phylogram", "cladogram", "fan", "radial"))
        ),

        conditionalPanel(
          condition = "input.plot_type == 'Kirchenplot' | input.plot_type == 'HD Graph' | input.plot_type == 'interactive HD Graph'",
          selectInput("palette", "color palette:",
                      choices = c("rainbow", "heat", "topo.colors"))
        ),

        conditionalPanel(
          condition = "input.plot_type != '' & input.plot_type != 'Read Dist.' & input.plot_type != 'HD Graph' & input.plot_type != 'interactive HD Graph'",
          checkboxInput("error_corr", "error correction", FALSE)
        ),

        conditionalPanel(
          condition = "input.plot_type != 'SeqLogo' & input.plot_type != 'Tree' & input.plot_type != '' & input.plot_type != 'Read Dist.'",
          checkboxInput("loga", "log values", FALSE)
        ),
        conditionalPanel(
          condition = "input.plot_type == 'HD Graph' | input.plot_type == 'interactive HD Graph'",
          checkboxInput("compl", "complete graph", FALSE)
        ),
        conditionalPanel(
          condition = "input.plot_type == 'Tree'",
          checkboxInput("show_tip_label", "show tip label", FALSE)
        ),
        br(),
        actionButton(inputId = "new", label = "New file"),
        actionButton(inputId = "save", label = "Save", style = ifelse(flag, "color: black", "color: lightgray")),
        actionButton(inputId = "exit", label = "Exit")
      )
    }

    G_and_T_UI <- function() {
      if(fileTrigger$depend() != 0) {
        span(
         conditionalPanel(
           condition = "input.plot_type != 'interactive HD Graph'",
           plotOutput("final_plot", width = "100%", click = "plot_click")
         ),
         conditionalPanel(
           condition = "input.plot_type == 'interactive HD Graph'",
           visNetwork::visNetworkOutput("final_plot2", width = "100%")
         ),

          tabsetPanel(
            tabPanel("overview", dataTableOutput("table_overview")),
            tabPanel('barcode list', dataTableOutput("table_dat")),
            tabPanel('barcode list (EC)', dataTableOutput("table_EC")),
            tabPanel('click history', dataTableOutput("table_clickHist")),
            tabPanel('source code', dataTableOutput("table_sCodeHist"))
          )
        )
      }
    }

    observe({
      if(!(is.null(input$plot_type)) & fileTrigger$depend() != 0) {
          if(input$plot_type != "interactive HD Graph") {
            output$final_plot <- renderPlot({
              finalPlot()
            })
          } else {
            output$final_plot2 <- visNetwork::renderVisNetwork({
              finalPlot2()
            })
          }
      }
    })

    #####################
    ### reactive part
    ####################

    oriBCs <- reactive({
      if(fileTrigger$depend() == 0 | identical(isolate(history$old()), input$file1$datapath)) {
         return(NULL)
      } else {
          if(!is.null(input$file1$name)){
            end <- unlist(strsplit(input$file1$name, split = "[.]"))[2]
            if(end == "csv") {
              s <- ";"
            } else {
              if(end == "txt") {
                s <- ""
              } else {
                warning("invalid file format!")
                return(NULL)
              }
            }
            isolate(history$update(input$file1$datapath))
            isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "oriBCs", s = s, gDir = givenDir))))
            return(as.character(unlist(read.table(input$file1$datapath, header = FALSE, sep = s, fill = TRUE))))
          }
      }
    })

    observeEvent(input$new, {
      options("genBaRcode-info" = "")
      fileTrigger$reset()
    })

    observeEvent(input$exit, {
      stopApp()
    })

    observeEvent(input$save, {
      if(flag) {
          if(input$plot_type == "interactive HD Graph") {
            visNetwork::visSave(graph = finalPlot2(), file = paste(givenDir, input$plot_type, "_", slot(BC()$dat, "label"), ".html", sep = ""))
          } else {
            if(input$plot_type == "SeqLogo") {
                grDevices::png(paste(givenDir, input$plot_type, "_", slot(BC()$dat, "label"), ".png", sep = ""))
                  plotSeqLogo(as.character(slot(dat(), "reads")$"barcode"))
                grDevices::dev.off()
            } else {
              if(input$plot_type == "Tree") {
                if(input$tree_est == "Neighbor-Joining") { tEst <- "NJ" }
                if(input$tree_est == "Unweighted Pair Group Method (UPGMA)") { tEst <- "UPGMA" }

                grDevices::png(paste(givenDir, input$plot_type, "_", slot(BC()$dat, "label"), ".png", sep = ""))
                  plotClusterTree(dat(), type = input$tree_style, tree_est = tEst, tipLabel = input$show_tip_label)
                grDevices::dev.off()
              } else {
                    ggplot2::ggsave(filename = paste(input$plot_type, "_", slot(BC()$dat, "label"), ".png", sep = ""),
                           plot = finalPlot(),
                           device = "png",
                           path = givenDir)
              }
            }
          }
      }
    })

    observeEvent(input$go, {
      req(input$f_name, input$mm, input$minReads)
      fileTrigger$trigger()
      BC()
    })

    BC <- reactive({

      if(input$backbone == "BC32-GFP") {
        bcp <- "ACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANNCTTNNCGANNCTTNNGGANNCTANNACTNNCGANN"
      }
      if(input$backbone == "BC32-Venus") {
        bcp <- "CGANNAGANNCTTNNCGANNCTANNGGANNCTTNNCGANNAGANNCTTNNCGANNCTANNGGANNCTTNNCGANNAGANN"
      }
      if(input$backbone == "BC32-eBFP") {
        bcp <- "CTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNN"
      }
      if(input$backbone == "BC32-T-Sapphire") {
        bcp <- "CAGNNATCNNCTTNNCGANNGGANNCTANNCTTNNCAGNNATCNNCTTNNCGANNGGANNCTANNCTTNNCAGNNATCNN"
      }
      if(input$backbone == "BC16") {
            bcp <- "AGATCNNTAGNNTCCNNAAGNNTCGNNAAGNNTCGNNAGTNNTAGAT"
      }
      if(input$ownDesign != "") {
              bcp <- input$ownDesign
      }

      if(length(input$f_name) == 1) {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "read", bcp = bcp, gDir = givenDir))))
        tmp <- get_dat_single(bcp, givenDir)
        if(sum(dim(slot(tmp, "reads"))) == 0 & fileTrigger$depend() > 0) {
          options("genBaRcode-info" = "...no barcode sequences detectable...")
          fileTrigger$reset()
        } else {
            list(
              dat = tmp,
              datEC = errorCorrection(tmp, maxHD = input$maxHD, save_it = flag)
            )
        }
      } else {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "timeS", bcp = bcp, gDir = givenDir))))
        list(
          dat = get_timeS(bcp, givenDir),
          datEC = get_timeS_EC(bcp, givenDir)
        )
      }
    })

    ##################
    # basic functions
    ####################

    get_dat_single <- function(bcp, gDir) {

      processingRawData(file_name = input$f_name, source_dir = gDir,
                        results_dir = gDir,
                        mismatch = input$mm,
                        bc_pattern = bcp,
                        quality_filtering = input$quaFil, min_score = 30, min_reads = input$minReads, unix = FALSE,
                        save_it = flag)

    }

    get_timeS <- function(bcp, gDir) {
      dat_list <- vector("list", length(input$f_name))
      for(i in 1:length(input$f_name)) {

        dat_list[[i]] <- processingRawData(file_name = input$f_name[i], source_dir = gDir,
                                           results_dir = gDir,
                                           mismatch = input$mm,
                                           bc_pattern = bcp,
                                           quality_filtering = input$quaFil, min_score = 30, min_reads = input$minReads, unix = FALSE,
                                           save_it = flag)
      }
      return(generateTimeSeriesData(dat_list))
    }

    get_timeS_EC <- function(bcp, gDir) {
      dat_list <- vector("list", length(input$f_name))
      for(i in 1:length(input$f_name)) {

        dat_list[[i]] <- processingRawData(file_name = input$f_name[i], source_dir = gDir,
                                           results_dir = gDir,
                                           mismatch = input$mm,
                                           bc_pattern = bcp,
                                           quality_filtering = input$quaFil, min_score = 30, min_reads = input$minReads, unix = FALSE,
                                           save_it = flag)
        if(input$maxHD > 0) {
          dat_list[[i]] <- errorCorrection(dat_list[[i]], maxHD = input$maxHD, save_it = flag)
        }
      }
      return(generateTimeSeriesData(dat_list))
    }

    sCode_snippets <- function(type, gDir = "", bcp = "", oriBCs = NULL, s = "") {

      if(type == "timeS") {
        return(
          matrix(
              c(
                paste("nameList <- c(", paste("'", input$f_name, collapse =", ", "'", sep = ""), ")"),
                "dat_list <- vector('list', length(nameList))",
                "for(i in 1:length(nameList)) {",
                paste("dat_list[[i]] <- processingRawData(file_name = nameList[i], source_dir = '", gDir, "', results_dir = '", gDir,
            "', mismatch = ", input$mm, ", bc_pattern = '", bcp, "', quality_filtering = ", input$quaFil, ", min_score = 30, min_reads = ", input$minReads, ", unix = FALSE)"),
            "}",
            "dat <- generateTimeSeriesData(dat_list)",
            "plotTimeSeries(dat)"
          ), ncol = 1, nrow = 7)
        )
      } else {
            return(
              switch(type,
                     read = isolate(paste("BC_dat <- processingRawData(file_name = '", input$f_name, "', source_dir = '", gDir,
                                          "', results_dir = '", gDir, "', mismatch = ", input$mm, ", bc_pattern = '", bcp,
                                          "', quality_filtering = ", input$quaFil, ", min_score = 30, min_reads = ", input$minReads, ", unix = FALSE)", sep = "")),
                     ec = isolate(paste("BC_dat_EC <- errorCorrection(BC_dat, maxHD = ", input$maxHD, ")", sep = "")),
                     oriBCs = paste("oriBCs <- as.character(unlist(read.table('", input$file1$datapath, "', header = FALSE, sep = ", s, ", fill = TRUE)))", sep = ""),
                     kirchenplot = isolate(paste("generateKirchenplotBCdis(BC_dat, ", ifelse(identical(isolate(history$old()), input$file1$datapath), " ", "oriBCs = oriBCs,"), " loga = ", input$loga, ", col_type = ", input$palette, ")", sep = "")),
                     HD_Graph = isolate(paste("ggplotHammDistGraph(BC_dat, minHD = 1, ", ifelse(identical(isolate(history$old()), input$file1$datapath), " ", "oriBCs = oriBCs,"), " loga = ", input$loga, ", lay = ", input$graph_layout, ", complete = ", input$compl, ", col_type = ", input$palette, ")", sep = "")),
                     tree = isolate(paste("plotClusterTree(BC_dat, type = ", input$tree_style, ", tipLabel = ", input$show_tip_label, ")", sep = ""))
              )
            )
      }

  }
})


