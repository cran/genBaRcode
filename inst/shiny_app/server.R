
## FV1.0

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

    if(is.null(getOption("genBaRcode-shinyDir"))) {
      stopApp("Please use the 'genBaRcode_app' function to start the app!")
    }

    seqL <- NULL

    givenDir <- genBaRcode:::.testDirIdentifier(getOption("genBaRcode-shinyDir"))
    flag <- givenDir != paste0(system.file("extdata", package = "genBaRcode"), .Platform$file.sep)

    fileTrigger <- makeReactiveTrigger()
    history <- makeHistory()
    genHist <- makeGeneralHistory()

    ############################
    ### UI functional linkage
    ############################

    # output$selection <- renderUI({
    #     if (fileTrigger$depend() == 0) {
    #       processing_UI()
    #     } else {
    #       options("genBaRcode-info" = "")
    #       plot_UI()
    #     }
    # })

    output$selection <- renderUI({
      if(fileTrigger$depend() == 0) {
        processing_UI_choose_file()
      } else {
        options("genBaRcode-info" = "")
        plot_UI()
      }
    })

    output$parameters <- renderUI({
        if(fileTrigger$depend() == 0 & !is.null(input$fileType)) {
          if(input$fileType == "fastq") {
            processing_UI_fastq()
          } else {
            if(input$fileType == "csv") {
              processing_UI_csv()
            } else {
              NULL
            }
          }
        } else {
          NULL
        }
    })

    output$end <- renderUI({
      if(fileTrigger$depend() == 0) {
        processing_UI_end()
      } else {
        NULL
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
      } # else {
    #
    #       ind <- NULL
    #
    #       if(input$plot_type == "Kirchenplot") {
    #           ind <- which(((ggplot2::ggplot_build(finalPlot())$data[[1]]$xmin <= input$plot_click$x) + (ggplot2::ggplot_build(finalPlot())$data[[1]]$xmax >= input$plot_click$x)) == 2)
    #       } else {
    #         if(input$plot_type == "HD Graph") {
    #           ind_x <- abs(ggplot2::ggplot_build(finalPlot())$data[[2]]$x - input$plot_click$x)
    #           ind_y <- abs(ggplot2::ggplot_build(finalPlot())$data[[2]]$y - input$plot_click$y)
    #           ind <- which.min(ind_x + ind_y)
    #         } else {
    #           if(input$plot_type == "Read Dist.") {
    #             ind <- cap <- NULL
    #           } else {
    #             if(input$plot_type == "Tree") {
    #               ind <- which(((ggplot2::ggplot_build(finalPlot())$data[[1]]$xmin <= input$plot_click$x) + (ggplot2::ggplot_build(finalPlot())$data[[1]]$xmax >= input$plot_click$x)) == 2)
    #             }
    #           }
    #         }
    #       }
    #       if(!is.null(ind)) {
    #         cap <- paste("- ", slot(dat(), "reads")$barcode[ind], ": ", slot(dat(), "reads")$read_count[ind], " -", sep = "")
    #         isolate(genHist$update(rbind(c(as.character(slot(dat(), "reads")$barcode[ind]),
    #                                           as.character(slot(dat(), "reads")$read_count[ind])), genHist$show())))
    #         return(cap)
    #       } else {
    #         return("")
    #       }
    #
    #   }
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
              tmp <- BC()$datEC#[1:input$maxBCs_EC, ]
            } else {
              tmp <- BC()$dat#[1:input$maxBCs, ]
            }
        }
      #}
      return(tmp)
    })

    finalPlot2 <- reactive({
      if(input$plot_type == "interactive HD Graph" & fileTrigger$depend() != 0) {
        return(plotDistanceVisNetwork(dat(), minDist = 1, loga = TRUE, oriBCs(), complete = input$compl, col_type = input$palette))
      }
    })

    finalPlot <-  reactive({

      if(fileTrigger$depend() == 0) {
        return(ggplot2::ggplot() + ggplot2::theme_minimal())
      }

      if(input$plot_type == "Kirchenplot") {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "kirchenplot", gDir = givenDir))))
        return(generateKirchenplot(dat(), ori_BCs = oriBCs(), loga = input$loga, col_type = input$palette))
      }

      if(input$plot_type == "Read Frequencies") {
        return(plotReadFrequencies(dat(), b = input$bins, show_it = FALSE))
      }

      if(input$plot_type == "SeqLogo - BCs") {
         return(plotSeqLogo(as.character(slot(dat(), "reads")$"barcode")))
      }

      if(input$plot_type == "SeqLogo - NGS reads") {
        if(is.null(seqL)) {

          ending <- strsplit(input$f_name, split="[.]", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
            if(ending == "fasta") {
              seqL <<- ShortRead::readFasta(dirPath = givenDir, pattern = input$f_name)
            } else {
              seqL <<- ShortRead::readFastq(dirPath = givenDir, pattern = input$f_name)
            }
          }

          l <- nchar(as.character(ShortRead::sread(seqL)[1]))
          return(plotSeqLogo(as.character(ShortRead::sread(seqL))) + ggplot2::scale_x_continuous(breaks = c(1, round(l/2), l)))
      }

      if(input$plot_type == "Tree") {
        if(input$tree_est == "Neighbor-Joining") { tEst <- "NJ" }
        if(input$tree_est == "Unweighted Pair Group Method (UPGMA)") { tEst <- "UPGMA" }

        return(plotClusterGgTree(dat(), type = input$tree_style, tree_est = tEst))
      }

      if(input$plot_type == "HD Graph") {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "HD_Graph", gDir = givenDir))))
        return(ggplotDistanceGraph(dat(), minDist = 1, loga = TRUE, oriBCs(), lay = input$graph_layout, complete = input$compl, col_type = input$palette) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
      }

      if(input$plot_type == "Time Series") {
        isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "timeS_plot"))))
          return(plotTimeSeries(dat()[[2]]))
      }

       if(input$plot_type == "Venn Diagram") {
         isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "timeS_Venn"))))
         return(plotVennDiagramm(dat()[[1]]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
       }

    })

    ####################
    ### table creation
    ####################

    output$table_overview <- renderDataTable({
      if(fileTrigger$depend() == 0) {
        return(NULL)
      }

      if((length(input$f_name) == 0) | is.null(input$f_name)) {
        d <- NULL
      }

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
                        data = c(as.character(dim(slot(BC()$dat, "reads"))[1]),
                                 as.character(dim(slot(BC()$datEC, "reads"))[1]),
                                 as.character(summary(slot(BC()$dat, "reads")[, 1])[c(1, 3, 4, 6)]),
                                 as.character(summary(slot(BC()$datEC, "reads")[, 1])[c(1, 3, 4, 6)]))
        )
      }
      if(length(input$f_name) > 1) {

        d <- data.frame(feature = c("number of barcodes",
                                    "number of barcodes (EC)",
                                    paste("time point", 1:dim(BC()$dat[[2]])[2], "(reads)")),
                        data = c(as.character(dim(BC()$dat[[2]])[1]),
                                 as.character(dim(BC()$datEC[[2]])[1]),
                                 as.character(colSums(BC()$dat[[2]])))
        )
      }
      return(d)
    })

    output$table_dat <- renderDataTable({
      if(length(input$f_name) == 1) {
        return(slot(BC()$dat, "reads")[, 2:3])
      }
      if(length(input$f_name) > 1) {
        tmp <- data.frame(BC = row.names(BC()$dat[[2]]))
        tmp <- cbind(tmp, as.data.frame(BC()$dat[[2]]))
        colnames(tmp) <- c("BC", paste("tp", 1:dim(BC()$dat[[2]])[2]))
        return(tmp)
      }
    })

    output$table_EC <- renderDataTable({

      if(length(input$f_name) == 1) {
        return(slot(BC()$datEC, "reads")[, 2:3])
      }
      if(length(input$f_name) > 1) {
        tmp <- data.frame(BC = row.names(BC()$datEC[[2]]))
        tmp <- cbind(tmp, as.data.frame(BC()$datEC[[2]]))
        colnames(tmp) <- c("BC", paste("tp", 1:dim(BC()$datEC[[2]])[2]))
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

    processing_UI_choose_file <- function() {
      span(
        br(),
        radioButtons("fileType", label = "file type", choices = c("csv", "fastq"), selected = "fastq", inline = TRUE)
      )
    }

    processing_UI_csv <- function() {
      span(
        br(),
        selectInput("f_name", "choose file", choices = c("", list.files(givenDir, pattern = input$fileType)), multiple = TRUE),
        fileInput("file1", "choose known BCs file", accept = c("text/csv",
                                                               "text/comma-separated-values,text/plain",
                                                               "csv")),

        conditionalPanel(
          condition = "input.backbone == 'own Design'",
          textInput("ownDesign", label = "backbone", value = "", width = NULL, placeholder = "N indicates variable positions e.g. ACTNNGCANN")
        ),

        fluidRow(
          column(8,
                 numericInput("minReads", "min. reads:", 3, min = 0)),
          column(8,
                 numericInput("maxHD", "maxHD EC", 1, min = 0))
        )
      )
    }

    processing_UI_fastq <- function(fname) {
      span(
        br(),
        selectInput("f_name", "choose file", choices = c("", list.files(givenDir, pattern = input$fileType)), multiple = TRUE),
        fileInput("file1", "choose known BCs file", accept = c(".fastq", ".fasta")),

        selectInput("backbone", "choose backbone", choices = c("", "BC16-GFP",
                                                                    "BC16-Venus",
                                                                    "BC16-mCherry",
                                                                    "BC16-Cerulean",
                                                                    "BC32-GFP",
                                                                    "BC32-Venus",
                                                                    "BC32-eBFP",
                                                                    "BC32-T-Sapphire",
                                                                    "Own Design")),

        conditionalPanel(
          condition = "input.backbone == 'own Design'",
          textInput("ownDesign", label = "backbone", value = "", width = NULL, placeholder = "N indicates variable positions e.g. ACTNNGCANN")
        ),

        fluidRow(
          column(8,
                 numericInput("mm", "mismatches", 0, min = 0)),
          column(8,
                 numericInput("minReads", "min. reads", 3, min = 0)),
          column(8,
                 numericInput("maxHD", "maxHD EC", 1, min = 0))
        ),
        checkboxInput("quaFil", "quality filtering", TRUE)
      )
    }

    processing_UI_end <- function() {
      span(
        br(),
        actionButton(inputId = "go", label = "Go"),
        actionButton(inputId = "exit", label = "Exit")
      )
    }


# Original UI
#
#     processing_UI <- function(fname) {
#       span(
#         br(),
#         selectInput("f_name", "choose file:", choices = c("",
#                                                           list.files(givenDir, pattern = ".fast"),
#                                                           list.files(givenDir, pattern = ".csv")), multiple = TRUE),
#         fileInput("file1", "choose known BCs file", accept = c("text/csv",
#                                                                "text/comma-separated-values,text/plain",
#                                                                ".csv")),
#
#         selectInput("backbone", "choose backbone:", choices = c("", "BC16-GFP",
#                                                                 "BC16-Venus",
#                                                                 "BC16-mCherry",
#                                                                 "BC16-Cerulean",
#                                                                 "BC32-GFP",
#                                                                 "BC32-Venus",
#                                                                 "BC32-eBFP",
#                                                                 "BC32-T-Sapphire",
#                                                                 "Own Design")),
#
#         conditionalPanel(
#           condition = "input.backbone == 'own Design'",
#           textInput("ownDesign", label = "backbone", value = "", width = NULL, placeholder = "N indicates variable positions e.g. ACTNNGCANN")
#         ),
#
#         fluidRow(
#           column(8,
#                  numericInput("mm", "mismatches:", 0, min = 0)),
#           column(8,
#                  numericInput("minReads", "min. reads:", 3, min = 0)),
#           column(8,
#                  numericInput("maxHD", "maxHD EC", 1, min = 0))
#         ),
#         checkboxInput("quaFil", "quality filtering", TRUE),
#
#         actionButton(inputId = "go", label = "Go"),
#         actionButton(inputId = "exit", label = "Exit")
#       )
#     }

    plot_UI <- function() {

      if(length(input$f_name) == 1) {
        cs <- c("", "Kirchenplot", "Read Frequencies", "HD Graph", "interactive HD Graph", "Tree", "SeqLogo - BCs")
        if(input$fileType == "fastq") {
          cs <- c(cs, "SeqLogo - NGS reads")
        }
      } else {
        cs <- c("", "Time Series", "Venn Diagram")
      }

      span(
        fluidRow(
          column(8,
                 selectInput("plot_type", "plot type:",
                             choices = cs)
          ),
          column(4,
                 conditionalPanel(
                   condition = "(((!(input.error_corr)) |
                                   input.plot_type == 'HD Graph' |
                                   input.plot_type == 'interactive HD Graph')) & input.f_name == 1",
                   numericInput("maxBCs", "max. BCs:",
                                value = ifelse(length(input$f_name) == 1, dim(slot(BC()$dat, "reads"))[1], dim(BC()$dat[[2]])[1]),
                                min = 0,
                                max = ifelse(length(input$f_name) == 1, dim(slot(BC()$dat, "reads"))[1], dim(BC()$dat[[2]])[1])
                   )
                 ),
                 conditionalPanel(
                   condition = "(input.error_corr & input.plot_type != 'HD Graph' & input.plot_type != 'interactive HD Graph') & input.f_name == 1",
                   numericInput("maxBCs_EC", "max. BCs:",
                                value = ifelse(length(input$f_name) == 1, dim(slot(BC()$datEC, "reads"))[1], dim(BC()$datEC[[2]])[1]),
                                min = 0,
                                max = ifelse(length(input$f_name) == 1, dim(slot(BC()$datEC, "reads"))[1], dim(BC()$datEC[[2]])[1])
                   )
                 )
           )
        ),
        conditionalPanel(
          condition = "input.plot_type == 'Read Frequencies'",
          numericInput("bins", "bins:", value = 30, min = 0, max = 100)
        ),

        conditionalPanel(
          condition = "input.plot_type == 'HD Graph'",
          selectInput("graph_layout", "graph layout:",
                      choices = c("fruchtermanreingold", "kamadakawai"))
        ),

        conditionalPanel(
          condition = "input.plot_type == 'Tree'",
          selectInput("tree_est", "tree estimation alg:",
                      choices = c("Neighbor-Joining", "Unweighted Pair Group Method (UPGMA)")),
          selectInput("tree_style", "tree layout:",
                      choices = c('rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle', 'daylight'))
        ),
        if(!is.null(oriBCs())) {
            conditionalPanel(
              condition = "input.plot_type == 'Kirchenplot' | input.plot_type == 'HD Graph' | input.plot_type == 'interactive HD Graph'",
              selectInput("palette", "color palette:",
                          choices = c("rainbow", "heat", "topo.colors"))
            )
        },

        conditionalPanel(
          condition = "input.plot_type != 'SeqLogo - NGS' & input.plot_type != '' & input.plot_type != 'Read Frequencies' & input.plot_type != 'HD Graph' & input.plot_type != 'interactive HD Graph'",
          checkboxInput("error_corr", "error correction", FALSE)
        ),

        conditionalPanel(
          condition = "input.plot_type == 'Kirchenplot'",
          checkboxInput("loga", "log values", FALSE)
        ),
        conditionalPanel(
          condition = "input.plot_type == 'HD Graph' | input.plot_type == 'interactive HD Graph'",
          checkboxInput("compl", "complete graph", FALSE)
        ),
        # conditionalPanel(
        #   condition = "input.plot_type == 'Tree'",
        #   checkboxInput("show_tip_label", "show tip label", FALSE)
        # ),
        br(),
        fluidRow(
          column(2, actionButton(inputId = "new", label = "New", width = '55px'), offset = 0),
          conditionalPanel(
                  condition = "input.plot_type == 'interactive HD Graph'",
                  column(2, actionButton(inputId = "save", label = "Save", width = '55px', style = ifelse(flag, "color: black", "color: lightgray")), offset = 0)
          ),
          column(2, actionButton(inputId = "exit2", label = "Exit", width = '55px'), offset = 0)
        )
      )
    }

    G_and_T_UI <- function() {
      if(fileTrigger$depend() != 0) {
        span(
         conditionalPanel(
           condition = "input.plot_type != 'interactive HD Graph'",
           #plotOutput("final_plot", width = "100%", click = "plot_click")
           plotly::plotlyOutput("final_plot", width = "100%")
         ),
         conditionalPanel(
           condition = "input.plot_type == 'interactive HD Graph'",
           visNetwork::visNetworkOutput("final_plot2", width = "100%")
         ),
          br(),
          tabsetPanel(
            tabPanel("overview", dataTableOutput("table_overview")),
            tabPanel('barcode list', dataTableOutput("table_dat")),
            tabPanel('barcode list (EC)', dataTableOutput("table_EC")),
            #tabPanel('click history', dataTableOutput("table_clickHist")),
            tabPanel('source code', dataTableOutput("table_sCodeHist"))
          )
        )
      }
    }

    observe({
      if(!(is.null(input$plot_type)) & fileTrigger$depend() != 0) {
          if(input$plot_type == "interactive HD Graph") {
              output$final_plot2 <- visNetwork::renderVisNetwork({
                finalPlot2()
              })
           } else {
                if(input$plot_type != "") {
                # output$final_plot <- renderPlot({
                #   finalPlot()
                # })
                  output$final_plot <- plotly::renderPlotly({
                        p <- suppressMessages(plotly::ggplotly(finalPlot()))
                        p$elementId <- NULL
                        if(input$plot_type != "SeqLogo - NGS" & input$plot_type != "SeqLogo" & input$plot_type != "Tree" & input$plot_type != "Venn Diagram" & input$plot_type != "Time Series") {
                          p$x$data[[ifelse(input$plot_type == "HD Graph" & length(p$x$data) > 1, 2, 1)]]$text <- paste(
                            paste("BC:", methods::slot(BC_dat, "reads")$barcode),
                            paste("reads:", methods::slot(BC_dat, "reads")$read_count), sep = "<br />")
                          p$x$layout$margin$l <- p$x$layout$margin$l + 15
                        } else {
                          p$x$data[[1]]$text[seq(1, length(methods::slot(BC_dat, "reads")$barcode) * 3, 1)] <- unlist(lapply(
                                                    paste(paste("BC:", methods::slot(BC_dat, "reads")$barcode, "<br />"),
                                                    paste("reads:", methods::slot(BC_dat, "reads")$read_count)), function(x) {
                                                                                                                        c(x, x, NA)
                                                                                                                  }))
                        }
                        p
                    })
                }
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
    observeEvent(input$exit2, {
      stopApp()
    })

    observeEvent(input$save, {
      if(flag) {
          if(input$plot_type == "interactive HD Graph") {
            visNetwork::visSave(graph = finalPlot2(), file = paste(givenDir, input$plot_type, "_", slot(BC()$dat, "label"), ".html", sep = ""))
            }
      }
    })

    observeEvent(input$go, {
      req(input$f_name)
      if(unlist(strsplit(input$f_name, split = "[.]"))[2] != "csv") {
        req(input$mm, input$minReads, input$backbone)
      }
      fileTrigger$trigger()
      BC()
    })

    BC <- reactive({

      bcp <- "not defined"
      cpus <- parallel::detectCores() - floor(parallel::detectCores()/2)

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
      if(input$backbone == "BC16-GFP") {
            bcp <- "AGATCNNTAGNNTCCNNAAGNNTCGNNAAGNNTCGNNAGTNNTAGAT"
      }
      if(input$backbone == "BC16-Venus") {
        bcp <- "CTANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNGAT"
      }
      if(input$backbone == "BC16-mCherry") {
        bcp <- "CTANNCAGNNATCNNCTTNNCGANNGGANNCTANNCTTNNGAT"
      }
      if(input$backbone == "BC16-Cerulean") {
        bcp <- "CTANNCACNNAGANNCTTNNCGANNCTANNGGANNCTTNNGAT"
      }

      if(input$ownDesign != "") {
              bcp <- input$ownDesign
      }

      if(fileTrigger$depend() != 0) {
        if(length(input$f_name) == 1) {
          withProgress(message = 'Data processing', value = 0, {
            if(unlist(strsplit(input$f_name, split = "[.]"))[2] != "csv") {
                incProgress(1/3, detail = "Barcode extraction...")
                  tmp <- get_dat_single(bcp, givenDir, cpus)
                incProgress(2/3, detail = "Error Correction...")
                  tmp_EC <- errorCorrection(tmp, maxDist = input$maxHD, save_it = flag, cpus = cpus)
                incProgress(3/3, detail = "Finished...")
                isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "read", gDir = givenDir))))
            } else {
              incProgress(1/3, detail = paste0("Reading ", input$f_name, "..."))
                tmp <- readBCdat(path = givenDir, label = unlist(strsplit(input$f_name, split = "[.]"))[1], mask = bcp, file_name = input$f_name, s = ";")
              incProgress(2/3, detail = "Error Correction...")
                tmp_EC <- errorCorrection(tmp, maxDist = input$maxHD, save_it = FALSE, cpus = cpus)
              incProgress(3/3, detail = "Finished...")
              isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "readCsv", bcp = bcp, gDir = givenDir))))
            }
          if(sum(dim(slot(tmp, "reads"))) == 0 & fileTrigger$depend() > 0) {
            options("genBaRcode-info" = "...no barcode sequences detectable...")
            fileTrigger$reset()
          } else {
              list(
                dat = tmp,
                datEC = tmp_EC
              )
          }
         })
        } else {
          withProgress(message = 'Data processing', value = 0, {
            if(unlist(strsplit(input$f_name, split = "[.]"))[2] != "csv") {
            incProgress(1/3, detail = "Barcode extraction...")
              tmp <- get_timeS(bcp, givenDir, cpus)
            incProgress(2/3, detail = "Error Correction...")
              tmp_EC <- get_timeS_EC(tmp, givenDir, cpus)
            incProgress(3/3, detail = "Finished...")
            } else {
             tmp <- tmp_EC <- list()
             incProgress(1/3, detail = "Reading files...")
              tmp <- TMP <- list()
              for(csvFiles in 1:length(input$f_name)) {
                incProgress(1/3 + (1/3 / (length(input$f_name) * csvFiles)), detail = paste("Reading ", input$f_name[csvFiles],"...", sep = ""))
                TMP[[csvFiles]] <- readBCdat(path = givenDir, label = unlist(strsplit(input$f_name[csvFiles], split = "[.]"))[1], mask = bcp, file_name = input$f_name[csvFiles], s = ";")
              }
              tmp[[1]] <- TMP
              tmp[[2]] <- generateTimeSeriesData(TMP)
             incProgress(2/3, detail = "Error Correction...")
              TMP <- errorCorrection(TMP, maxDist = input$maxHD, save_it = FALSE, cpus = cpus)
              tmp_EC[[1]] <- TMP
              tmp_EC[[2]] <- generateTimeSeriesData(TMP)
             incProgress(3/3, detail = "Finished...")
            }
          })
          isolate(genHist$update_sCode(rbind(genHist$show_sCode(), sCode_snippets(type = "timeS", bcp = bcp, gDir = givenDir))))
          list(dat = tmp, datEC = tmp_EC)
        }
      }
    })

    ##################
    # basic functions
    ####################

    get_dat_single <- function(bcp, gDir, cpus) {

        processingRawData(file_name = input$f_name, source_dir = gDir,
                          results_dir = gDir,
                          mismatch = input$mm,
                          bc_backbone = bcp,
                          min_score = ifelse(input$quaFil, 30, 0), min_reads = input$minReads,
                          save_it = flag, seqLogo = TRUE, cpus = cpus)
    }

    get_timeS <- function(bcp, gDir, cpus) {

      dat_list <- vector("list", 2)
      dat_list[[1]] <- processingRawData(file_name = input$f_name, source_dir = gDir,
                                                                   results_dir = gDir,
                                                                   mismatch = input$mm,
                                                                   bc_backbone = bcp,
                                                                   min_score = ifelse(input$quaFil, 30, 0), min_reads = input$minReads,
                                                                   save_it = flag, cpus = cpus)
      dat_list[[2]] <- generateTimeSeriesData(dat_list[[1]])

      return(dat_list)

    }

    get_timeS_EC <- function(dat_list, gDir, cpus) {

      if(input$maxHD > 0) {
             dat_list[[1]] <- errorCorrection(dat_list[[1]], maxDist = input$maxHD,
                                                             save_it = flag,
                                                             cpus = cpus)
             dat_list[[2]] <- generateTimeSeriesData(dat_list[[1]])
      }

      return(dat_list)

    }

    sCode_snippets <- function(type, gDir = "", bcp = "", oriBCs = NULL, s = "") {

      # if(length(grep("timeS", type)) > 0) {
      #
      #     cmd <- matrix(
      #         c(
      #           paste("nameList <- c(", paste("'", input$f_name, collapse =", ", "'", sep = ""), ")"),
      #           paste("dat_list <- processingRawData(file_name = nameList, source_dir = '", gDir, "', results_dir = '", gDir,
      #           "', mismatch = ", input$mm, ", bc_backbone = '", bcp, "', min_score = ", ifelse(input$quaFil, 30, 0),", min_reads = ", input$minReads)
      #         ), ncol = 1, nrow = 2)
      #
      #     if(type == "timeS") {
      #       cmd <- rbind(cmd, matrix(c("dat <- generateTimeSeriesData(dat_list)", "plotTimeSeries(dat)", ncol = 1, nrow = 2)))
      #     }
      #     if(type == "timeS_Venn") {
      #       cmd <- rbind(cmd, "plotVennDiagramm(dat_list)")
      #     }
      #   return()
      # } else {
            return(
              switch(type,
                     timeS = matrix(
                       c(
                         paste("nameList <- c(", paste("'", input$f_name, collapse =", ", "'", sep = ""), ")"),
                         paste("dat_list <- processingRawData(file_name = nameList, source_dir = '", gDir, "', results_dir = '", gDir,
                               "', mismatch = ", input$mm, ", bc_backbone = '", bcp, "', min_score = ", ifelse(input$quaFil, 30, 0),", min_reads = ", input$minReads)
                       ), ncol = 1, nrow = 2),
                     timeS_plot = matrix(c("dat <- generateTimeSeriesData(dat_list)", "plotTimeSeries(dat)", ncol = 1, nrow = 2)),
                     timeS_Venn = "plotVennDiagramm(dat_list)",
                     read = isolate(paste("BC_dat <- processingRawData(file_name = '", input$f_name, "', source_dir = '", gDir,
                                          "', results_dir = '", gDir, "', mismatch = ", input$mm, ", bc_backbone = '", bcp,
                                          "', min_score = ", ifelse(input$quaFil, 30, 0),", min_reads = ", input$minReads, ", seqLogo = TRUE)", sep = "")),
                     readCsv = isolate(paste("BC_dat <- readBCdat(path = '", gDir, "', label = 'csvFile', mask = 'not_defined', file_name = '", input$f_name, "', s = ';')")),
                     ec = isolate(paste("BC_dat_EC <- errorCorrection(BC_dat, maxDist = ", input$maxHD, ")", sep = "")),
                     oriBCs = paste("oriBCs <- as.character(unlist(read.table('", input$file1$datapath, "', header = FALSE, sep = ", s, ", fill = TRUE)))", sep = ""),
                     kirchenplot = isolate(paste("generateKirchenplot(BC_dat, ", ifelse(identical(isolate(history$old()), input$file1$datapath), " ", "oriBCs = oriBCs,"), " loga = ", input$loga, ", col_type = ", input$palette, ")", sep = "")),
                     HD_Graph = isolate(paste("ggplotDistanceGraph(BC_dat, minDist = 1, ", ifelse(identical(isolate(history$old()), input$file1$datapath), " ", "oriBCs = oriBCs,"), " loga = ", input$loga, ", lay = '", input$graph_layout, "', complete = ", input$compl, ifelse(is.null(oriBCs()), "", paste0(", col_type = ", input$palette)), ")", sep = "")),
                     tree = isolate(paste("plotClusterTree(BC_dat, type = ", input$tree_style, ", tipLabel = ", input$show_tip_label, ")", sep = ""))
              )
            )
#      }

  }
})


