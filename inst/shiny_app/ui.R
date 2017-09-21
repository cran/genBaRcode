

shinyUI(
  pageWithSidebar(

    headerPanel("genBaRcode Package"),

    # SideBar
    sidebarPanel(
      uiOutput("selection")
    ),

    # MainPanel
    mainPanel(
      h5(textOutput("caption"), align = "center"),

      uiOutput("G_and_T")
    )
  )
)
