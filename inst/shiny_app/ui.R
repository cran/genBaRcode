

shinyUI(

  pageWithSidebar(

    headerPanel(""),

    # SideBar
    sidebarPanel(
      uiOutput("selection"),
      uiOutput("parameters"),
      uiOutput("end")
    ),

    # MainPanel
    mainPanel(
      h5(textOutput("caption"), align = "center"),

      uiOutput("G_and_T")
    )
  )
)
