shinyUI(
  fluidPage(theme = "metaboshiny.css",
            ECharts2Shiny::loadEChartsLibrary(),
            shinyalert::useShinyalert(), # enable shinyalert (necessary for on close message)
            shinyjqui::includeJqueryUI(),
    uiOutput("currUI")
  )
)